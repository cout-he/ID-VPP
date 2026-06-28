"""Data-center load model based on second-order ETP theory (Phase 1).

Implements Zhao et al. (JMPCE 2025), Section III, equations (6)-(14):

    P_DC  = P_IT + P_cool + P_others                         (6)
    P_IT  = n_on*P_idle + (P_peak-P_idle)*u*n_on             (7)
    u     = E / (n_on*mu)                                    (8)
    ETP ODE / analytical recursion                           (9)-(13)
    Q_DC / P_cool inversion                                  (14)

The thermal core is the analytical recursion of Pan et al. (AEPS 2023),
eq. (5), which is provably bounded (its modal multipliers e^{r*dt} have
modulus < 1), so the simulation is numerically stable.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .parameters import ETPParams, ITParams


# --------------------------------------------------------------------------- #
# IT-equipment power  (Zhao eq. 7-8)
# --------------------------------------------------------------------------- #
def p_it_from_utilisation(it: ITParams, u: float, n_on: int | None = None) -> float:
    """IT power [kW] from CPU utilisation u in [0, u_max].  Eq. (7)."""
    n_on = it.n_servers if n_on is None else n_on
    u = float(np.clip(u, 0.0, it.u_max))
    return n_on * it.P_idle + (it.P_peak - it.P_idle) * u * n_on


def p_it_from_workload(it: ITParams, workload: float, n_on: int | None = None) -> float:
    """IT power [kW] from processed workload E (eq. 8 -> eq. 7)."""
    n_on = it.n_servers if n_on is None else n_on
    u = workload / (n_on * it.mu)
    return p_it_from_utilisation(it, u, n_on)


def p_it_full_load(it: ITParams) -> float:
    """IT power at full utilisation (all servers, u = u_max)."""
    return p_it_from_utilisation(it, it.u_max, it.n_servers)


# --------------------------------------------------------------------------- #
# Second-order ETP recursion matrices  (Zhao eq. 11-13 / Pan eq. 5)
# --------------------------------------------------------------------------- #
@dataclass
class ETPRecursion:
    """Pre-computed constant A, B matrices of the analytical recursion.

        [Tema_t; Tems_t] = A @ [Tema_{t-1}; Tems_{t-1}] + B*(Temo_t + Q_DC_t*Ra)
    """

    A: np.ndarray   # 2x2
    B: np.ndarray   # 2-vector
    r1: float
    r2: float
    M: float
    Z: float
    dt_h: float

    @property
    def modal_multipliers(self) -> tuple[float, float]:
        """Discrete-time modal multipliers e^{r*dt}; both must have |.| < 1
        for the recursion to be stable / temperatures to stay bounded
        (Pan et al. 2023, Section 1.3.2)."""
        return math.exp(self.r1 * self.dt_h), math.exp(self.r2 * self.dt_h)

    @property
    def spectral_radius(self) -> float:
        return max(abs(m) for m in self.modal_multipliers)


def build_recursion(etp: ETPParams) -> ETPRecursion:
    """Build the analytical A, B matrices (eq. 12-13)."""
    Ca, Cs, Ra, Rs, dt = etp.Ca, etp.Cs, etp.Ra, etp.Rs, etp.dt_h
    b, c = etp.b, etp.c

    disc = b * b - 4.0 * c
    if disc <= 0:
        raise ValueError(f"b^2-4c = {disc:.4g} <= 0; complex roots not supported here.")
    sq = math.sqrt(disc)
    r1 = (-b + sq) / 2.0          # dominant (slow) mode, eq. (4)/(13)
    r2 = (-b - sq) / 2.0          # fast mode

    e1 = math.exp(r1 * dt)
    e2 = math.exp(r2 * dt)
    Z = (e1 - e2) / sq            # eq. (13)
    M = (e1 + e2) / 2.0           # eq. (13)

    inv_CaRs = 1.0 / (Ca * Rs)
    inv_CsRs = 1.0 / (Cs * Rs)
    inv_CaRa = 1.0 / (Ca * Ra)

    A = np.array([
        [M - (b / 2.0) * Z + inv_CsRs * Z,            inv_CaRs * Z],
        [inv_CsRs * Z,                                 M + (b / 2.0) * Z - inv_CsRs * Z],
    ])
    B = np.array([
        1.0 - M - (b / 2.0) * Z + inv_CaRa * Z,
        1.0 - M - (b / 2.0) * Z,
    ])
    return ETPRecursion(A=A, B=B, r1=r1, r2=r2, M=M, Z=Z, dt_h=dt)


# --------------------------------------------------------------------------- #
# Forward thermal recursion  (Zhao eq. 11)
# --------------------------------------------------------------------------- #
def step_forward(rec: ETPRecursion, etp: ETPParams,
                 temp_prev: np.ndarray, Temo: float, Q_DC: float) -> np.ndarray:
    """One step of [Tema, Tems] given outdoor temp and net thermal power Q_DC."""
    return rec.A @ temp_prev + rec.B * (Temo + Q_DC * etp.Ra)


def simulate_forward(rec: ETPRecursion, etp: ETPParams,
                     temp0: np.ndarray, Temo: np.ndarray, Q_DC: np.ndarray) -> np.ndarray:
    """Forward-integrate the ETP model.

    Returns array of shape (T+1, 2): [Tema, Tems] including the initial state.
    Temo, Q_DC are length-T sequences.
    """
    T = len(Temo)
    out = np.empty((T + 1, 2))
    out[0] = temp0
    for t in range(T):
        out[t + 1] = step_forward(rec, etp, out[t], Temo[t], Q_DC[t])
    return out


# --------------------------------------------------------------------------- #
# Inverse: required thermal power & cooling power  (Zhao eq. 14)
# --------------------------------------------------------------------------- #
def q_dc_required(rec: ETPRecursion, etp: ETPParams,
                  Tema_t: float, temp_prev: np.ndarray, Temo_t: float) -> float:
    """Net thermal power Q_DC needed to reach target Tema_t this step  (eq. 14).

    Inverts eq. (11) row 1 for Q_DC:
        Tema_t = A11*Tema_{t-1} + A12*Tems_{t-1} + B11*(Temo_t + Q_DC*Ra)
    => Q_DC = ((Tema_t - A11*Tema_{t-1} - A12*Tems_{t-1})/B11 - Temo_t)/Ra
    This is algebraically identical to Zhao eq. (14) line 1.
    """
    A11, A12 = rec.A[0, 0], rec.A[0, 1]
    B11 = rec.B[0]
    rhs = (Tema_t - A11 * temp_prev[0] - A12 * temp_prev[1]) / B11
    return (rhs - Temo_t) / etp.Ra


def p_cool_from_q_dc(etp: ETPParams, it: ITParams, P_IT: float, Q_DC: float,
                     others_ratio: float = 0.10) -> tuple[float, float]:
    """Solve cooling power and auxiliary power from eq. (10) + the 10% rule.

    Eq. (10):  Q_DC = k_IT*P_IT + P_others - k_cop*P_cool
    Aux rule:  P_others = others_ratio*(P_IT + P_cool)
    Solving simultaneously:
        P_cool = ((k_IT + r)*P_IT - Q_DC) / (k_cop - r)
    Returns (P_cool, P_others), both clipped at >= 0.
    """
    r = others_ratio
    P_cool = ((etp.k_IT + r) * P_IT - Q_DC) / (etp.k_cop - r)
    P_cool = max(P_cool, 0.0)
    P_others = r * (P_IT + P_cool)
    return P_cool, P_others


def total_dc_power(P_IT: float, P_cool: float, P_others: float) -> float:
    """Total DC load, eq. (6)."""
    return P_IT + P_cool + P_others


# --------------------------------------------------------------------------- #
# Steady-state helper (sanity / sizing)
# --------------------------------------------------------------------------- #
def steady_state_temp(etp: ETPParams, Temo: float, Q_DC: float) -> float:
    """Closed-form steady-state indoor temp: Tema_ss = Temo + Ra*Q_DC.

    (Set dTema/dt = dTems/dt = 0 in eq. (9).)  Useful to size the cooling that
    holds a setpoint, and to cross-check the recursion's fixed point.
    """
    return Temo + etp.Ra * Q_DC


def q_dc_for_setpoint(etp: ETPParams, Temo: float, Tema_set: float) -> float:
    """Constant Q_DC whose steady state equals a desired setpoint."""
    return (Tema_set - Temo) / etp.Ra
