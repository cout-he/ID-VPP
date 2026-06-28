"""Parameter definitions for the DCVPP data-center load model (Phase 1).

All parameters follow the notation of:
  Zhao et al., "Optimal Scheduling of Data Center Virtual Power Plant in
  Electricity-Carbon Joint Market Under Uncertainty", JMPCE (2025).
Section III, equations (6)-(14).

The second-order ETP thermal parameters are migrated from:
  Pan Dihan et al., "Identification method for air-conditioning-building
  electrothermal coupling system based on second-order equivalent thermal
  parameter analytical solution", Automation of Electric Power Systems,
  47(11):77-87, 2023.  (eq. (1) ODE == Zhao eq. (9); eq. (5) recursion ==
  Zhao eq. (12)).

Symbol map between the two papers:
    Pan  Cm, Rm   <->   Zhao  Cs, Rs
    Pan  Cop      <->   Zhao  k_cop
    Ca, Ra        identical in both papers.

IMPORTANT — unit of time
------------------------
With heat capacity C in [kWh/degC] and thermal resistance R in [degC/kW],
the product C*R is in HOURS.  Hence b = 1/(Ca*Ra)+... is in [1/h] and the
recursion exponent b*dt/2 requires dt expressed in HOURS (dt = 1.0 for the
1-hour sampling step used by Zhao et al.).  Do NOT pass dt = 3600 here:
that would drive every modal term e^{r*dt} to underflow.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# Room-scale ETP parameters identified by Pan Dihan et al. (AEPS 2023), Table 3,
# standard operating point (COP 2.9).  These are the migration baseline.
PAN_ROOM = dict(Ca=0.110, Cs=1.657, Ra=7.14, Rs=5.51)


@dataclass
class ETPParams:
    """Second-order ETP thermal parameters (Zhao eq. 9-13).

    Calibration logic (decoupled, see ETPParams.calibrated):
      * Thermal resistances R set the COOLING POWER / PUE.  Steady state is
        Tema_ss = Temo + Ra*Q_DC (eq. 9 with derivatives = 0), so Q_DC and
        thus P_cool depend on Ra ALONE -- heat capacities never enter the
        steady-state power.  Ra/Rs are kept at Pan/1000 (PUE ~ 1.59, validated).
      * Heat capacities C set the TIME CONSTANT.  With R fixed, every 1/(C*R)
        term scales by 1/f when C->f*C, so the characteristic roots scale by
        1/f and the dominant time constant tau1 scales LINEARLY by f.  We pick
        f so tau1 lands in a data-center-realistic 2-4 h band (default 3 h),
        instead of the ~21 h that a room-preserving (f=1000/1000) scaling left.
    """

    Ca: float = 15.60      # indoor-air heat capacity            [kWh/degC]
    Cs: float = 235.0      # indoor-solid heat capacity (Pan Cm) [kWh/degC]
    Ra: float = 0.007      # air <-> outdoor thermal resistance  [degC/kW]
    Rs: float = 0.0055     # solid <-> air thermal resistance (Pan Rm) [degC/kW]
    k_cop: float = 4.0     # cooling-system coefficient of performance (k^cop)
    k_IT: float = 1.0      # IT heat conversion factor (k^IT, ~1: electricity->heat)
    dt_h: float = 1.0      # recursion step in HOURS (sampling step = 1 h)

    # ---- derived characteristic-equation quantities (Zhao eq. 13) ----
    @property
    def b(self) -> float:
        return 1.0 / (self.Ca * self.Ra) + 1.0 / (self.Ca * self.Rs) + 1.0 / (self.Cs * self.Rs)

    @property
    def c(self) -> float:
        return 1.0 / (self.Ca * self.Cs * self.Ra * self.Rs)

    @property
    def roots(self) -> tuple[float, float]:
        """Characteristic roots r1 (slow), r2 (fast) of r^2+br+c=0 [1/h]."""
        disc = self.b ** 2 - 4.0 * self.c
        sq = math.sqrt(disc)
        return (-self.b + sq) / 2.0, (-self.b - sq) / 2.0

    @property
    def tau1_h(self) -> float:
        """Dominant (slow) thermal time constant -1/r1 [h]."""
        return -1.0 / self.roots[0]

    @property
    def tau2_h(self) -> float:
        """Fast thermal time constant -1/r2 [h]."""
        return -1.0 / self.roots[1]

    @classmethod
    def calibrated(cls, target_tau1_h: float = 3.0,
                   Ra: float = 0.007, Rs: float = 0.0055,
                   Ca_seed: float = 110.0, Cs_seed: float = 1657.0,
                   **kw) -> "ETPParams":
        """Build params with R fixed (PUE anchor) and C scaled so the dominant
        time constant equals target_tau1_h.  Because tau1 scales linearly with
        the C-factor when R is fixed, one analytic step is exact.

        Ca_seed/Cs_seed default to the Pan-room x1000 capacities; only their
        ratio matters (it is preserved -> second-order structure intact).
        """
        seed = cls(Ca=Ca_seed, Cs=Cs_seed, Ra=Ra, Rs=Rs, **kw)
        f = target_tau1_h / seed.tau1_h
        return cls(Ca=Ca_seed * f, Cs=Cs_seed * f, Ra=Ra, Rs=Rs, **kw)


@dataclass
class ITParams:
    """IT-equipment power parameters (Zhao eq. 7-8).

    Sized so that the full-load IT power (~2.5 MW) plus cooling/auxiliary
    lands the total DC load near the 4 MW cap reported by Zhao et al.
    (Case Studies, Section V.A).
    """

    n_servers: int = 5000      # total number of servers n
    P_idle: float = 0.30       # idle power per server  [kW]
    P_peak: float = 0.50       # peak power per server  [kW]
    u_max: float = 1.0         # max CPU utilization per server
    mu: float = 1.0            # average service rate per server [workload/h]


@dataclass
class TempLimits:
    """Indoor temperature admissible band (Zhao eq. 41-42)."""

    Tema_min: float = 18.0      # [degC]
    Tema_max: float = 27.0      # [degC]
    dTema_max: float = 3.0      # max hourly variation [degC]


# Convenience: R anchored for PUE, C calibrated for a 3 h dominant time constant.
DEFAULT_ETP = ETPParams.calibrated(target_tau1_h=3.0)
DEFAULT_IT = ITParams()
DEFAULT_LIMITS = TempLimits()
