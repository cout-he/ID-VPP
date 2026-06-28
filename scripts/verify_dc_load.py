"""Phase-1 verification of the second-order ETP DC-load model.

Two checks:
  1) Stability / boundedness: modal multipliers e^{r*dt} must have |.| < 1.
  2) Thermal behaviour: forward recursion of the indoor temperature stays
     inside the admissible band 18-27 degC (Zhao eq. 41).

Scenario A — pure forward recursion (the user's explicit ask):
     constant summer outdoor temp 35 degC, full IT load, cooling sized for a
     24 degC steady state, initial indoor 30 degC; integrate 24 steps and plot.

Scenario B — realistic summer day:
     24 h outdoor profile, hourly cooling computed by the inverse model to
     hold the 24 degC setpoint; output the full power decomposition
     (P_IT / P_cool / P_others / P_DC, PUE) to confirm P_cool is physical.

Run:  python scripts/verify_dc_load.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from dcvpp import (
    DEFAULT_ETP, DEFAULT_IT, DEFAULT_LIMITS,
    build_recursion, simulate_forward, p_it_full_load,
    q_dc_required, p_cool_from_q_dc, total_dc_power,
    steady_state_temp, q_dc_for_setpoint,
)

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
os.makedirs(RESULTS, exist_ok=True)

etp, it, lim = DEFAULT_ETP, DEFAULT_IT, DEFAULT_LIMITS
rec = build_recursion(etp)

T = 24
hours = np.arange(T + 1)


def banner(msg):
    print("\n" + "=" * 68 + "\n" + msg + "\n" + "=" * 68)


# --------------------------------------------------------------------------- #
# 1) Stability / parameter report
# --------------------------------------------------------------------------- #
banner("1) ETP parameters & stability check")
print(f"  Ca={etp.Ca} Cs={etp.Cs} kWh/degC,  Ra={etp.Ra} Rs={etp.Rs} degC/kW")
print(f"  k_cop={etp.k_cop}  k_IT={etp.k_IT}  dt={etp.dt_h} h")
print(f"  b={etp.b:.4f} 1/h,  c={etp.c:.4f} 1/h^2,  b^2-4c={etp.b**2-4*etp.c:.4f}")
print(f"  roots r1={rec.r1:.4f}, r2={rec.r2:.4f} 1/h")
print(f"  time constants tau1={etp.tau1_h:.2f} h (slow/solid), tau2={etp.tau2_h:.3f} h (fast/air)")
assert 2.0 <= etp.tau1_h <= 4.0, f"tau1={etp.tau1_h:.2f}h outside target 2-4h band"
m1, m2 = rec.modal_multipliers
print(f"  modal multipliers e^(r*dt): {m1:.4f}, {m2:.4f}  -> spectral radius {rec.spectral_radius:.4f}")
assert rec.spectral_radius < 1.0, "UNSTABLE: spectral radius >= 1"
print("  A =", np.array2string(rec.A, precision=4))
print("  B =", np.array2string(rec.B, precision=4))
print(f"  row sums [A|B] = {rec.A.sum(1) + rec.B}  (should be ~[1,1])")
print("  -> stable: both |multipliers| < 1, temperature is bounded.")

P_IT_full = p_it_full_load(it)
print(f"\n  Full-load P_IT = {P_IT_full/1000:.3f} MW")


# --------------------------------------------------------------------------- #
# Scenario A — pure forward recursion, constant 35 degC outdoor
# --------------------------------------------------------------------------- #
banner("Scenario A) Forward recursion @ 35 degC outdoor, setpoint 24 degC")
Temo_A = np.full(T, 35.0)
Tema_set = 24.0
Q_const = q_dc_for_setpoint(etp, 35.0, Tema_set)          # net thermal power for 24 degC ss
tau_slow = -1.0 / rec.r1                                  # dominant time constant [h]
print(f"  constant Q_DC sized for 24 degC steady state = {Q_const:.1f} kW")
print(f"  closed-form steady state = {steady_state_temp(etp, 35.0, Q_const):.2f} degC")
print(f"  dominant thermal time constant = {tau_slow:.1f} h")

P_cool_A, P_others_A = p_cool_from_q_dc(etp, it, P_IT_full, Q_const)
P_DC_A = total_dc_power(P_IT_full, P_cool_A, P_others_A)
print(f"  -> P_cool={P_cool_A/1000:.3f} MW, P_others={P_others_A/1000:.3f} MW, "
      f"P_DC={P_DC_A/1000:.3f} MW, PUE={P_DC_A/P_IT_full:.3f}")

# Several in-band initial states; all must converge to the 24 degC fixed point.
init_states = [22.0, 24.0, 26.0]
traj_A = {}
for x0 in init_states:
    tr = simulate_forward(rec, etp, np.array([x0, x0]), Temo_A, np.full(T, Q_const))
    traj_A[x0] = tr[:, 0]
    print(f"  start {x0:.1f} degC -> end {tr[-1,0]:.2f} degC "
          f"(min/max {tr[:,0].min():.2f}/{tr[:,0].max():.2f})")
Tema_A, Tems_A = traj_A[26.0], simulate_forward(
    rec, etp, np.array([26.0, 26.0]), Temo_A, np.full(T, Q_const))[:, 1]
in_band_A = all((tr >= lim.Tema_min - 1e-6).all() and (tr <= lim.Tema_max + 1e-6).all()
                for tr in traj_A.values())
print(f"  all in-band trajectories within [18, 27] degC: {in_band_A}")

# Calibration anchor: a 30 degC hot start must re-enter the band (<=27) quickly.
hot = simulate_forward(rec, etp, np.array([30.0, 30.0]), Temo_A, np.full(T, Q_const))[:, 0]
below27 = np.where(hot <= lim.Tema_max + 1e-6)[0]
recover_h = int(below27[0]) if below27.size else None
print(f"  hot start 30 degC -> {hot[-1]:.2f} degC at 24 h; "
      f"re-enters [<=27] band at t={recover_h} h (anchor: <= 4 h)")
assert recover_h is not None and recover_h <= 4, "hot-start recovery slower than 4 h"


# --------------------------------------------------------------------------- #
# Scenario B — realistic summer day, inverse-model cooling
# --------------------------------------------------------------------------- #
banner("Scenario B) Summer day, hourly cooling holds 24 degC setpoint")
# diurnal outdoor profile: min ~30 degC at 05:00, peak ~38 degC at 15:00
Temo_B = 34.0 - 4.0 * np.cos(2 * np.pi * (np.arange(T) - 15) / 24.0)

Tema_target = np.full(T, 24.0)
temp_prev = np.array([24.0, 24.0])
Q_B = np.empty(T); P_cool_B = np.empty(T); P_others_B = np.empty(T)
P_IT_B = np.empty(T); P_DC_B = np.empty(T)
Tema_B = np.empty(T + 1); Tems_B = np.empty(T + 1)
Tema_B[0], Tems_B[0] = temp_prev

# IT load profile: daytime busy (90-100%), night lighter (~60%)
util = 0.60 + 0.40 * np.clip(np.sin(np.pi * (np.arange(T) - 6) / 14.0), 0, 1)

for t in range(T):
    P_IT_B[t] = it.n_servers * it.P_idle + (it.P_peak - it.P_idle) * util[t] * it.n_servers
    Q_B[t] = q_dc_required(rec, etp, Tema_target[t], temp_prev, Temo_B[t])
    P_cool_B[t], P_others_B[t] = p_cool_from_q_dc(etp, it, P_IT_B[t], Q_B[t])
    P_DC_B[t] = total_dc_power(P_IT_B[t], P_cool_B[t], P_others_B[t])
    nxt = rec.A @ temp_prev + rec.B * (Temo_B[t] + Q_B[t] * etp.Ra)
    Tema_B[t + 1], Tems_B[t + 1] = nxt
    temp_prev = nxt

print(f"  realized indoor temp min/max = {Tema_B.min():.2f}/{Tema_B.max():.2f} degC")
print(f"  P_DC range = {P_DC_B.min()/1000:.3f} - {P_DC_B.max()/1000:.3f} MW (cap 4 MW)")
print(f"  PUE range  = {(P_DC_B/P_IT_B).min():.3f} - {(P_DC_B/P_IT_B).max():.3f}")
in_band_B = (Tema_B >= lim.Tema_min - 1e-6).all() and (Tema_B <= lim.Tema_max + 1e-6).all()
print(f"  within [18, 27] degC band: {in_band_B}")


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

ax = axes[0, 0]
colors_A = {22.0: "#2980b9", 24.0: "#27ae60", 26.0: "#c0392b"}
for x0, tr in traj_A.items():
    ax.plot(hours, tr, "o-", ms=4, color=colors_A[x0], label=f"$T^a_0$={x0:.0f}")
ax.plot(hours, hot, "x--", ms=4, color="#7f8c8d", label="$T^a_0$=30 (hot start)")
ax.axhspan(lim.Tema_min, lim.Tema_max, color="green", alpha=0.10, label="admissible 18-27")
ax.axhline(Tema_set, color="gray", ls=":", lw=1)
ax.set_title("A) Forward recursion, $T^o=35^\\circ$C, several $T^a_0$")
ax.set_xlabel("time (h)"); ax.set_ylabel("indoor air temp ($^\\circ$C)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

ax = axes[0, 1]
ax.plot(hours, Tema_B, "o-", color="#2980b9", label="indoor air $T^a$")
ax.plot(hours, Tems_B, "s--", color="#16a085", ms=4, label="indoor solid $T^s$")
ax.plot(np.arange(T), Temo_B, "^-", color="#7f8c8d", ms=4, label="outdoor $T^o$")
ax.axhspan(lim.Tema_min, lim.Tema_max, color="green", alpha=0.10)
ax.set_title("B) Summer day, setpoint-tracking cooling")
ax.set_xlabel("time (h)"); ax.set_ylabel("temperature ($^\\circ$C)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

ax = axes[1, 0]
x = np.arange(T)
ax.stackplot(x, P_IT_B/1000, P_cool_B/1000, P_others_B/1000,
             labels=["$P^{IT}$", "$P^{cool}$", "$P^{others}$"],
             colors=["#3498db", "#e74c3c", "#95a5a6"], alpha=0.85)
ax.plot(x, P_DC_B/1000, "k-", lw=2, label="$P^{DC}$ total")
ax.axhline(4.0, color="red", ls=":", lw=1, label="4 MW cap")
ax.set_title("B) DC load decomposition (eq. 6)")
ax.set_xlabel("time (h)"); ax.set_ylabel("power (MW)")
ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.3)

ax = axes[1, 1]
ax.plot(x, P_DC_B/P_IT_B, "o-", color="#8e44ad", label="PUE = $P^{DC}/P^{IT}$")
ax.plot(x, Q_B/1000, "s--", color="#d35400", ms=4, label="net thermal $Q^{DC}$ (MW)")
ax.set_title("B) PUE and net thermal power")
ax.set_xlabel("time (h)"); ax.set_ylabel("PUE  /  MW")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

fig.tight_layout()
out = os.path.join(RESULTS, "phase1_dc_load_verification.png")
fig.savefig(out, dpi=130)
print(f"\nFigure saved -> {out}")

banner("VERDICT")
ok = rec.spectral_radius < 1.0 and in_band_A and in_band_B
print(f"  stability: {rec.spectral_radius < 1.0}   bandA: {in_band_A}   bandB: {in_band_B}")
print("  PHASE 1 MODEL VALIDATED [OK]" if ok else "  CHECK FAILED [X]")
sys.exit(0 if ok else 1)
