# Reproduction notes — interpretation & modelling decisions

Decisions where the source paper (Zhao et al., JMPCE 2025) is silent or
ambiguous and we had to choose. These are *our* readings, recorded so the
final report can flag them and so they can be checked against the
Supplementary Material when available.

---

## D-1. `P_others` (auxiliary load): 10% of *electrical* power, not thermal

**Paper text (Sec. III-B, after eq. 14):** *"[P_others] is assumed to be 10%
of the combined power consumption of IT equipment and cooling systems."*

**Ambiguity.** Eq. (10) is a *thermal* balance:
`Q_DC = k_IT·P_IT + P_others − k_cop·P_cool`. The cooling term enters as the
*thermal* power `k_cop·P_cool`, whereas `P_cool` itself is the *electrical*
draw of the chiller. So "combined power consumption of IT and cooling" could
mean either:
- (a) electrical: `P_others = 0.10·(P_IT + P_cool)`, or
- (b) thermal: `P_others = 0.10·(k_IT·P_IT + k_cop·P_cool)`.

**Our decision: (a), electrical.** "Power consumption" in the DC-load
literature, and in eq. (6) `P_DC = P_IT + P_cool + P_others`, denotes the
*electrical* load on the grid. `P_others` is itself an electrical term in the
total load, so taking 10% of the electrical IT+cooling draw is the consistent
reading.

**Consequence / implementation.** `P_others` and `P_cool` are mutually
dependent (`P_others` needs `P_cool`, and `P_cool` via eq. 10 needs
`P_others`). We solve the two equations simultaneously
(`p_cool_from_q_dc`), giving the closed form
`P_cool = ((k_IT + r)·P_IT − Q_DC) / (k_cop − r)`, with `r = 0.10`.
Denominator `k_cop − r = 3.9` here; **watch it if `k_cop` is ever set small.**

**To verify:** cross-check against Supplementary Material Table S-II once
obtained (it lists DC parameters and may pin down the convention).

---

## D-2. `Δt` is in **hours**, not seconds

With `C` in kWh/°C and `R` in °C/kW, the product `C·R` has units of **hours**,
so `b` (eq. 13) is in 1/h and the recursion step must be `Δt = 1 h`. Using
`Δt = 3600 s` underflows every modal term `e^{r·Δt}` to 0 and silently
produces garbage. (Confirmed by the stability/temperature anchors.)

---

## D-3. ETP parameter calibration — decoupled R / C scaling

Room-scale parameters of Pan Dihan et al. (AEPS 2023, Table 3) are migrated to
MW DC scale with **different** factors for R and C, because they govern
different quantities:
- **R sets PUE.** Steady state `Tema = Temo + Ra·Q_DC` (eq. 9, derivatives 0)
  → cooling power depends on `Ra` alone. `Ra,Rs = Pan/1000` → full-load
  PUE ≈ 1.59.
- **C sets the time constant.** With R fixed, `τ₁` scales linearly with the
  C-factor. Pan's room gives `τ₁ ≈ 21 h`; scaling C by ≈ Pan-room ×143 (i.e.
  ×0.142 vs. the earlier ×1000) brings `τ₁ → 3 h`, a DC-realistic value. The
  `Cs/Ca` ratio (≈15) is preserved → second-order structure intact.

Built by `ETPParams.calibrated(target_tau1_h=3.0)`. Retarget within 2–4 h via
the argument.

---

## D-4. Setpoint-tracking (inverse) vs. optimiser (forward) — scope of Phase 1

`q_dc_required` + `p_cool_from_q_dc` run the model *inverse*: given a desired
indoor temperature, find the cooling. Phase-1 Scenario B pins the indoor temp
at 24 °C, so its flat temperature curve is **imposed, not emergent** — it only
proves forward↔inverse self-consistency. In the paper's optimisation
(eqs. 41–42) the indoor temperature is a *decision variable* the optimiser
moves within the band to save energy. That genuine test arrives in **Phase 5**,
where temperature becomes a Gurobi variable.

**Chiller-off handling (Scenario C).** When the inverse demands `P_cool < 0`
(cold outdoor / low load: the room sheds heat fast enough to need no cooling),
we do **not** clip `P_cool` to 0 while keeping the setpoint `Q_DC` — that would
violate the eq.-(10) heat balance. Instead the chiller turns off (`P_cool = 0`)
and the step uses the *natural* net heat `Q_DC = k_IT·P_IT + P_others`; the
indoor temperature then **drifts** (Scenario C: down to ~16 °C, below the
comfort floor). This out-of-band drift is correct model behaviour — it is the
heating/DR trade-off space the Phase-5 optimiser will resolve, not a defect.

---

## D-5. Utilisation / feasibility is reported, never silently clipped

Eq. (39) `0 ≤ u ≤ u_max` is a **feasibility constraint**, not a clamp.
`p_it_from_utilisation` computes with the true `u` and emits a `RuntimeWarning`
when out of range, rather than clipping — clipping would disguise an infeasible
dispatch (one needing more active servers or a load shift) as feasible and
corrupt later optimisation results.
