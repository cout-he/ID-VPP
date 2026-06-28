# DCVPP Reproduction — Zhao et al. (JMPCE 2025)

Reproducing *"Optimal Scheduling of Data Center Virtual Power Plant in
Electricity-Carbon Joint Market Under Uncertainty"* in stages.

## Phase 1 — Second-order ETP data-center load model ✅

Implements Section III, eq. (6)-(14): IT power, the second-order ETP thermal
recursion, cooling-power inversion, and the total DC load.

```
dcvpp/
  parameters.py     # ETP / IT / temperature-limit parameter sets
  dc_load_model.py  # eq. (6)-(14): P_IT, A/B recursion, Q_DC & P_cool, P_DC
scripts/
  verify_dc_load.py # stability + temperature-band verification, makes the figure
results/
  phase1_dc_load_verification.png
```

### Run
```bash
python scripts/verify_dc_load.py
```

### Parameter source
Thermal parameters migrated from Pan Dihan et al. (AEPS 2023), whose eq. (1)
ODE and eq. (5) analytical recursion are identical (up to symbol map
`Cm→Cs, Rm→Rs, Cop→k_cop`) to Zhao eq. (9) and (12).

**Decoupled calibration** (`ETPParams.calibrated`): `R` and `C` are scaled by
*different* factors because they control different things:
- **Resistance `R` sets PUE.** Steady state `Tema_ss = Temo + Ra·Q_DC`, so the
  cooling power / PUE depends on `Ra` alone (capacities cancel). `Ra,Rs` are
  Pan/1000 → full-load PUE ≈ 1.59. ✅
- **Capacity `C` sets the time constant.** With `R` fixed, `τ₁` scales linearly
  with the `C`-factor. Pan's room has `τ₁ ≈ 21 h`; scaling `C` by ×0.142
  (≈ Pan-room ×143, not ×1000) brings `τ₁` to a DC-realistic **3 h**. The
  `Cs/Ca` ratio (≈15, the solid-mass second-order structure) is preserved.

| param | value | param | value |
|---|---|---|---|
| Ca | 15.6 kWh/°C | Ra | 0.007 °C/kW |
| Cs | 235 kWh/°C | Rs | 0.0055 °C/kW |
| k_cop | 4.0 | k_IT | 1.0 |

> ⚠ Earlier draft used `Ca=110, Cs=1657` (C and R both ×1000 → `τ₁=21 h`,
> unphysically sluggish, would over-state DC demand-response flexibility in
> later phases). Recalibrated above; retarget via
> `ETPParams.calibrated(target_tau1_h=...)`.

### ⚠ Unit-of-time gotcha
With `C` in kWh/°C and `R` in °C/kW, `C·R` is in **hours**, so the recursion
exponent `b·Δt/2` needs **Δt = 1 h**, *not* 3600 s. Using 3600 underflows every
modal term to 0. The 1-hour sampling step therefore enters as `dt_h = 1.0`.

### Verification result
- Stability: modal multipliers `e^{r·Δt}` = 0.717, 0.000 → spectral radius
  0.717 < 1 ⇒ temperature provably bounded; `[A|B]` row sums = 1.
- Time constants: `τ₁=3.0 h` (slow solid mode), `τ₂=0.05 h` (fast air mode).
- Scenario A (forward recursion, T°=35°C): all in-band initial states converge
  to the 24°C fixed point; a 30°C hot start re-enters [≤27°C] in **1 h**
  (calibration anchor: ≤4 h).
- Scenario B (summer day, setpoint-tracking cooling): indoor held at 24°C,
  `P_DC` 3.3–3.8 MW (< 4 MW cap), PUE 1.51–1.68 — physically realistic.
- Scenario C (cold winter, 30% load): the inverse demands negative cooling, so
  the chiller turns **off** (`P_cool=0`) and the indoor temp drifts naturally
  (24 → 16.3 °C). Exiting the band here is *correct* model behaviour — the
  heat-balance-consistent natural drift, not a clip — and is the DR/heating
  trade-off space the Phase-5 optimiser resolves.

### Modelling decisions & caveats
See [`reproduction_notes.md`](reproduction_notes.md) for points where the paper
is ambiguous and we chose a reading, notably:
- `P_others` = 10% of **electrical** IT+cooling power (→ simultaneous solve of
  `P_cool`/`P_others`; watch the `k_cop−0.1` denominator if `k_cop` is small).
- Chiller-off uses natural recursion, never a `P_cool` clip (heat-balance
  consistency).
- Utilisation out of `[0, u_max]` (eq. 39) is **reported as a warning, never
  silently clipped** — clipping would disguise an infeasible dispatch.
- Scenario B's flat 24 °C is *imposed* (inverse model); emergent temperature
  behaviour is a Phase-5 test.

## Next phases (planned)
2. Workload model (eq. 4-5) + DR constraints (eq. 38-45).
3. Uncertainty: Gaussian copula + two-sided superquantile reserve (eq. 15-24).
4. Electricity-carbon joint-market scheduling MILP (eq. 25-49, YALMIP/Gurobi → here likely Pyomo/Gurobi).
