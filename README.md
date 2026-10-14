# Hybrid Energy Storage: Physics-Based Battery Models, FMU and Energy Management

A workflow from cell physics to system design for a hybrid energy storage system (HESS):

| Step | Question | Method |
|---|---|---|
| 1. Calibration | Does the physics-based model match measured cell data? | PyBaMM SPMe calibrated against experimental 1C and 5C discharges (Ecker 2015), two-stage least squares |
| 2. Reduction and screening | Which chemistry should act as the power device, and can a fast model replace the physics? | Virtual pulse tests on PyBaMM models of a Li-ion (Chen2020) and a Na-ion (Chayambuka2022) cell; one-RC equivalent-circuit model checked against the physics on a drive-cycle profile |
| 3. EMS, sizing and FMU | How large should the supercapacitor bank be, and how should power be split? | Battery-supercapacitor HESS for frequency containment reserve (FCR); filter-based EMS with state-of-charge feedback; sizing sweep; the plant packaged as an FMI 2.0 FMU and run in a co-simulation loop |

```bash
pip install -r requirements.txt
python step1_calibration.py      # about 1 min, downloads the Ecker data through pybamm.DataLoader
python step2_reduced_models.py   # about 1 min, writes models/ecm_*.json
python step3_ems_sizing.py       # about 1 min, uses fmu/HessFmu.fmu
# rebuild the FMU after changing the plant or the models:
pythonfmu build -f fmu/HessFmu.py hess_plant.py ecm.py models --dest fmu
```

## Step 1: calibration against experimental data

Ecker et al. measured constant-current discharges of a Kokam pouch cell; PyBaMM ships the 1C and 5C curves and the matching Ecker2015 parameter set. The calibration uses each curve for what it can identify:

1. **Capacity balance from 1C** (close to equilibrium): negative electrode active material volume fraction.
2. **Transport and kinetics from 5C**: negative and positive particle diffusivities and the exchange-current density.

| Parameters | RMSE at 1C (mV) | RMSE at 5C (mV) |
|---|---:|---:|
| Literature (Ecker2015) | 36.9 | 81.5 |
| After stage 1 | 15.3 | 107.3 |
| **After stage 2** | **15.0** | **23.0** |

Stage 1 alone improves 1C but worsens 5C, which is why the second stage is fitted on high-rate data. Fitted multipliers: active material fraction 0.992, negative diffusivity 1.31, positive diffusivity 0.93, exchange-current density 0.85. With only two curves, both are used in fitting, so there is no fully independent validation set; a third rate or a pulse test would provide one.

![Calibration](figures/fig1_calibration.png)

## Step 2: reduced models and technology screening

A virtual pulse test (C/20 discharge for the open-circuit voltage; 2C pulses of 10 s with 60 s rest every 10% SOC) gives a one-RC equivalent-circuit model (ECM) per chemistry. Resistances are stored per unit capacity (Ω·Ah), so packs of any size can be built from them.

| Cell (PyBaMM parameter set) | R0 at 50% SOC (Ω·Ah) | R1 at 50% SOC (Ω·Ah) | τ at 50% SOC (s) | ECM vs physics, drive cycle (mV RMSE) |
|---|---:|---:|---:|---:|
| Li-ion, Chen2020 (SPMe) | 0.10 | 0.11 | 42 | 14.8 |
| Na-ion, Chayambuka2022 (DFN) | 0.19 | 0.24 | 128 | not used |

- **Screening result.** In these published parameter sets, the Na-ion cell has about twice the capacity-normalised resistance of the Li-ion cell, so it would not work as the power device of a Li-ion/Na-ion HESS. The power device in step 3 is therefore a supercapacitor bank. A Na-ion cell designed for power would need its own parameter set before this conclusion changes.
- **Reduced Li-ion model.** On a WLTC-shaped current profile (peak 2C), the ECM follows the SPMe within 14.8 mV RMS. The Na-ion one-RC model was not adequate on the same test, so it is used only for the screening above.

![ECM validation](figures/fig2_ecm_validation.png)

## Step 3: energy management, sizing and FMU

**Application.** A stationary 1 MW / 1 MWh Li-ion battery delivering 1 MW of FCR, using the reduced model from step 2. FCR makes the battery follow every frequency fluctuation, so it sees many fast current reversals. The frequency signal is synthetic (2 h at 1 s; slow and fast Ornstein-Uhlenbeck components with standard deviations of 20 and 6 mHz).

**Supercapacitor model.** Strings of 150 cells of 3000 F, 2.7 V and 0.29 mΩ ESR, typical datasheet values, used between half and full rated voltage.

**Energy management.** The battery follows the low-pass-filtered FCR request (time constant τ); the supercapacitor bank supplies the rest, with a feedback term that pulls its state of charge back to 50%.

| Bank (strings / usable kWh) | τ (s) | Battery current reversals per hour | Battery RMS C-rate | FCR error (% RMS) | Losses over 2 h (kWh) |
|---|---:|---:|---:|---:|---:|
| none (battery only) | | 269 | 0.097 | 0.00 | 0.79 |
| 5 / 1.7 | 10 | 50 | 0.093 | 0.00 | 1.10 |
| 5 / 1.7 | 30 | 28 | 0.091 | 0.03 | 1.32 |
| 5 / 1.7 | 120 | 19 | 0.080 | 1.31 | 1.35 |
| 20 / 6.8 | 30 | 19 | 0.083 | 0.00 | 0.83 |
| 40 / 13.7 | 120 | 14 | 0.065 | 0.00 | 0.58 |

(Full sweep in `results/step3_sizing.csv`.)

- Even the smallest bank, 1.7 kWh, removes 81 to 93% of the battery's current reversals, the micro-cycling that FCR imposes on it.
- A longer filter protects the battery more but empties a small bank: at τ = 120 s the 1.7 kWh bank runs out and FCR delivery suffers (1.3% error). Bank size and EMS tuning have to be chosen together.
- Small banks add losses (the supercapacitor works hard near half voltage), while large banks reduce them below the battery-only case. Energy efficiency, battery stress and bank cost pull in different directions, which is the trade-off a sizing study has to settle.

![EMS and sizing](figures/fig3_ems_sizing.png)

**FMU.** `fmu/HessFmu.fmu` is the battery-supercapacitor plant as an FMI 2.0 co-simulation FMU (inputs: battery and supercapacitor power set-points; outputs: delivered power, SOCs, voltages, losses). Step 3 runs the selected design (5 strings, τ = 30 s) through the FMU with fmpy as the co-simulation master, the role Simulink's FMU Import block would play, and the result matches the direct Python run exactly (`results/step3_fmu_check.csv`). The FMU is built with pythonfmu, so it needs a Python installation with numpy on the host; it has been checked with fmpy, not inside Simulink.

## Limitations

- Calibration uses two published discharge curves of one cell; no independent validation set.
- The battery is represented by its reduced ECM in step 3: no thermal model and no ageing model. Current reversals and RMS C-rate are used as stress indicators, not as a degradation prediction.
- The frequency signal is synthetic; measured grid frequency would make the sizing result specific to a market.
- The EMS is a filter with SOC feedback; an optimisation-based EMS (for example model predictive control) is the natural next step.

## Repository layout

```
step1_calibration.py      PyBaMM SPMe calibration against Ecker 2015 data
step2_reduced_models.py   pulse tests on PyBaMM models -> ECM parameters (models/*.json)
step3_ems_sizing.py       HESS EMS, sizing sweep and FMU co-simulation check
ecm.py                    ECM cell and pack, supercapacitor bank
hess_plant.py             battery + supercapacitor plant (used directly and inside the FMU)
fmu/HessFmu.py            pythonfmu wrapper; fmu/HessFmu.fmu is the built FMU
plotstyle.py              figure style
figures/, results/
```

## Data and references

- Experimental discharge curves: M. Ecker et al., J. Electrochem. Soc. 162 (2015) A1836, distributed with PyBaMM.
- Parameter sets: Chen et al., J. Electrochem. Soc. 167 (2020) 080534 (Li-ion); Chayambuka et al., Electrochim. Acta (2022), physics-based modelling of sodium-ion batteries (Na-ion); both as implemented in PyBaMM.
- PyBaMM: V. Sulzer et al., J. Open Res. Softw. 9 (2021) 14.
