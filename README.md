# IUFM analysis pipeline

R script and input data for:

Lee H, Lee Y-J. Effect of the Filtered Volume Ratio on Ultrafiltration-Based Protein Binding Assessment: Experimental Characterization, Phenomenological Modeling, and Feasibility of Zero-Volume Extrapolation. *Pharmaceuticals* (full citation to be added on publication).

`IUFM_pipeline.R` computes the model-derived results of the article and the condition means of Supplementary Table S1 from the device-level measurements. It fits the filtrate volume (Eq. 4) and the membrane-associated protein deposition (Eq. 9), simulates the integrated ultrafiltration model (IUFM; Eqs. 8–12), estimates the filtrate-directed release fraction φ, determines the admissible range of the release rate constant k_rel, and extrapolates the filtrate concentration to zero volume against FVR² (Eq. 14).

## Requirements

R. The results in the article were obtained with R 4.3.3 (see `sessionInfo.txt`). Only base R is used; no packages need to be installed.

## Usage

From a terminal in the repository folder:

```
Rscript IUFM_pipeline.R
```

Or, in R or RStudio, set the working directory to the repository folder and run `source("IUFM_pipeline.R")`.

The run takes about a minute. Results are written to `iufm_output/`, which is created if necessary. The script uses no random numbers, so repeated runs give identical results.

## Input files (`iufm_input/`)

| File | Content |
|---|---|
| `in_filtrate.csv` | One row per filter device: initial, filtrate, and retentate drug concentrations (ng/mL), and the weighed filtrate and retentate (mg) |
| `in_protein.csv` | Pre-UF and retentate protein amounts of each device (mg; Section 4.2.5) |
| `in_membrane_drug.csv` | Drug extracted from each filter membrane (ng; Section 4.2.4) and the FVR of the condition (%) |
| `in_binding_R1.csv` | Protein concentrations (M), association constants (M⁻¹), numbers of binding sites, and drug molecular weights (Table 1) |
| `in_density.csv` | Densities used to convert the weighed masses to volumes (g/mL; Sections 4.2.3 and 4.2.5) |
| `in_dataset_map.csv` | Dataset codes (`hsa_cbz`, `pl_cbz`, `pl_warf`, `pl_lido`), matrices, and labels |
| `in_settings.csv` | Analysis settings (see below) |

In `in_filtrate.csv`, the columns `V_filtrate_uL` and `V_retentate_uL` contain the weighed masses in mg. The script divides them by the densities in `in_density.csv` to obtain volumes in µL.

One retentate protein amount (warfarin, 5 min, device 3) is missing and is omitted from the protein analysis.

| Setting | Value | Meaning |
|---|---|---|
| `V0_uL` | 300 | Sample volume loaded into each device (µL) |
| `k_rel` | 0.1 | Release rate constant (min⁻¹) |
| `integration_steps` | 2000 | Time steps over 200 min (step 0.1 min) |
| `analysis_window_pct` | 35 | Upper FVR limit of the zero-volume extrapolation (%) |
| `conf_level` | 0.95 | Confidence level of all intervals |

## Output files (`iufm_output/`)

| File | Content | Article |
|---|---|---|
| `out_device_values.csv` | FVR, conventional f_u, mass balance f_u, and recovery of each device (Eqs. 1–3) | Figures 1, 2, and 3A–C; Table S2 |
| `out_tableS1.csv` | Condition means and SDs | Table S1 |
| `out_parameters.csv` | F_max and k_f; σ_max and k_dep with 95% CIs; k_rel; φ with 95% CI | Table 1; Section 2.5 |
| `out_deposition_points.csv` | Membrane-associated protein fraction (mean, SD) and its Eq. 9 fit | Sections 2.3 and 2.4; Figures 3E and 4D |
| `out_simulation_summary.csv` | Simulated and measured recovery at the highest FVR, with the 95% CI of the measured mean | Section 2.5; Figure 5E |
| `out_membrane_check.csv` | Simulated membrane pool and extracted drug | Section 2.5 |
| `out_krel_range.csv` | Admissible range of k_rel | Section 4.3.5 |
| `out_zv_extrapolation.csv` | Zero-volume extrapolation of the measured data for FVR ≤ 25%, 35%, and 40% | Section 2.6; Figure 6; Note S2 |
| `out_zv_bias.csv` | Deviation of the extrapolated simulated profiles from C_ret,u(0) for k_rel = 0.014, 0.1, and 0.28 min⁻¹ | Sections 2.6 and 4.3.5; Note S2 |
| `out_baseline_shift.csv` | Vertical offset of the simulated profiles in Figure 5A–D | Section 4.3.5 |
| `sessionInfo.txt` | R session information of the run | |

## Notes

- To repeat the time-step check of Section 4.3.5, set `integration_steps` to 20000 (step 0.01 min) and compare `out_zv_bias.csv` with that of the default run. This run takes about ten times longer.
- The statistical tests and the regressions of membrane-associated protein and drug on FVR (Section 4.2.8) are not part of this script, and the script does not draw the figures.
