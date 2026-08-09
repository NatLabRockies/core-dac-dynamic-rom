# core-dac-dynamic-rom
Developed at NLR for dynamic modeling, surrogate learning, and optimization of cyclic CO2 separation systems.
Physics-based dynamic reduced-order model (ROM), surrogate modeling, and optimization workflow for cyclic atmospheric CO2 separation under time-varying environmental and operating conditions.

## Repository structure

- `src/`
  - `TVSA_CorePhysics_v2.m`: core physics-based dynamic ROM.
  - `utils_c_hzn.py`: Python utility functions for surrogate model training, OMLT/Pyomo formulation, optimization, and visualization.

- `scripts/`
  - `TVSA_Fig2.m`: MATLAB script for reproducing Figure 2.
  - `TVSA_FigS4.m`: MATLAB script for reproducing Supplementary Figure S4.
  - `TVSA_make_LHS_13D_dataset.m`: MATLAB script for generating the 120,000-case cycle-resolved ROM dataset.
  - `tvsa_endpoints_ext2048_rand5000.mat`: periodic-state results used to define the feasible initial CO2 and H2O sorbent-state domain for dataset generation.
  - `main_c_hzn.py`: Python driver script for surrogate model training and optimization.

## How to run

### MATLAB ROM scripts

The core model is defined in `src/TVSA_CorePhysics_v2.m`.

To reproduce MATLAB-based results, ensure that both the src/ and scripts/ folders are added to the MATLAB path, and run the following scripts:

- `scripts/TVSA_Fig2.m`
- `scripts/TVSA_FigS4.m`
- `scripts/TVSA_make_LHS_13D_dataset.m`
  
These scripts internally call the core physics model to generate cycle-resolved results. Before running `TVSA_make_LHS_13D_dataset.m`, ensure that the `src/` folder containing `TVSA_CorePhysics_v2.m` is added to the MATLAB path. The required periodic-state file `tvsa_endpoints_ext2048_rand5000.mat` is provided in the same `scripts/` folder as the dataset-generation script. The script generates 120,000 cycle-resolved cases from 30,000 Latin-hypercube samples, each paired with four feasible initial sorbent states.

### Python surrogate and optimization workflow

The Python workflow is driven by:

- `scripts/main_c_hzn.py`

Utility functions are defined in:

- `src/utils_c_hzn.py`

The Python workflow includes:
- loading ROM-generated training data,
- formatting data for PyTorch,
- training or loading a neural-network surrogate,
- formulating an OMLT/Pyomo optimization model,
- testing optimized operation under time-varying environmental conditions,
- visualizing baseline and optimized trajectories.

## Required Python packages

The Python scripts require:

- `numpy`
- `pandas`
- `matplotlib`
- `scikit-learn`
- `torch`
- `pyomo`
- `pyomolayers`
- `omlt`

An IPOPT-compatible solver is required for the Pyomo optimization step.

## Notes

- The 120,000-case dataset is not stored directly in this repository; it can be regenerated using `scripts/TVSA_make_LHS_13D_dataset.m` together with the provided periodic-state results file.
- File paths in the Python scripts may need to be updated for the local computing environment before execution.
- The repository is currently maintained as part of software record `SWR-25-178`.
