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
  - `TVSA_generate_200k_dataset.m`: MATLAB script for generating the ROM-based simulation dataset.
  - `main_c_hzn.py`: Python driver script for surrogate model training and optimization.

## How to run

### MATLAB ROM scripts

The core model is defined in `src/TVSA_CorePhysics_v2.m`.

To reproduce MATLAB-based results, run the following scripts from the repository root or ensure that the `src/` folder is added to the MATLAB path:

- `scripts/TVSA_Fig2.m`
- `scripts/TVSA_FigS4.m`
- `scripts/TVSA_generate_200k_dataset.m`

These scripts internally call the core physics model to generate cycle-resolved results.

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

- Large datasets and generated result files are not included in this repository.
- File paths in the Python scripts may need to be updated for the local computing environment before execution.
- The repository is currently maintained as part of software record `SWR-25-178`.
