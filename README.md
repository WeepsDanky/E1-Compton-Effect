# Monte Carlo Simulation for Cs-137 Compton Scattering

This repository contains all the data, code, and notebooks used to iteratively refine and analyze Monte Carlo simulations for Cs-137 Compton scattering experiments. The project demonstrates the systematic improvement of simulations by incorporating progressively detailed detector modeling, geometry adjustments, and sensitivity analyses.

---

## Repository Structure

### 1. **Code**
- **`notebooks/compton_helpers.py`**: Python helper functions for data processing, simulation analysis, and visualization.
- **`notebooks/study.ipynb`**: Jupyter notebook for conducting and visualizing the mini-studies, including baseline, detector refinement, and geometry sensitivity analysis.

---

### 2. **Data**
#### **Simulation Data**
Located in `data/simulation/`, this directory contains:
- **Study 1**: Simplified baseline simulation data at different angles (`angle_0_deg.csv`, `angle_30_deg.csv`, etc.).
- **Study 2**: Data from simulations incorporating energy-dependent detector resolution (`simulation_angle_*deg.csv`) and residual analysis (`mini_study_2_residual.png`).
- **Study 3**: Geometry refinement data, including environmental and backscatter effects (`simulation_angle_*deg.csv`) and sensitivity metrics for lead thickness and alignment uncertainty (`*_variation.png`).

#### **Calibration Data**
Located in `data/calibration/`, this directory includes:
- Calibration data files for gamma-ray sources (e.g., Am-241, Co-57, Cs-137, Na-22).
- Spectral and experimental metadata files (`.csv`, `.txt`, `.labx`).

#### **Experimental Data**
Located in `data/experiment/`, this directory contains:
- Raw experimental spectra for Cs-137 and other sources across various configurations.
- Data for experiments with and without aluminum shielding at angles like 30°, 45°, 60°, etc.
- Spectra images and CSV files (`data/experiment/*`).

---

## Features

1. **Monte Carlo Simulation**:
   - Custom Python implementation of Monte Carlo photon tracking.
   - Iterative refinement for detector resolution, backscatter, and environmental effects.

2. **Sensitivity Analysis**:
   - Impact of lead shielding and angular alignment uncertainties on the simulation.

3. **Calibration & Experimental Validation**:
   - Real-world calibration data for validating simulation results.
   - Experimentally measured Compton spectra for comparison.
