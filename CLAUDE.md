# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview
This is a MATLAB-based MRI simulation and data fitting project focused on inhomogeneous Magnetization Transfer (ihMT) imaging. The codebase performs steady-state simulations of MRI pulse sequences and fits experimental Z-spectrum data from various tissue samples.

## High-Level Architecture

### Core Simulation Engine
The project implements a multi-pool exchange model for MRI signal simulation:
- **Free water pool**: Bulk water protons
- **Semisolid pools**: Macromolecular components with dipolar coupling
- **Exchange modeling**: Chemical exchange between pools using Bloch-McConnell equations

Key simulation functions use matrix exponential methods to solve coupled differential equations for magnetization evolution under RF pulses.

### Data Processing Pipeline
1. **Experimental data loading**: Reads .mat files containing Z-spectrum measurements
2. **Parameter optimization**: Uses fmincon to fit tissue parameters to experimental data
3. **Visualization**: Generates comparative plots of simulated vs experimental Z-spectra

### Key Tissue Models
- Hair conditioner (hc): Test phantom with known properties
- Egg white: Protein phantom for validation
- White matter (ic): Brain tissue modeling
- BSA: Bovine serum albumin phantom

## Common Development Tasks

### Running Simulations
```matlab
% Initialize tissue parameters
tissuepars = init_tissue('hc');  % Options: 'hc', 'ic', 'egg', 'BSA'

% Run steady-state simulation
[Mz_free] = new_Dualcase_ssSPGR_ihMT_integrate(b1_band, dt, delta, tissuepars);
```

### Data Fitting
```matlab
% Main fitting scripts for different samples
Fit_hc_select_power    % Fit hair conditioner data
Fit_egg_select_power   % Fit egg white data
```

### Testing
```matlab
% Run consistency tests
datafittingconsistencytest
steadystate_confirm_tests
```

## Directory Structure

- **src/**: Core MATLAB functions for simulation and fitting
  - Integration functions: `*_integrate.m` files handle Bloch equation solving
  - Fitting functions: `Fit_*.m` files for parameter optimization
  - Utility functions: Lineshape calculations, tissue initialization
  
- **Experimental Data/**: .mat files with measured Z-spectra
- **experiments_shaihan/**: Additional analysis scripts
- **Simulation figures**: Multiple directories with .fig files for different experiments

## Key Parameters

When modifying simulations, these are the critical parameters:
- **R1/R2**: Relaxation rates for each pool
- **M0s**: Semisolid pool fraction
- **f**: Fraction of dipolar-coupled semisolid
- **k**: Exchange rate
- **T2s**: Semisolid T2 time (affects lineshape)
- **B1**: RF pulse amplitude (in μT)
- **delta**: Frequency offset (in Hz)

## Important Implementation Details

- The code supports both 1-band and 2-band RF pulse schemes
- Matrix exponential calculations use `expm()` for numerical stability
- Lineshapes can be SuperLorentzian or Gaussian
- Normalization is performed as Mz/M0 to match experimental protocols