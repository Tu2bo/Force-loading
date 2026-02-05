This repository contains the experimental data and computational procedures associated with the paper "Force Loading on Molecular Clutches Governs the Stability of Cell Lamellipodia."

MainTime.m
Main simulation script. Runs stochastic Monte Carlo simulations and mean-field (ODE) simulations of the BR-integrated motor–clutch model, and analyzes force loading rate, spreading radius, and force per bond.

MC_BR_motor_clutch_model.m
Stochastic Monte Carlo implementation of the BR-integrated motor–clutch model.

MF_BR_motor_clutch_model.m
Mean-field (ODE) approximation of the model.
State variables: total adhesion force, number of bound clutches, and cell spreading radius.

Combined_Analysis.m
Post-processing and statistical analysis of oscillation data, including amplitude, period, and power spectral density (PSD) analysis.

All_Data_Combined.xlsx
Data file used by Combined_Analysis.m.
Must be placed in the same directory as the analysis script.

Usage

Run MainTime.m to perform simulations.
Run Combined_Analysis.m to analyze oscillation data.
