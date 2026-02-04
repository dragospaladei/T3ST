# T3ST GUI

T3ST GUI is a Tkinter-based application for defining simulation parameters,
applying scenario presets, configuring parameter sweeps, and exporting
simulation input files.

## Requirements
- Python 3.8+
- numpy
- matplotlib
- sv-ttk
- pyyaml

## Running
python Define_Simulations.py

## Config Files
parameters.yaml — parameter definitions
scenarios.yaml — scenario presets

## Export
Files are written to Sims/Sim_XX.dat

## Parameters.yaml
Supports legacy list or enriched schema with groups, defaults, ranges, enums.

## Scenarios.yaml
Supports dict-based or legacy list-based scenarios.

## Sweeps
Up to 3 sweeps, Cartesian product expansion.
