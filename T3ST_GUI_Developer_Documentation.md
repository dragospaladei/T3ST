# T3ST GUI — Developer Documentation

## Architecture Overview
The application is a single-file Tkinter GUI organized into logical layers:

1. Configuration discovery & YAML parsing
2. Data model (ParamDef, ParamWidget, SweepWidgets)
3. UI construction (scrollable layout)
4. Validation & formatting
5. Export / load pipeline

## Key Classes

SimulationGUI
- Main application controller
- Owns state, widgets, and callbacks
- Coordinates validation, export, and preview

ParamDef
- Immutable parameter definition
- Derived from parameters.yaml

ParamWidget
- Binds ParamDef to a Tk widget
- Handles setting values and export validation

SweepWidgets
- UI bundle for one sweep definition

## Configuration Discovery
discover_config_dir():
- Prefers ./config/
- Falls back to script directory

require_yaml_files():
- Fails fast with user-friendly error dialogs

## Validation Strategy
- Live numeric validation during typing
- Full validation only on export
- Aggregated error reporting (no fail-fast)

## Sweep Expansion
- Sweeps parsed into (index, ndarray, name)
- Recursive Cartesian product writer
- Values formatted per parameter kind

## Export Format Guarantees
- Header always written
- Floats formatted with 8 significant digits
- Ints written as integers
- Total row count written first

## Extending the GUI
Common extension points:
- New parameter UI types via ParamDef.ui
- Additional sweep dimensions (extend sweeps list)
- Alternate exporters (CSV/HDF5)
- Headless batch mode (reuse validation + export logic)

## Known Design Decisions
- Single-file GUI for portability
- YAML schema tolerant to legacy formats
- Tkinter chosen for zero-runtime dependencies

## Non-Goals
- Editing YAML files from GUI
- Multi-row .dat editing
- Real-time simulation execution
