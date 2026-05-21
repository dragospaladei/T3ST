# T3ST Simulation Framework

T3ST is a simulation framework for test-particle ion transport in turbulent tokamak-like magnetic equilibria. Parameter setup is handled by a Python GUI, execution is handled by the Fortran code, and output can be post-processed with Python or Mathematica tooling.

## Directory Layout

```text
Parent_Folder/
├── build/                  # Generated Fortran .o and .mod files
├── config/                 # GUI and run-generation configuration
│   ├── parameters.yaml
│   ├── scenarios.yaml
│   ├── databases.yaml
│   └── selections.yaml
├── data/                   # Simulation outputs
│   ├── Sim_001/
│   │   ├── Run_00001/
│   │   └── ...
│   └── DB_001/
│       ├── Run_00001/
│       └── ...
├── G_EQDSK/                # Magnetic equilibrium data and documentation
│   ├── MAST-U/
│   ├── TCV/
│   └── WEST/
├── scripts/                # Build/run helpers and post-processing scripts
│   ├── set_ifx
│   ├── compile
│   ├── gui
│   ├── t3st
│   ├── run_all.sh
│   ├── process_data
│   ├── Process_data.nb
│   └── Nice_photo.nb
├── Sims/                   # GUI-generated input matrices
│   ├── Sim_001.dat
│   ├── DB_001.dat
│   └── ...
└── src/
    ├── fortran/            # Fortran source and Makefile
    └── gui/                # Python GUI
```

## Configuration Model

- `config/parameters.yaml` defines the canonical parameter list, defaults, UI groups, and help text.
- `config/scenarios.yaml` defines named baseline scenarios.
- `config/databases.yaml` defines parameter spaces for larger scans.
- `config/selections.yaml` defines reusable selections from those databases, such as full, random, slice, and stride selections.

A `Sim_XXX.dat` or `DB_XXX.dat` file contains a parameter matrix. Each row in that matrix becomes one `Run_000XX` output folder.

## GUI Workflow

From the project root:

```bash
./scripts/gui
```

Use the GUI to choose a scenario, edit parameters, define manual sweeps or database selections, and export a `Sims/Sim_XXX.dat` or `Sims/DB_XXX.dat` file.

You can also start the GUI directly:

```bash
python3 src/gui/T3ST_GUI.py
```

## Build Workflow

From the project root:

```bash
source scripts/set_ifx
source scripts/compile
```

`scripts/compile` delegates to `src/fortran/Makefile`.

Build outputs are written to:

```text
build/*.o
build/*.mod
build/T3ST
```

You can also build directly with:

```bash
cd src/fortran
make
```

To clean generated build artifacts:

```bash
cd src/fortran
make clean
```

## Run One Simulation

After compiling:

```bash
./scripts/t3st
```

`scripts/t3st` loads the Intel runtime with `scripts/set_ifx`, changes to the project root, checks that `build/T3ST` exists, and then runs the executable.

In interactive mode, answer the prompts:

```text
1    # Sim, or 2 for DB
1    # file number, e.g. Sim_001.dat or DB_001.dat
```

The wrapper also accepts non-interactive command-line arguments:

```bash
./scripts/t3st 1 1
./scripts/t3st 2 1
```

The first command runs `Sims/Sim_001.dat`. The second command runs `Sims/DB_001.dat`.

You can still run the executable directly. If this is a new terminal, load the Intel runtime first:

```bash
source scripts/set_ifx
./build/T3ST 1 1
```

The Fortran code locates the project root automatically, so the executable can find `Sims/`, `data/`, and `G_EQDSK/` without relying on the old fixed folder-depth assumptions.

## Script Reference

The current helper scripts are:

```text
scripts/set_ifx       Loads Intel oneAPI compiler runtime from ~/intel/oneapi/compiler/2024.1/env/vars.sh and sets unlimited stack size.
scripts/compile       Builds the Fortran executable by running make in src/fortran.
scripts/gui           Starts the Python Tk GUI from the project root.
scripts/t3st          Loads Intel runtime and runs build/T3ST, forwarding optional arguments.
scripts/run_all.sh    Runs a range of Sim_XXX.dat or DB_XXX.dat files.
scripts/process_data  Opens scripts/Process_data.nb with wolframnb, mathematica, or xdg-open.
```

The Mathematica notebooks in `scripts/` are:

```text
scripts/Process_data.nb
scripts/Nice_photo.nb
```

`scripts/process_data` is only a launcher for `Process_data.nb`; it does not run a Python post-processing pipeline.

## Run Multiple Simulations

Use `run_all.sh`:

```bash
./scripts/run_all.sh [file_type] [first] [last]
```

Examples:

```bash
./scripts/run_all.sh 1 1 5
./scripts/run_all.sh 2 1 3
```

The first command runs `Sim_001.dat` through `Sim_005.dat`. The second runs `DB_001.dat` through `DB_003.dat`.

`file_type` values:

```text
1 = Sim
2 = DB
```

## Output Layout

Simulation output is written under:

```text
data/Sim_XXX/Run_00001/
data/DB_XXX/Run_00001/
```

Each run folder contains `parameters.dat` plus binary output files such as `Xtraj.bin`, `Ytraj.bin`, `Obs_full.bin`, and related trajectory, final-state, and transport files.

## Post-Processing Note

Some older Python and Mathematica post-processing tools still reference legacy names such as `file_01.dat` and `file_11.bin`. Current Fortran output uses semantic names like `parameters.dat`, `Xtraj.bin`, and `Obs_full.bin`, so post-processing notebooks or scripts may need to be checked before use with newly generated runs.

To open the current Mathematica post-processing notebook:

```bash
./scripts/process_data
```
