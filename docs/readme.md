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
│   ├── Sim_01/
│   │   ├── Run_00001/
│   │   └── ...
│   └── DB_01/
│       ├── Run_00001/
│       └── ...
├── G_EQDSK/                # Magnetic equilibrium data and documentation
│   ├── MAST-U/
│   ├── TCV/
│   └── WEST/
├── scripts/                # Build/run helpers and post-processing scripts
│   ├── set_ifx
│   ├── compile
│   ├── run_all.sh
│   ├── process_data.py
│   └── plot_gui.py
├── Sims/                   # GUI-generated input matrices
│   ├── Sim_01.dat
│   ├── DB_01.dat
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

A `Sim_XX.dat` or `DB_XX.dat` file contains a parameter matrix. Each row in that matrix becomes one `Run_000XX` output folder.

## GUI Workflow

From the project root:

```bash
python3 src/gui/T3ST_GUI.py
```

Use the GUI to choose a scenario, edit parameters, define manual sweeps or database selections, and export a `Sims/Sim_XX.dat` or `Sims/DB_XX.dat` file.

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
./build/T3ST
```

If this is a new terminal, load the Intel runtime first:

```bash
source scripts/set_ifx
./build/T3ST
```

When prompted:

```text
1    # Sim, or 2 for DB
1    # file number, e.g. Sim_01.dat or DB_01.dat
```

The Fortran code locates the project root automatically, so the executable can find `Sims/`, `data/`, and `G_EQDSK/` without relying on the old fixed folder-depth assumptions.

The executable also accepts non-interactive command-line arguments:

```bash
./build/T3ST 1 1
./build/T3ST 2 1
```

The first command runs `Sims/Sim_01.dat`. The second command runs `Sims/DB_01.dat`.

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

The first command runs `Sim_01.dat` through `Sim_05.dat`. The second runs `DB_01.dat` through `DB_03.dat`.

`file_type` values:

```text
1 = Sim
2 = DB
```

## Output Layout

Simulation output is written under:

```text
data/Sim_XX/Run_00001/
data/DB_XX/Run_00001/
```

Each run folder contains `parameters.dat` plus binary output files such as `Xtraj.bin`, `Ytraj.bin`, `Obs_full.bin`, and related trajectory, final-state, and transport files.

## Post-Processing Note

Some older Python and Mathematica post-processing tools still reference legacy names such as `file_01.dat` and `file_11.bin`. Current Fortran output uses semantic names like `parameters.dat`, `Xtraj.bin`, and `Obs_full.bin`, so post-processing scripts may need to be updated before use with newly generated runs.
