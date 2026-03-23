"""
Define_Simulations.py (extended, refactored)

Adds database-driven run generation on top of the existing manual sweep workflow.

New concepts:
- databases.yaml  : defines large parameter spaces ("run databases")
- selections.yaml : defines reusable subsets of those databases

Supported selection methods:
- full
- random
- slice
- stride

Main idea:
- Manual mode: current GUI values + up to 3 sweeps
- Database mode: choose a database + selection, then export selected runs

Dependencies: tkinter, numpy, matplotlib, sv_ttk, pyyaml
"""

from __future__ import annotations

import itertools
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import tkinter.font as tkFont

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

import sv_ttk
import yaml


# -----------------------------
# UI layout constants
# -----------------------------

UI_TOTAL_COLUMNS = 10
UI_PARAM_COLS_PER_ROW = 6
UI_ENTRY_WIDTH = 8
UI_COMBO_WIDTH = 5
TOOLTIP_FONT_FAMILY = "tahoma"
TOOLTIP_FONT_SIZE = 11
TOOLTIP_WRAP_PX = 420

MAX_PREVIEW_ROWS = 20
MAX_EXPORT_ROWS = 100_000


# -----------------------------
# Config discovery / YAML loading
# -----------------------------

def discover_config_dir() -> Path:
    script_dir = Path(__file__).resolve().parent
    cfg = script_dir / "config"
    return cfg if cfg.is_dir() else script_dir


CONFIG_DIR = discover_config_dir()


def load_yaml(path: Path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def require_yaml_files(config_dir: Path, filenames: List[str]) -> Dict[str, Path]:
    missing, resolved = [], {}
    for name in filenames:
        p = config_dir / name
        if not p.exists():
            missing.append(str(p))
        else:
            resolved[name] = p
    if missing:
        raise FileNotFoundError("Missing required config file(s):\n\n" + "\n".join(missing))
    return resolved


def format_sweep_enum_value(x: float) -> int:
    if not np.isfinite(x):
        raise ValueError("not finite")
    xi = int(round(x))
    if abs(x - xi) > 1e-12:
        raise ValueError(f"expected integer-like value, got {x}")
    return xi


# -----------------------------
# Tooltip
# -----------------------------

class ToolTip:
    def __init__(self, widget, text: str):
        self.widget = widget
        self.text = text
        self.tip = None
        self._after_id = None
        self.widget.bind("<Enter>", self._schedule_show)
        self.widget.bind("<Leave>", self.hide)

    def _schedule_show(self, event=None):
        self._cancel_scheduled()
        self._after_id = self.widget.after(250, self.show)

    def _cancel_scheduled(self):
        if self._after_id is not None:
            try:
                self.widget.after_cancel(self._after_id)
            except Exception:
                pass
            self._after_id = None

    def show(self, event=None):
        if self.tip:
            return
        x = self.widget.winfo_pointerx() + 15
        y = self.widget.winfo_pointery() + 10
        self.tip = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        label = tk.Label(
            tw, text=self.text, justify="left", anchor="w",
            wraplength=TOOLTIP_WRAP_PX, background="#ffffe0",
            relief="solid", borderwidth=1,
            font=(TOOLTIP_FONT_FAMILY, TOOLTIP_FONT_SIZE, "normal"),
        )
        label.pack(ipadx=6, ipady=4)
        tw.update_idletasks()
        w, h = tw.winfo_reqwidth(), tw.winfo_reqheight()
        scr_w, scr_h = tw.winfo_screenwidth(), tw.winfo_screenheight()
        tw.wm_geometry(f"+{min(x, scr_w - w - 10)}+{min(y, scr_h - h - 10)}")

    def hide(self, event=None):
        self._cancel_scheduled()
        if self.tip:
            self.tip.destroy()
            self.tip = None


# -----------------------------
# Core data models
# -----------------------------

@dataclass(frozen=True)
class Choice:
    value: int
    label: str


@dataclass
class ParamDef:
    name: str
    kind: str  # "float" | "int" | "enum" | "string"
    group: str = "Ungrouped"
    units: str = ""
    help: str = ""
    default: Optional[Any] = None
    min: Optional[float] = None
    max: Optional[float] = None
    ui: str = ""
    choices: Optional[List[Choice]] = None


@dataclass
class ParamWidget:
    pdef: ParamDef
    label: ttk.Label
    widget: tk.Widget
    var: Optional[tk.Variable] = None
    enum_label_to_value: Optional[Dict[str, int]] = None
    enum_value_to_label: Optional[Dict[int, str]] = None

    def set_value(self, value: Any):
        if value is None:
            value = ""

        if self.pdef.kind == "enum":
            label = ""
            try:
                vi = int(float(value))
                if self.enum_value_to_label and vi in self.enum_value_to_label:
                    label = self.enum_value_to_label[vi]
            except Exception:
                vstr = str(value)
                if self.enum_label_to_value and vstr in self.enum_label_to_value:
                    label = vstr
            if isinstance(self.widget, ttk.Combobox):
                self.widget.set(label)
            return

        if self.pdef.ui == "checkbox" and isinstance(self.var, tk.IntVar):
            try:
                self.var.set(int(float(value)))
            except Exception:
                self.var.set(0)
            return

        if isinstance(self.widget, ttk.Entry):
            self.widget.config(state=tk.NORMAL)
            self.widget.delete(0, tk.END)
            self.widget.insert(0, str(value))

    def get_export_value(self) -> Tuple[Optional[str], Optional[str]]:
        name = self.pdef.name

        if self.pdef.kind == "enum":
            if not isinstance(self.widget, ttk.Combobox):
                return None, f"Parameter {name}: internal enum widget error"
            label = self.widget.get().strip()
            if label == "":
                return None, f"Parameter {name}: empty"
            if not self.enum_label_to_value or label not in self.enum_label_to_value:
                return None, f"Parameter {name}: invalid choice {label!r}"
            return str(int(self.enum_label_to_value[label])), None

        if self.pdef.ui == "checkbox" and isinstance(self.var, tk.IntVar):
            v = int(self.var.get())
            err = validate_range(float(v), self.pdef.min, self.pdef.max)
            if err:
                return None, f"Parameter {name}: {err}"
            return str(v), None

        if not isinstance(self.widget, ttk.Entry):
            return None, f"Parameter {name}: unsupported widget"

        raw = self.widget.get().strip()
        if raw == "":
            return None, f"Parameter {name}: empty"
        if self.pdef.kind in ("int", "float"):
            if not is_float_like(raw):
                return None, f"Parameter {name}: not numeric ({raw!r})"
            try:
                x = float(raw)
                if not np.isfinite(x):
                    return None, f"Parameter {name}: not finite"
                rng_err = validate_range(x, self.pdef.min, self.pdef.max)
                if rng_err:
                    return None, f"Parameter {name}: {rng_err}"
                return format_value(raw, self.pdef.kind), None
            except Exception as e:
                return None, f"Parameter {name}: {e}"

        return raw, None


@dataclass
class SweepWidgets:
    label: ttk.Label
    param: ttk.Combobox
    points: ttk.Entry
    init: ttk.Entry
    final: ttk.Entry

    def attach_tooltips(self):
        ToolTip(self.param, "Parameter name to sweep")
        ToolTip(self.points, "Number of points in the sweep (integer ≥ 1)")
        ToolTip(self.init, "Initial value for this parameter")
        ToolTip(self.final, "Final value for this parameter")


@dataclass(frozen=True)
class DatabaseAxis:
    name: str
    mode: str
    min: float
    max: float
    nvals: int


@dataclass
class DatabaseDef:
    name: str
    description: str
    base_scenario: str
    varying_parameters: List[DatabaseAxis]

    def axis_map(self) -> Dict[str, DatabaseAxis]:
        return {a.name: a for a in self.varying_parameters}


@dataclass
class SelectionDef:
    name: str
    description: str
    database: str
    method: str
    n_samples: Optional[int] = None
    seed: Optional[int] = None
    fixed_indices: Optional[Dict[str, int]] = None
    fixed_values: Optional[Dict[str, float]] = None
    free_parameters: Optional[List[str]] = None
    strides: Optional[Dict[str, int]] = None


# -----------------------------
# Parsing helpers
# -----------------------------

def _as_float_or_none(x) -> Optional[float]:
    if x is None or x == "":
        return None
    try:
        return float(x)
    except Exception:
        return None


def _require_field(d: dict, key: str, context: str) -> str:
    """Extract and validate a required non-empty string field from a dict."""
    val = str(d.get(key, "")).strip()
    if not val:
        raise ValueError(f"{context}: missing required field '{key}'")
    return val


def parse_parameters_yaml(obj) -> List[ParamDef]:
    if obj is None:
        return []

    # Legacy: plain list of names
    if isinstance(obj, list) and all(isinstance(x, str) for x in obj):
        return [ParamDef(name=s, kind="float") for s in obj]

    if isinstance(obj, dict) and isinstance(obj.get("parameters"), list):
        out: List[ParamDef] = []
        for i, p in enumerate(obj["parameters"]):
            if not isinstance(p, dict) or "name" not in p:
                raise ValueError(f"parameters.yaml: item #{i} must be a dict with a 'name' field")
            name = str(p["name"]).strip()
            kind = str(p.get("kind", "float")).strip()

            choices: Optional[List[Choice]] = None
            if kind == "enum":
                choices_raw = p.get("choices")
                if not isinstance(choices_raw, list) or not choices_raw:
                    raise ValueError(f"parameters.yaml: enum param {name!r} must provide non-empty 'choices' list")
                for ch in choices_raw:
                    if not isinstance(ch, dict) or "value" not in ch or "label" not in ch:
                        raise ValueError(
                            f"parameters.yaml: enum param {name!r} choices must be list of {{value, label}} dicts"
                        )
                choices = [Choice(value=int(ch["value"]), label=str(ch["label"])) for ch in choices_raw]

            out.append(ParamDef(
                name=name,
                kind=kind,
                group=str(p.get("group", "Ungrouped") or "Ungrouped"),
                units=str(p.get("units", "") or ""),
                help=str(p.get("help", "") or ""),
                default=p.get("default"),
                min=_as_float_or_none(p.get("min")),
                max=_as_float_or_none(p.get("max")),
                ui=str(p.get("ui", "") or ""),
                choices=choices,
            ))
        return out

    raise ValueError("parameters.yaml: unsupported format (expected legacy list or enriched dict with 'parameters')")


def parse_databases_yaml(obj) -> List[DatabaseDef]:
    if not isinstance(obj, dict) or not isinstance(obj.get("databases"), list):
        raise ValueError("databases.yaml: expected a dict with key 'databases' as a list")

    out: List[DatabaseDef] = []
    for i, db in enumerate(obj["databases"]):
        if not isinstance(db, dict):
            raise ValueError(f"databases.yaml: database item #{i} must be a dict")
        ctx = f"databases.yaml[{i}]"
        name = _require_field(db, "name", ctx)
        base_scenario = _require_field(db, "base_scenario", f"databases.yaml: database {name!r}")

        varying_raw = db.get("varying_parameters")
        if not isinstance(varying_raw, list) or not varying_raw:
            raise ValueError(f"databases.yaml: database {name!r} must define non-empty varying_parameters")

        axes: List[DatabaseAxis] = []
        for j, a in enumerate(varying_raw):
            if not isinstance(a, dict):
                raise ValueError(f"databases.yaml: database {name!r} axis #{j} must be a dict")
            aname = _require_field(a, "name", f"databases.yaml: database {name!r} axis #{j}")
            mode = str(a.get("mode", "")).strip()
            if mode != "linspace":
                raise ValueError(
                    f"databases.yaml: database {name!r} axis {aname!r}: only mode='linspace' is supported"
                )
            try:
                axes.append(DatabaseAxis(
                    name=aname, mode=mode,
                    min=float(a["min"]), max=float(a["max"]), nvals=int(a["nvals"])
                ))
            except Exception:
                raise ValueError(
                    f"databases.yaml: database {name!r} axis {aname!r} requires numeric min/max and integer nvals"
                )
            if axes[-1].nvals < 1:
                raise ValueError(f"databases.yaml: database {name!r} axis {aname!r}: nvals must be >= 1")

        out.append(DatabaseDef(name=name, description=str(db.get("description", "") or ""),
                               base_scenario=base_scenario, varying_parameters=axes))
    return out


def parse_selections_yaml(obj) -> List[SelectionDef]:
    if not isinstance(obj, dict) or not isinstance(obj.get("selections"), list):
        raise ValueError("selections.yaml: expected a dict with key 'selections' as a list")

    VALID_METHODS = ("full", "random", "slice", "stride")
    out: List[SelectionDef] = []
    for i, s in enumerate(obj["selections"]):
        if not isinstance(s, dict):
            raise ValueError(f"selections.yaml: selection item #{i} must be a dict")
        name = _require_field(s, "name", f"selections.yaml[{i}]")
        database = _require_field(s, "database", f"selections.yaml: selection {name!r}")
        method = str(s.get("method", "")).strip()
        if method not in VALID_METHODS:
            raise ValueError(
                f"selections.yaml: selection {name!r}: unsupported method {method!r} "
                f"(supported: {', '.join(VALID_METHODS)})"
            )

        def _opt_int(key):
            v = s.get(key)
            if v is None:
                return None
            try:
                return int(v)
            except Exception:
                raise ValueError(f"selections.yaml: selection {name!r}: {key} must be an integer")

        def _opt_dict_int(key):
            v = s.get(key)
            if v is None:
                return None
            if not isinstance(v, dict):
                raise ValueError(f"selections.yaml: selection {name!r}: {key} must be a dict")
            return {str(k): int(val) for k, val in v.items()}

        def _opt_dict_float(key):
            v = s.get(key)
            if v is None:
                return None
            if not isinstance(v, dict):
                raise ValueError(f"selections.yaml: selection {name!r}: {key} must be a dict")
            return {str(k): float(val) for k, val in v.items()}

        def _opt_list_str(key):
            v = s.get(key)
            if v is None:
                return None
            if not isinstance(v, list):
                raise ValueError(f"selections.yaml: selection {name!r}: {key} must be a list")
            return [str(x) for x in v]

        out.append(SelectionDef(
            name=name,
            description=str(s.get("description", "") or ""),
            database=database,
            method=method,
            n_samples=_opt_int("n_samples"),
            seed=_opt_int("seed"),
            fixed_indices=_opt_dict_int("fixed_indices"),
            fixed_values=_opt_dict_float("fixed_values"),
            free_parameters=_opt_list_str("free_parameters"),
            strides=_opt_dict_int("strides"),
        ))

    return out


def is_float_like(text: str) -> bool:
    try:
        float(text)
        return text.strip() != ""
    except Exception:
        return False


def format_value(text: str, kind: str) -> str:
    if kind == "int":
        x = float(text)
        if not np.isfinite(x):
            raise ValueError("not finite")
        xi = int(round(x))
        if abs(x - xi) > 1e-12:
            raise ValueError(f"expected integer-like value, got {text!r}")
        return str(xi)
    x = float(text)
    if not np.isfinite(x):
        raise ValueError("not finite")
    return f"{x:.8g}"


def validate_range(x: float, pmin: Optional[float], pmax: Optional[float]) -> Optional[str]:
    if pmin is not None and x < pmin:
        return f"must be ≥ {pmin}"
    if pmax is not None and x > pmax:
        return f"must be ≤ {pmax}"
    return None


# -----------------------------
# Database / selection logic
# -----------------------------

def build_axis_values(axis: DatabaseAxis) -> np.ndarray:
    if axis.mode == "linspace":
        return np.linspace(axis.min, axis.max, axis.nvals)
    if axis.mode == "logspace":
        if axis.min <= 0 or axis.max <= 0:
            raise ValueError(f"logspace axis {axis.name!r} requires min > 0 and max > 0")
        return np.logspace(np.log10(axis.min), np.log10(axis.max), axis.nvals)
    raise ValueError(f"Unsupported axis mode: {axis.mode!r}")


def nearest_axis_index(axis: DatabaseAxis, target: float) -> int:
    return int(np.argmin(np.abs(build_axis_values(axis) - target)))


def resolved_fixed_indices(
    db: DatabaseDef, sel: SelectionDef,
    scenario_map: Dict[str, Dict[str, Any]], parameter_defs: List[ParamDef],
) -> Dict[str, int]:
    out: Dict[str, int] = {}
    amap = db.axis_map()
    defaults = {p.name: p.default for p in parameter_defs}

    if sel.free_parameters:
        free = set(sel.free_parameters)
        base = scenario_map.get(db.base_scenario, {}) or {}
        for pname, axis in amap.items():
            if pname not in free:
                base_value = base.get(pname, defaults.get(pname))
                if base_value is not None:
                    out[pname] = nearest_axis_index(axis, float(base_value))

    if sel.fixed_indices:
        out.update(sel.fixed_indices)

    if sel.fixed_values:
        for pname, target in sel.fixed_values.items():
            if pname in amap:
                out[pname] = nearest_axis_index(amap[pname], float(target))

    return out


def stride_indices(nvals: int, step: int) -> List[int]:
    if step < 1:
        raise ValueError("stride must be >= 1")
    idx = list(range(0, nvals, step))
    if idx[-1] != nvals - 1:
        idx.append(nvals - 1)
    return idx


def database_axes_in_order(db: DatabaseDef) -> List[DatabaseAxis]:
    return list(db.varying_parameters)


def database_total_runs(db: DatabaseDef) -> int:
    total = 1
    for axis in db.varying_parameters:
        total *= axis.nvals
    return total


def unflatten_index(run_id: int, dims: List[int]) -> Tuple[int, ...]:
    if run_id < 0:
        raise ValueError("run_id must be >= 0")
    coords, rem = [], run_id
    for d in dims:
        coords.append(rem % d)
        rem //= d
    if rem != 0:
        raise ValueError("run_id out of range for these dimensions")
    return tuple(coords)


def _selection_ranges(
    db: DatabaseDef, sel: SelectionDef,
    scenario_map: Dict[str, Dict[str, Any]], parameter_defs: List[ParamDef],
) -> List:
    """
    Return a list of index iterables (one per axis) that define the selected runs.
    Used by both selection_count and selection_index_iterator.
    Returns None for random selections (handled separately).
    """
    axes = database_axes_in_order(db)

    if sel.method == "full":
        return [range(a.nvals) for a in axes]

    if sel.method == "slice":
        fixed = resolved_fixed_indices(db, sel, scenario_map, parameter_defs)
        ranges = []
        for axis in axes:
            if axis.name in fixed:
                idx = int(fixed[axis.name])
                if idx < 0 or idx >= axis.nvals:
                    raise ValueError(
                        f"Selection {sel.name!r}: fixed index for {axis.name!r} must be in [0, {axis.nvals - 1}]"
                    )
                ranges.append([idx])
            else:
                ranges.append(range(axis.nvals))
        return ranges

    if sel.method == "stride":
        strides = sel.strides or {}
        return [stride_indices(a.nvals, int(strides.get(a.name, 1))) for a in axes]

    return None  # random — caller handles it


def selection_count(
    db: DatabaseDef, sel: SelectionDef,
    scenario_map: Dict[str, Dict[str, Any]], parameter_defs: List[ParamDef],
) -> int:
    if sel.method == "random":
        if sel.n_samples is None or sel.n_samples < 1:
            raise ValueError(f"Selection {sel.name!r}: random requires n_samples >= 1")
        total = database_total_runs(db)
        if sel.n_samples > total:
            raise ValueError(
                f"Selection {sel.name!r}: n_samples={sel.n_samples} exceeds total database size {total}"
            )
        return sel.n_samples

    ranges = _selection_ranges(db, sel, scenario_map, parameter_defs)
    if ranges is None:
        raise ValueError(f"Unsupported selection method: {sel.method}")
    total = 1
    for r in ranges:
        total *= len(r) if not isinstance(r, range) else len(r)
    return total


def selection_index_iterator(
    db: DatabaseDef, sel: SelectionDef,
    scenario_map: Dict[str, Dict[str, Any]], parameter_defs: List[ParamDef],
) -> Iterable[Tuple[int, ...]]:
    if sel.method == "random":
        if sel.n_samples is None or sel.n_samples < 1:
            raise ValueError(f"Selection {sel.name!r}: random requires n_samples >= 1")
        total = database_total_runs(db)
        if sel.n_samples > total:
            raise ValueError(
                f"Selection {sel.name!r}: n_samples={sel.n_samples} exceeds total database size {total}"
            )
        dims = [a.nvals for a in database_axes_in_order(db)]
        rng = np.random.default_rng(sel.seed)
        picked = rng.choice(total, size=sel.n_samples, replace=False)
        return (unflatten_index(int(rid), dims) for rid in picked.tolist())

    ranges = _selection_ranges(db, sel, scenario_map, parameter_defs)
    if ranges is None:
        raise ValueError(f"Unsupported selection method: {sel.method}")
    return itertools.product(*ranges)


# -----------------------------
# Main GUI
# -----------------------------

class SimulationGUI:
    CONFIG_FILE = "gui_config.json"
    REQUIRED_YAMLS = ["parameters.yaml", "scenarios.yaml", "databases.yaml", "selections.yaml"]

    def __init__(self, root: tk.Tk):
        self.root = root
        self.dark_mode = False
        self._frame_window_id = None

        self.load_config()

        self.root.title("T3ST: GUI Studio")
        screen_w = self.root.winfo_screenwidth()
        screen_h = self.root.winfo_screenheight()
        w, h = int(screen_w * 0.84), int(screen_h * 0.88)
        self.root.geometry(f"{w}x{h}+{(screen_w - w) // 2}+{(screen_h - h) // 2}")

        sv_ttk.set_theme("dark" if self.dark_mode else "light")

        try:
            paths = require_yaml_files(CONFIG_DIR, self.REQUIRED_YAMLS)
            params_raw     = load_yaml(paths["parameters.yaml"]) or {}
            scenarios_raw  = load_yaml(paths["scenarios.yaml"]) or {}
            databases_raw  = load_yaml(paths["databases.yaml"]) or {}
            selections_raw = load_yaml(paths["selections.yaml"]) or {}
        except Exception as e:
            messagebox.showerror("Startup error", str(e))
            raise

        self.params: List[ParamDef] = self._load_parsed("parameters.yaml", parse_parameters_yaml, params_raw)
        if not self.params:
            messagebox.showerror("Config error", "parameters.yaml contains no parameters.")
            raise ValueError("Empty parameters.yaml")

        self.labels: List[str] = [p.name for p in self.params]
        self.param_by_name: Dict[str, ParamDef] = {p.name: p for p in self.params}

        self.scenarios: Dict[str, object] = {str(k): v for k, v in (scenarios_raw or {}).items()}
        self.valid_scenarios = list(self.scenarios.keys()) or ["(none)"]

        self.databases: List[DatabaseDef] = self._load_parsed("databases.yaml", parse_databases_yaml, databases_raw)
        self.selections: List[SelectionDef] = self._load_parsed("selections.yaml", parse_selections_yaml, selections_raw)

        self.database_by_name: Dict[str, DatabaseDef] = {d.name: d for d in self.databases}
        self.selection_by_name: Dict[str, SelectionDef] = {s.name: s for s in self.selections}

        warnings = self._collect_config_warnings()

        self.param_widgets: List[ParamWidget] = []
        self.param_rows: List[Tuple[ttk.Label, tk.Widget]] = []
        self.sweeps: List[SweepWidgets] = []

        self.mode_var = tk.StringVar(value="manual")
        self.database_var = tk.StringVar()
        self.selection_var = tk.StringVar()
        self.db_scenario_var = tk.StringVar()  # scenario override for database mode

        # Override widgets for database mode: {param_name: ParamWidget}
        self._db_override_widgets: Dict[str, ParamWidget] = {}
        self._db_override_frame: Optional[ttk.Frame] = None

        self.setup_fonts()
        self.vcmd = (self.root.register(self.validate_numeric), "%P", "%W")

        # Undo / redo stacks — each entry is a snapshot dict {name: raw_value}
        self._undo_stack: List[Dict[str, Any]] = []
        self._redo_stack: List[Dict[str, Any]] = []

        self.create_widgets()

        self.root.after(50, self.on_scenario_change)
        self.root.after(50, self.initialize_database_controls)

        if warnings:
            messagebox.showwarning(
                "Configuration warnings",
                "\n".join(warnings[:40]) + ("\n\n… (more)" if len(warnings) > 40 else "")
            )

    def _load_parsed(self, filename: str, parser, raw):
        """Parse a YAML object, showing a config error dialog on failure."""
        try:
            return parser(raw)
        except Exception as e:
            messagebox.showerror("Config error", f"Invalid {filename}:\n\n{e}")
            raise

    def _collect_config_warnings(self) -> List[str]:
        """Run all cross-reference validations and return combined warning strings."""
        warnings: List[str] = []
        labels_set = set(self.labels)
        scenario_set = set(self.scenarios.keys())

        # Scenario warnings
        for scen_name, values in self.scenarios.items():
            if isinstance(values, dict):
                extras = [k for k in values if k not in labels_set]
                if extras:
                    warnings.append(f"Scenario {scen_name!r}: unknown parameter key(s): {extras}")
            elif isinstance(values, list):
                if len(values) != len(self.labels):
                    warnings.append(
                        f"Scenario {scen_name!r}: list has {len(values)} values, expected {len(self.labels)}"
                    )
            else:
                warnings.append(f"Scenario {scen_name!r}: unsupported type {type(values).__name__}")

        # Database warnings
        for db in self.databases:
            if db.base_scenario not in scenario_set:
                warnings.append(
                    f"Database {db.name!r}: base_scenario {db.base_scenario!r} not found in scenarios.yaml"
                )
            axis_names = [a.name for a in db.varying_parameters]
            if len(axis_names) != len(set(axis_names)):
                warnings.append(f"Database {db.name!r}: duplicate varying parameter names")
            for name in axis_names:
                if name not in labels_set:
                    warnings.append(f"Database {db.name!r}: varying parameter {name!r} not found in parameters.yaml")

        # Selection warnings
        for sel in self.selections:
            if sel.database not in self.database_by_name:
                warnings.append(f"Selection {sel.name!r}: database {sel.database!r} not found in databases.yaml")
                continue
            db = self.database_by_name[sel.database]
            axis_names = {a.name for a in db.varying_parameters}
            axis_map = db.axis_map()

            if sel.method == "random" and (sel.n_samples is None or sel.n_samples < 1):
                warnings.append(f"Selection {sel.name!r}: random selection requires n_samples >= 1")

            for key, idx_map in [("fixed_indices", sel.fixed_indices), ("fixed_values", sel.fixed_values),
                                  ("free_parameters", {p: None for p in (sel.free_parameters or [])}),
                                  ("strides", sel.strides)]:
                if not idx_map:
                    continue
                for pname in idx_map:
                    if pname not in axis_names:
                        warnings.append(f"Selection {sel.name!r}: {key} refers to unknown axis {pname!r}")
                    elif key == "fixed_indices":
                        idx = int(idx_map[pname])
                        if idx < 0 or idx >= axis_map[pname].nvals:
                            warnings.append(f"Selection {sel.name!r}: fixed index for {pname!r} out of range")
                    elif key == "strides" and int(idx_map[pname]) < 1:
                        warnings.append(f"Selection {sel.name!r}: stride for {pname!r} must be >= 1")

        return warnings

    # -----------------------------
    # Config
    # -----------------------------

    def load_config(self):
        cfg_path = Path(self.CONFIG_FILE)
        if cfg_path.exists():
            try:
                cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
                self.dark_mode = bool(cfg.get("dark_mode", False))
            except Exception:
                self.dark_mode = False

    def save_config(self):
        Path(self.CONFIG_FILE).write_text(
            json.dumps({"dark_mode": self.dark_mode}), encoding="utf-8"
        )

    # -----------------------------
    # UI setup
    # -----------------------------

    def setup_fonts(self):
        tkFont.nametofont("TkDefaultFont").configure(size=10, weight="bold")

    def configure_numeric_entry_styles(self):
        style = ttk.Style()
        fg_normal = "white" if self.dark_mode else "black"
        style.configure("Numeric.TEntry", foreground=fg_normal)
        style.configure("NumericInvalid.TEntry", foreground="red")
        style.map("TButton", background=[("active", "#4b9ce2")], foreground=[("active", "white")])
        style.configure("Section.TLabel", font=("Helvetica", 12, "bold"))
        style.configure("Param.TLabel", font=("Consolas", 10, "bold"))
        style.configure("GroupHeader.TLabel", font=("Helvetica", 11, "bold"))

    def create_widgets(self):
        self.root.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)

        self.canvas = tk.Canvas(self.root, highlightthickness=0)
        self.scrollbar = ttk.Scrollbar(self.root, orient="vertical", command=self.canvas.yview)
        self.canvas.grid(row=0, column=0, sticky="nsew")
        self.scrollbar.grid(row=0, column=1, sticky="ns")

        self.scrollable_frame = ttk.Frame(self.canvas)
        self.scrollable_frame.bind(
            "<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )
        self._frame_window_id = self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfig(self._frame_window_id, width=e.width))
        self.canvas.configure(yscrollcommand=self.scrollbar.set)
        self.canvas.bind_all("<MouseWheel>", self._on_mousewheel)
        self.canvas.bind_all("<Button-4>", lambda e: self.canvas.yview_scroll(-1, "units"))
        self.canvas.bind_all("<Button-5>", lambda e: self.canvas.yview_scroll(1, "units"))

        self.configure_numeric_entry_styles()
        self.build_layout()

        # Status bar — sits below the canvas, always visible
        self._status_var = tk.StringVar(value="Ready.")
        status_bar = ttk.Label(
            self.root, textvariable=self._status_var,
            anchor="w", relief="sunken", padding=(6, 2)
        )
        status_bar.grid(row=1, column=0, columnspan=2, sticky="ew")

        # Keyboard shortcuts
        self.root.bind_all("<Control-z>", lambda e: self.undo())
        self.root.bind_all("<Control-y>", lambda e: self.redo())

    def _on_mousewheel(self, event):
        delta = -1 if event.delta > 0 else 1 if event.delta < 0 else 0
        if delta:
            self.canvas.yview_scroll(delta, "units")

    def validate_numeric(self, value_if_allowed: str, widget_name: str):
        widget = self.root.nametowidget(widget_name)
        if value_if_allowed in ("", "-", "+", ".", "-.", "+."):
            widget.configure(style="Numeric.TEntry")
            return True
        try:
            float(value_if_allowed)
            widget.configure(style="Numeric.TEntry")
        except ValueError:
            widget.configure(style="NumericInvalid.TEntry")
            self.root.bell()
        return True

    def setup_sweep_autocomplete(self, cb: ttk.Combobox):
        all_values = list(self.labels)
        cb._all_values = all_values  # type: ignore[attr-defined]
        cb["values"] = all_values

        def refresh_dropdown(filtered):
            cb["values"] = filtered if filtered else all_values

        def on_keyrelease(event):
            nav_keys = {"Up", "Down", "Left", "Right", "Home", "End", "Prior", "Next", "Return", "Escape", "Tab"}
            if event.keysym in nav_keys:
                return
            typed = cb.get().strip().lower()
            refresh_dropdown(all_values if not typed else [x for x in all_values if typed in x.lower()])
            try:
                cb.event_generate("<Down>")
            except Exception:
                pass

        cb.bind("<KeyRelease>", on_keyrelease)
        cb.bind("<FocusIn>", lambda e: refresh_dropdown(all_values))

    # -----------------------------
    # Layout builders
    # -----------------------------

    def build_layout(self):
        row = 0
        sf = self.scrollable_frame

        # Banner
        banner_frame = ttk.Frame(sf)
        banner_frame.grid(row=row, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="ew", pady=(5, 4))
        ttk.Label(banner_frame, text="T3ST: GUI Studio",
                  font=("Helvetica", 18, "bold"), anchor="center").pack(fill="x", pady=(10, 2))
        ttk.Label(banner_frame, text="Test Particle Simulation Setup  ·  Manual sweeps  ·  Database-driven runs",
                  font=("Helvetica", 10, "normal"), anchor="center", foreground="gray").pack(fill="x", pady=(0, 8))
        row += 1

        ttk.Separator(sf, orient="horizontal").grid(
            row=row, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="ew", pady=(0, 6)
        )
        row += 1

        self.mode_notebook = ttk.Notebook(sf)
        self.mode_notebook.grid(row=row, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="nsew", padx=2, pady=(4, 10))
        row += 1

        self.manual_tab = ttk.Frame(self.mode_notebook)
        self.database_tab = ttk.Frame(self.mode_notebook)
        self.mode_notebook.add(self.manual_tab, text="Manual")
        self.mode_notebook.add(self.database_tab, text="Database")
        self.mode_notebook.bind("<<NotebookTabChanged>>", self.on_tab_change)

        self._build_manual_tab()
        self._build_database_tab()
        self._build_action_buttons(sf, row)

    def _build_manual_tab(self):
        tab = self.manual_tab
        ttk.Label(tab, text="Parameter Set", style="Section.TLabel").grid(
            row=0, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="w", padx=2, pady=(8, 0)
        )

        # Row 1: Scenario + Manual sim number
        ttk.Label(tab, text="Scenario:").grid(row=1, column=0, sticky="e", padx=(2, 2), pady=4)
        self.scenario_var = tk.StringVar()
        self.scenario_dropdown = ttk.Combobox(
            tab, textvariable=self.scenario_var, values=self.valid_scenarios, state="readonly", width=22
        )
        self.scenario_dropdown.grid(row=1, column=1, padx=(0, 12), pady=4, sticky="w")
        self.scenario_dropdown.set(self.valid_scenarios[0])
        self.scenario_dropdown.bind("<<ComboboxSelected>>", self.on_scenario_change)
        ToolTip(self.scenario_dropdown, "Load a named scenario to pre-fill all parameters below.")

        ttk.Label(tab, text="Sim no.:").grid(row=1, column=2, sticky="e", padx=(2, 2), pady=4)
        self.manual_sim_no_var = tk.StringVar(value="1")
        self.manual_sim_no_entry = ttk.Spinbox(
            tab, from_=1, to=9999, textvariable=self.manual_sim_no_var, width=10
        )
        self.manual_sim_no_entry.grid(row=1, column=3, padx=(0, 12), pady=4, sticky="w")
        ToolTip(self.manual_sim_no_entry, "Output: Sims/Sim_NN.dat")

        # Row 2: Filter
        ttk.Label(tab, text="Filter:").grid(row=2, column=0, sticky="e")
        self.filter_var = tk.StringVar()
        filter_entry = ttk.Entry(tab, textvariable=self.filter_var, width=24)
        filter_entry.grid(row=2, column=1, sticky="w", padx=5, pady=(0, 8))
        ToolTip(filter_entry, "Type to filter parameter names (substring match).")
        self.filter_var.trace_add("write", lambda *_: self.apply_filter())

        self.params_frame = ttk.Frame(tab)
        self.params_frame.grid(row=3, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="ew")
        self._build_param_widgets()

        ttk.Separator(tab, orient="horizontal").grid(
            row=4, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="ew", pady=10
        )
        ttk.Label(tab, text="Sweep Configuration", style="Section.TLabel").grid(
            row=5, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="w", padx=2, pady=(10, 0)
        )
        self._build_sweep_rows(tab)

    def _build_param_widgets(self):
        """Populate self.params_frame with labelled input widgets for every parameter."""
        def tooltip_for(p: ParamDef) -> str:
            parts = []
            if p.help:     parts.append(p.help)
            if p.units:    parts.append(f"Units: {p.units}")
            if p.default is not None and p.default != "": parts.append(f"Default: {p.default}")
            if p.min is not None: parts.append(f"Min: {p.min}")
            if p.max is not None: parts.append(f"Max: {p.max}")
            if p.kind == "enum" and p.choices:
                opts = ", ".join(f"{c.value}={c.label}" for c in p.choices[:6])
                if len(p.choices) > 6:
                    opts += ", …"
                parts.append(f"Choices: {opts}")
            return "\n".join(parts) if parts else f"Enter value for {p.name}"

        # Group parameters
        group_order, group_to_params = [], {}
        for p in self.params:
            g = p.group or "Ungrouped"
            if g not in group_to_params:
                group_to_params[g] = []
                group_order.append(g)
            group_to_params[g].append(p)

        total_param_cols = UI_PARAM_COLS_PER_ROW * 2
        for col in range(total_param_cols):
            self.params_frame.columnconfigure(col, weight=0)

        # Track collapse state and per-group widget rows for show/hide
        self._group_collapsed: Dict[str, bool] = {}
        self._group_param_rows: Dict[str, List[Tuple[ttk.Label, tk.Widget]]] = {}
        self.param_group_headers: List[ttk.Label] = []

        cur_row = 0
        for g in group_order:
            self._group_collapsed[g] = False
            self._group_param_rows[g] = []

            # Clickable header button acting as collapse toggle
            hdr_text = f"▼  {g}"
            hdr = ttk.Button(
                self.params_frame, text=hdr_text,
                style="GroupHeader.TLabel",
                command=lambda grp=g: self._toggle_group(grp),
            )
            hdr.grid(row=cur_row, column=0, columnspan=total_param_cols, sticky="w", padx=4, pady=(8, 2))
            self.param_group_headers.append(hdr)
            cur_row += 1

            params_in_group = group_to_params[g]
            for i, p in enumerate(params_in_group):
                r = cur_row + (i // UI_PARAM_COLS_PER_ROW)
                c = (i % UI_PARAM_COLS_PER_ROW) * 2

                lbl = ttk.Label(self.params_frame, text=p.name, style="Param.TLabel")
                lbl.grid(row=r, column=c, padx=(2, 2), pady=4, sticky="e")
                ToolTip(lbl, tooltip_for(p))

                if p.kind == "enum":
                    values = [ch.label for ch in (p.choices or [])]
                    cb = ttk.Combobox(self.params_frame, values=values, width=UI_COMBO_WIDTH, state="readonly")
                    cb.grid(row=r, column=c + 1, padx=(2, 6), pady=4, sticky="w")
                    pw = ParamWidget(
                        pdef=p, label=lbl, widget=cb,
                        enum_label_to_value={ch.label: int(ch.value) for ch in (p.choices or [])},
                        enum_value_to_label={int(ch.value): ch.label for ch in (p.choices or [])},
                    )
                elif p.ui == "checkbox" and p.kind == "int":
                    v = tk.IntVar(value=0)
                    chk = ttk.Checkbutton(self.params_frame, variable=v)
                    chk.grid(row=r, column=c + 1, padx=(2, 6), pady=4, sticky="w")
                    pw = ParamWidget(pdef=p, label=lbl, widget=chk, var=v)
                else:
                    validate = "key" if p.kind in ("int", "float") else "none"
                    ent = ttk.Entry(
                        self.params_frame, width=UI_ENTRY_WIDTH,
                        validate=validate,
                        validatecommand=self.vcmd if validate == "key" else None,
                        style="Numeric.TEntry",
                    )
                    ent.grid(row=r, column=c + 1, padx=(2, 6), pady=4, sticky="w")
                    pw = ParamWidget(pdef=p, label=lbl, widget=ent)

                if p.default is not None and p.default != "":
                    pw.set_value(p.default)
                self.param_widgets.append(pw)
                self.param_rows.append((lbl, pw.widget))
                self._group_param_rows[g].append((lbl, pw.widget))

            cur_row += (len(params_in_group) + UI_PARAM_COLS_PER_ROW - 1) // UI_PARAM_COLS_PER_ROW

    def _toggle_group(self, group: str):
        """Collapse or expand a parameter group."""
        self._group_collapsed[group] = not self._group_collapsed[group]
        collapsed = self._group_collapsed[group]
        arrow = "▶" if collapsed else "▼"
        # Update header button text
        for hdr in self.param_group_headers:
            if hdr.cget("text").endswith(f"  {group}"):
                hdr.config(text=f"{arrow}  {group}")
                break
        for lbl, widget in self._group_param_rows.get(group, []):
            if collapsed:
                lbl.grid_remove()
                widget.grid_remove()
            else:
                lbl.grid()
                widget.grid()

    def _build_sweep_rows(self, tab):
        self.sweeps = []
        for si in range(1, 4):
            sweep_row = 5 + si  # rows 0-5 now used by header/scenario/filter/params/sep/sweep-label
            ttk.Label(tab, text=f"Param #{si}:").grid(row=sweep_row, column=0, sticky="e")
            cb = ttk.Combobox(tab, values=self.labels, width=18)
            cb.grid(row=sweep_row, column=1)
            self.setup_sweep_autocomplete(cb)

            ttk.Label(tab, text="Points:").grid(row=sweep_row, column=2, sticky="e")
            pts = ttk.Entry(tab, width=10, style="Numeric.TEntry")
            pts.grid(row=sweep_row, column=3)

            ttk.Label(tab, text="Init:").grid(row=sweep_row, column=4, sticky="e")
            init = ttk.Entry(tab, width=10, style="Numeric.TEntry")
            init.grid(row=sweep_row, column=5)

            ttk.Label(tab, text="Final:").grid(row=sweep_row, column=6, sticky="e")
            fin = ttk.Entry(tab, width=10, style="Numeric.TEntry")
            fin.grid(row=sweep_row, column=7)

            sw = SweepWidgets(label=ttk.Label(tab, text=f"Sweep {si}"), param=cb, points=pts, init=init, final=fin)
            sw.attach_tooltips()
            self.sweeps.append(sw)

    def _build_database_tab(self):
        tab = self.database_tab

        # Section header — tooltip carries the explanation
        _db_section_lbl = ttk.Label(tab, text="Database Selection  ⓘ", style="Section.TLabel", cursor="question_arrow")
        _db_section_lbl.grid(row=0, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="w", padx=2, pady=(8, 4))
        ToolTip(_db_section_lbl,
                "Pick a database (defines the parameter space axes) and a selection (how to sample it).\n"
                "Choose a base scenario to fill all non-varying parameters, then tweak individual values\n"
                "in the Parameter Overrides section below before exporting.")

        # Row 1: Database | Selection | run-count badge
        ttk.Label(tab, text="Database:").grid(row=1, column=0, sticky="e")
        self.database_dropdown = ttk.Combobox(
            tab, textvariable=self.database_var,
            values=sorted(self.database_by_name.keys()), state="readonly", width=28
        )
        self.database_dropdown.grid(row=1, column=1, padx=2, pady=4, sticky="w")
        self.database_dropdown.bind("<<ComboboxSelected>>", self.on_database_change)

        ttk.Label(tab, text="Selection:").grid(row=1, column=2, sticky="e")
        self.selection_dropdown = ttk.Combobox(
            tab, textvariable=self.selection_var, values=[], state="readonly", width=28
        )
        self.selection_dropdown.grid(row=1, column=3, padx=2, pady=4, sticky="w")
        self.selection_dropdown.bind("<<ComboboxSelected>>", self.on_selection_change)

        self._run_count_var = tk.StringVar(value="")
        self._run_count_badge = tk.Label(
            tab, textvariable=self._run_count_var,
            font=("Helvetica", 11, "bold"),
            foreground="white", background="#1a7abf",
            relief="flat", padx=10, pady=3,
            borderwidth=0,
        )
        self._run_count_badge.grid(row=1, column=4, padx=(12, 4), pady=4, sticky="w")

        # Row 2: Base scenario + DB sim number
        ttk.Label(tab, text="Base scenario:").grid(row=2, column=0, sticky="e")
        self.db_scenario_dropdown = ttk.Combobox(
            tab, textvariable=self.db_scenario_var,
            values=self.valid_scenarios, state="readonly", width=28
        )
        self.db_scenario_dropdown.grid(row=2, column=1, padx=2, pady=4, sticky="w")
        ToolTip(
            self.db_scenario_dropdown,
            "Supplies values for all parameters not swept by the database axes.\n"
            "Defaults to the scenario declared in databases.yaml; freely overridable here."
        )
        self.db_scenario_dropdown.bind("<<ComboboxSelected>>", self._on_db_scenario_change)

        ttk.Label(tab, text="Sim no.:").grid(row=2, column=2, sticky="e", padx=(2, 2), pady=4)
        self.db_sim_no_var = tk.StringVar(value="1")
        self.db_sim_no_entry = ttk.Spinbox(
            tab, from_=1, to=9999, textvariable=self.db_sim_no_var, width=10
        )
        self.db_sim_no_entry.grid(row=2, column=3, padx=(0, 12), pady=4, sticky="w")
        ToolTip(self.db_sim_no_entry, "Output: Sims/DB_NN.dat")

        # Live status labels (these are dynamic info, not explanatory — keep visible)
        self.db_info_label = ttk.Label(tab, text="Database: —", justify="left", wraplength=900)
        self.db_info_label.grid(row=3, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="w", padx=4, pady=(2, 0))

        self.sel_info_label = ttk.Label(tab, text="Selection: —", justify="left", wraplength=900)
        self.sel_info_label.grid(row=4, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="w", padx=4, pady=(0, 4))

        # --- Axis Inspector (collapsible) ---
        self._axis_inspector_collapsed = False
        self._axis_inspector_btn = ttk.Button(
            tab, text="▼  Axis Inspector", style="GroupHeader.TLabel",
            command=self._toggle_axis_inspector,
        )
        self._axis_inspector_btn.grid(row=5, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="ew", padx=2, pady=(8, 0))

        cols = ("Parameter", "Min", "Max", "N", "Step size")
        self._axis_inspector_frame = ttk.Frame(tab)
        self._axis_inspector_frame.grid(row=6, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="ew", padx=4, pady=(2, 6))
        self._axis_inspector_frame.columnconfigure(0, weight=1)

        self._axis_tree = ttk.Treeview(self._axis_inspector_frame, columns=cols, show="headings", height=4)
        for col in cols:
            self._axis_tree.heading(col, text=col)
            self._axis_tree.column(col, width=140, anchor="center")
        self._axis_tree.column("Parameter", width=200, anchor="w")
        self._axis_tree.grid(row=0, column=0, sticky="ew")

        _axis_scroll = ttk.Scrollbar(self._axis_inspector_frame, orient="vertical", command=self._axis_tree.yview)
        self._axis_tree.configure(yscrollcommand=_axis_scroll.set)
        _axis_scroll.grid(row=0, column=1, sticky="ns")

        # --- Parameter Overrides (collapsible) ---
        self._overrides_collapsed = False
        self._overrides_btn = ttk.Button(
            tab, text="▼  Parameter Overrides  ⓘ", style="GroupHeader.TLabel",
            command=self._toggle_overrides,
        )
        self._overrides_btn.grid(row=7, column=0, columnspan=6, sticky="ew", padx=2, pady=(8, 0))
        ToolTip(self._overrides_btn,
                "All simulation parameters are shown here.\n\n"
                "White label  — not part of the database; freely editable.\n"
                "Blue label   — database axis frozen to one value for this selection; read-only.\n"
                "Steel label  — database axis that varies run-to-run; shown for reference, read-only.\n"
                "Amber label  — you have edited this value away from the scenario default.")

        # "Copy to Manual" button on the same visual row as the collapse toggle
        ttk.Button(tab, text="Copy to Manual", command=self._copy_db_to_manual).grid(
            row=7, column=6, padx=(4, 2), pady=(8, 0), sticky="e"
        )

        # Collapsible body: canvas + scrollbar
        self._overrides_body_frame = ttk.Frame(tab)
        self._overrides_body_frame.grid(row=8, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="ew", padx=4, pady=(2, 6))
        self._overrides_body_frame.columnconfigure(0, weight=1)

        self._override_canvas = tk.Canvas(self._overrides_body_frame, highlightthickness=0, height=220)
        _ov_scroll = ttk.Scrollbar(self._overrides_body_frame, orient="vertical", command=self._override_canvas.yview)
        self._override_canvas.configure(yscrollcommand=_ov_scroll.set)
        self._override_canvas.grid(row=0, column=0, sticky="ew")
        _ov_scroll.grid(row=0, column=1, sticky="ns")

        self._db_override_frame = ttk.Frame(self._override_canvas)
        self._db_override_frame_id = self._override_canvas.create_window(
            (0, 0), window=self._db_override_frame, anchor="nw"
        )
        self._db_override_frame.bind(
            "<Configure>",
            lambda e: self._override_canvas.configure(scrollregion=self._override_canvas.bbox("all"))
        )
        self._override_canvas.bind(
            "<Configure>",
            lambda e: self._override_canvas.itemconfig(self._db_override_frame_id, width=e.width)
        )
        self._override_canvas.bind("<Enter>", self._bind_override_scroll)
        self._override_canvas.bind("<Leave>", self._unbind_override_scroll)

        # Scenario summary — tooltip on the label, keeps one compact info line visible
        _scen_lbl = ttk.Label(tab, text="Scenario summary  ⓘ", style="Section.TLabel", cursor="question_arrow")
        _scen_lbl.grid(row=9, column=0, columnspan=3, sticky="w", padx=2, pady=(6, 0))
        ToolTip(_scen_lbl, "Shows the active base scenario and which parameters are swept by the database axes.")
        self.database_scenario_info_label = ttk.Label(
            tab, text="—", justify="left", wraplength=900
        )
        self.database_scenario_info_label.grid(row=9, column=3, columnspan=UI_TOTAL_COLUMNS - 3, sticky="w", padx=4, pady=(6, 0))

    def _bind_override_scroll(self, event=None):
        self._override_canvas.bind_all("<MouseWheel>", self._on_override_mousewheel)
        self._override_canvas.bind_all("<Button-4>", lambda e: self._override_canvas.yview_scroll(-1, "units"))
        self._override_canvas.bind_all("<Button-5>", lambda e: self._override_canvas.yview_scroll(1, "units"))

    def _unbind_override_scroll(self, event=None):
        self._override_canvas.unbind_all("<MouseWheel>")
        self._override_canvas.unbind_all("<Button-4>")
        self._override_canvas.unbind_all("<Button-5>")
        # Restore main canvas scroll bindings
        self.canvas.bind_all("<MouseWheel>", self._on_mousewheel)
        self.canvas.bind_all("<Button-4>", lambda e: self.canvas.yview_scroll(-1, "units"))
        self.canvas.bind_all("<Button-5>", lambda e: self.canvas.yview_scroll(1, "units"))

    def _on_override_mousewheel(self, event):
        delta = -1 if event.delta > 0 else 1 if event.delta < 0 else 0
        if delta:
            self._override_canvas.yview_scroll(delta, "units")

    def _toggle_axis_inspector(self):
        self._axis_inspector_collapsed = not self._axis_inspector_collapsed
        if self._axis_inspector_collapsed:
            self._axis_inspector_frame.grid_remove()
            self._axis_inspector_btn.config(text="▶  Axis Inspector")
        else:
            self._axis_inspector_frame.grid()
            self._axis_inspector_btn.config(text="▼  Axis Inspector")

    def _toggle_overrides(self):
        self._overrides_collapsed = not self._overrides_collapsed
        if self._overrides_collapsed:
            self._overrides_body_frame.grid_remove()
            self._overrides_btn.config(text="▶  Parameter Overrides")
        else:
            self._overrides_body_frame.grid()
            self._overrides_btn.config(text="▼  Parameter Overrides")

    def _copy_db_to_manual(self):
        """Push current DB scenario + overrides into the Manual tab, then switch to it."""
        db_name = self.database_var.get().strip()
        db = self.database_by_name.get(db_name)
        scenario_name = self.db_scenario_var.get().strip()
        if db is None or not scenario_name:
            messagebox.showwarning("Copy to Manual", "Please select a database and scenario first.")
            return

        # Build merged value map: scenario → overrides (axis values not applied — those are db-specific)
        try:
            value_map = self.get_scenario_value_map(scenario_name)
        except Exception as e:
            messagebox.showerror("Copy to Manual", f"Could not load scenario: {e}")
            return

        for pname, pw in self._db_override_widgets.items():
            v, err = pw.get_export_value()
            if err:
                messagebox.showerror("Copy to Manual", f"Invalid override for {pname!r}: {err}")
                return
            pdef = self.param_by_name.get(pname)
            if pdef and pdef.kind == "enum":
                value_map[pname] = int(v)
            elif pdef and pdef.kind == "int":
                value_map[pname] = int(v)
            else:
                try:
                    value_map[pname] = float(v)
                except Exception:
                    value_map[pname] = v

        self._push_undo()
        for pw in self.param_widgets:
            val = value_map.get(pw.pdef.name)
            if val is not None:
                pw.set_value(val)

        # Switch to the Manual tab
        self.mode_notebook.select(self.manual_tab)
        self.set_status(f"Copied DB scenario '{scenario_name}' + overrides → Manual tab")

    def _rebuild_override_table(self, db: Optional[DatabaseDef], scenario_name: str,
                                sel: Optional[SelectionDef] = None):
        """Destroy and recreate parameter widgets for the given db + scenario + selection.

        Three categories of parameters are shown:
          • Free-axis params  (vary run-to-run for this selection) — steel-blue label, read-only, shows axis range
          • Frozen-axis params (axis exists but fixed to one value for this selection) — cornflower-blue label, read-only
          • Non-axis params   (not part of the database at all) — editable, amber highlight on change
        """
        frame = self._db_override_frame
        if frame is None:
            return

        for child in frame.winfo_children():
            child.destroy()
        self._db_override_widgets.clear()

        if db is None:
            ttk.Label(frame, text="(no database selected)").grid(row=0, column=0, sticky="w", padx=4, pady=2)
            return

        all_axis_names = {a.name for a in db.varying_parameters}
        axis_map = db.axis_map()

        # Determine which axes are frozen for the current selection
        frozen_axis_names: set = set()
        if sel is not None:
            try:
                fixed = resolved_fixed_indices(db, sel, self.scenarios, self.params)
                frozen_axis_names = set(fixed.keys())
            except Exception:
                pass
        free_axis_names = all_axis_names - frozen_axis_names

        # Fetch scenario values for pre-filling and change-detection
        scenario_values: Dict[str, Any] = {}
        if scenario_name in self.scenarios:
            raw = self.scenarios[scenario_name]
            if isinstance(raw, dict):
                scenario_values = raw
            elif isinstance(raw, list) and len(raw) == len(self.labels):
                scenario_values = dict(zip(self.labels, raw))

        COLS = UI_PARAM_COLS_PER_ROW
        total_cols = COLS * 2
        for col in range(total_cols):
            frame.columnconfigure(col, weight=0)

        # Colour scheme — consistent across light/dark themes
        amber_fg       = "#e07b00"   # user-modified override
        free_axis_fg   = "#4a9eda"   # freely varying axis (runs differ)
        frozen_axis_fg = "#7ab8e8"   # axis present but frozen to one value for this selection

        # All params in their defined order, grouped
        group_order: List[str] = []
        group_to_params: Dict[str, List[ParamDef]] = {}
        for p in self.params:
            g = p.group or "Ungrouped"
            if g not in group_to_params:
                group_to_params[g] = []
                group_order.append(g)
            group_to_params[g].append(p)

        cur_row = 0
        for g_idx, g in enumerate(group_order):
            if g_idx > 0:
                ttk.Separator(frame, orient="horizontal").grid(
                    row=cur_row, column=0, columnspan=total_cols, sticky="ew", pady=(4, 0)
                )
                cur_row += 1

            ttk.Label(frame, text=g, style="GroupHeader.TLabel").grid(
                row=cur_row, column=0, columnspan=total_cols, sticky="w", padx=4, pady=(4, 1)
            )
            cur_row += 1

            params_in_group = group_to_params[g]
            for i, p in enumerate(params_in_group):
                r = cur_row + (i // COLS)
                c = (i % COLS) * 2

                is_free_axis   = p.name in free_axis_names
                is_frozen_axis = p.name in frozen_axis_names
                is_axis        = is_free_axis or is_frozen_axis

                # Label colour
                if is_free_axis:
                    lbl_fg = free_axis_fg
                elif is_frozen_axis:
                    lbl_fg = frozen_axis_fg
                else:
                    lbl_fg = ""   # theme default

                lbl = ttk.Label(frame, text=p.name, style="Param.TLabel", foreground=lbl_fg)
                lbl.grid(row=r, column=c, padx=(2, 2), pady=1, sticky="e")

                # Tooltip — add axis context for axis params
                tip_parts = []
                if p.help:    tip_parts.append(p.help)
                if p.units:   tip_parts.append(f"Units: {p.units}")
                if p.default is not None and p.default != "":
                    tip_parts.append(f"Default: {p.default}")
                if p.min is not None: tip_parts.append(f"Min: {p.min}")
                if p.max is not None: tip_parts.append(f"Max: {p.max}")
                if is_free_axis:
                    ax = axis_map[p.name]
                    tip_parts.append(f"⟳ Freely varying axis: [{ax.min:.6g} … {ax.max:.6g}], {ax.nvals} pts")
                    tip_parts.append("Value changes each run — shown here for reference only.")
                elif is_frozen_axis:
                    ax = axis_map[p.name]
                    tip_parts.append(f"⬡ Frozen axis: fixed to one point in [{ax.min:.6g} … {ax.max:.6g}]")
                    tip_parts.append("Value is axis-controlled for this selection — shown read-only.")
                ToolTip(lbl, "\n".join(tip_parts) if tip_parts else f"Value for {p.name}")

                # Fill value: for axis params use the axis value if frozen, else scenario value
                if is_frozen_axis:
                    try:
                        fixed_idx = resolved_fixed_indices(db, sel, self.scenarios, self.params)[p.name]
                        fill_val = float(build_axis_values(axis_map[p.name])[fixed_idx])
                    except Exception:
                        fill_val = scenario_values.get(p.name, p.default if p.default is not None else "")
                elif is_free_axis:
                    # Show the scenario/default value as a reference — not exported
                    fill_val = scenario_values.get(p.name, p.default if p.default is not None else "")
                else:
                    fill_val = scenario_values.get(p.name, p.default if p.default is not None else "")

                if p.kind == "enum":
                    values_list = [ch.label for ch in (p.choices or [])]
                    cb = ttk.Combobox(frame, values=values_list, width=UI_COMBO_WIDTH,
                                      state="disabled" if is_axis else "readonly")
                    cb.grid(row=r, column=c + 1, padx=(2, 6), pady=1, sticky="w")
                    pw = ParamWidget(
                        pdef=p, label=lbl, widget=cb,
                        enum_label_to_value={ch.label: int(ch.value) for ch in (p.choices or [])},
                        enum_value_to_label={int(ch.value): ch.label for ch in (p.choices or [])},
                    )
                    if not is_axis:
                        def _on_enum_change(event, _lbl=lbl, _fill=fill_val, _cb=cb):
                            _lbl.config(foreground=amber_fg if _cb.get() != str(_fill) else "")
                        cb.bind("<<ComboboxSelected>>", _on_enum_change)

                elif p.ui == "checkbox" and p.kind == "int":
                    v = tk.IntVar(value=0)
                    chk = ttk.Checkbutton(frame, variable=v,
                                          state="disabled" if is_axis else "normal")
                    chk.grid(row=r, column=c + 1, padx=(2, 6), pady=1, sticky="w")
                    pw = ParamWidget(pdef=p, label=lbl, widget=chk, var=v)
                    if not is_axis:
                        try:
                            orig = int(float(fill_val)) if fill_val != "" else 0
                        except Exception:
                            orig = 0
                        def _on_check_change(*_, _lbl=lbl, _v=v, _orig=orig):
                            _lbl.config(foreground=amber_fg if _v.get() != _orig else "")
                        v.trace_add("write", _on_check_change)

                else:
                    validate = "key" if (p.kind in ("int", "float") and not is_axis) else "none"
                    ent = ttk.Entry(
                        frame, width=UI_ENTRY_WIDTH,
                        validate=validate,
                        validatecommand=self.vcmd if validate == "key" else None,
                        style="Numeric.TEntry",
                    )
                    ent.grid(row=r, column=c + 1, padx=(2, 6), pady=1, sticky="w")
                    if is_axis:
                        ent.config(state="readonly")
                    pw = ParamWidget(pdef=p, label=lbl, widget=ent)
                    if not is_axis:
                        orig_str = str(fill_val)
                        def _on_entry_change(*_, _lbl=lbl, _ent=ent, _orig=orig_str):
                            _lbl.config(foreground=amber_fg if _ent.get() != _orig else "")
                        ent.bind("<KeyRelease>", _on_entry_change)

                pw.set_value(fill_val)
                # Only register editable (non-axis) params in the export dict
                if not is_axis:
                    self._db_override_widgets[p.name] = pw

            cur_row += (len(params_in_group) + COLS - 1) // COLS

        self._override_canvas.yview_moveto(0)

    def _build_action_buttons(self, sf, row: int):
        ttk.Separator(sf, orient="horizontal").grid(
            row=row, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="ew", pady=(4, 0)
        )
        row += 1
        buttons = [
            ("Export Current Tab",  self.export_current_mode),
            ("Load Parameters",     self.load_parameters),
            ("Clear All Fields",    self.clear_all_fields),
            ("Reset to Scenario",   self.reset_to_scenario_defaults),
        ]
        for col, (label, cmd) in enumerate(buttons):
            ttk.Button(sf, text=label, command=cmd).grid(row=row, column=col, padx=4, pady=12)

        self.theme_toggle_btn = ttk.Button(
            sf, text="Switch to Dark Mode" if not self.dark_mode else "Switch to Light Mode",
            command=self.toggle_theme
        )
        self.theme_toggle_btn.grid(row=row, column=len(buttons), padx=(20, 4), pady=12)
        ttk.Button(sf, text="Quit", command=self.root.quit).grid(
            row=row, column=len(buttons) + 1, padx=4, pady=12
        )

    def initialize_database_controls(self):
        db_names = sorted(self.database_by_name.keys())
        if db_names:
            self.database_var.set(db_names[0])
            self.on_database_change()
        self.update_database_scenario_summary()
        self.on_mode_change()

    def apply_filter(self):
        q = self.filter_var.get().strip().lower()
        visible_by_group: Dict[str, bool] = {}
        for (lbl, entry), pw in zip(self.param_rows, self.param_widgets):
            g = pw.pdef.group or "Ungrouped"
            collapsed = self._group_collapsed.get(g, False)
            visible = (not q) or (q in lbl.cget("text").lower())
            visible_by_group[g] = visible_by_group.get(g, False) or visible
            # Only show if both the filter matches AND the group is not collapsed
            show = visible and not collapsed
            (lbl.grid if show else lbl.grid_remove)()
            (entry.grid if show else entry.grid_remove)()

        for hdr in getattr(self, "param_group_headers", []):
            # Extract group name from "▼  GroupName" / "▶  GroupName"
            hdr_text = hdr.cget("text")
            group_name = hdr_text[3:] if len(hdr_text) > 3 else hdr_text
            if visible_by_group.get(group_name, False) or not q:
                hdr.grid()
            else:
                hdr.grid_remove()

    def on_tab_change(self, event=None):
        current = self.mode_notebook.select()
        self.mode_var.set("database" if current == str(self.database_tab) else "manual")
        self.on_mode_change()

    def update_database_scenario_summary(self):
        if not hasattr(self, "database_scenario_info_label"):
            return
        if self.mode_var.get().strip() == "database":
            db_name = self.database_var.get().strip()
            if db_name in self.database_by_name:
                db = self.database_by_name[db_name]
                axis_names = [a.name for a in db.varying_parameters]
                chosen_scenario = self.db_scenario_var.get().strip() or db.base_scenario
                self.database_scenario_info_label.config(
                    text=f"Base scenario summary: {chosen_scenario} | "
                         f"varying parameters={', '.join(axis_names) or '-'}"
                )
                return
        scenario = self.scenario_var.get().strip()
        self.database_scenario_info_label.config(
            text=f"Base scenario summary: active scenario={scenario or '-'}"
        )

    def on_mode_change(self):
        if self.mode_var.get().strip() == "manual":
            self.sel_info_label.config(text="Selection info: manual tab active")
        else:
            self.on_selection_change()
        self.update_database_scenario_summary()

    # -----------------------------
    # Scenario / parameters
    # -----------------------------

    def clear_all_fields(self):
        self._push_undo()
        for pw in self.param_widgets:
            if pw.pdef.kind == "enum" and isinstance(pw.widget, ttk.Combobox):
                pw.widget.set("")
            elif pw.pdef.ui == "checkbox" and isinstance(pw.var, tk.IntVar):
                pw.var.set(0)
            elif isinstance(pw.widget, ttk.Entry):
                pw.widget.delete(0, tk.END)

    def on_scenario_change(self, event=None):
        scenario = self.scenario_var.get().strip()
        if scenario not in self.scenarios:
            return
        if event is not None:   # UI-triggered — push undo here
            self._push_undo()
        values = self.scenarios[scenario]

        if isinstance(values, dict):
            for pw in self.param_widgets:
                pw.set_value(values.get(pw.pdef.name, pw.pdef.default if pw.pdef.default is not None else ""))
            self.set_status(f"Loaded scenario: {scenario}")
            return

        if isinstance(values, list):
            if len(values) != len(self.labels):
                messagebox.showerror(
                    "Scenario error",
                    f"Scenario '{scenario}' has {len(values)} values, expected {len(self.labels)}."
                )
                return
            for i, pw in enumerate(self.param_widgets):
                pw.set_value(values[i])
            return

        messagebox.showerror("Scenario error", f"Scenario '{scenario}' has unsupported format: {type(values)}")

    def reset_to_scenario_defaults(self):
        self._push_undo()
        self.on_scenario_change()
        scenario = self.scenario_var.get().strip()
        self.set_status(f"Reset to scenario: {scenario}")

    # -----------------------------
    # Status bar
    # -----------------------------

    def set_status(self, msg: str):
        """Update the bottom status bar."""
        if hasattr(self, "_status_var"):
            self._status_var.set(msg)

    # -----------------------------
    # Undo / Redo
    # -----------------------------

    def _snapshot(self) -> Dict[str, Any]:
        """Capture current widget values as a {name: raw} dict."""
        snap = {}
        for pw in self.param_widgets:
            name = pw.pdef.name
            if pw.pdef.kind == "enum" and isinstance(pw.widget, ttk.Combobox):
                snap[name] = pw.widget.get()
            elif pw.pdef.ui == "checkbox" and isinstance(pw.var, tk.IntVar):
                snap[name] = pw.var.get()
            elif isinstance(pw.widget, ttk.Entry):
                snap[name] = pw.widget.get()
            else:
                snap[name] = ""
        return snap

    def _restore_snapshot(self, snap: Dict[str, Any]):
        for pw in self.param_widgets:
            val = snap.get(pw.pdef.name, "")
            pw.set_value(val)

    def _push_undo(self):
        """Save current state onto the undo stack and clear redo."""
        self._undo_stack.append(self._snapshot())
        if len(self._undo_stack) > 50:
            self._undo_stack.pop(0)
        self._redo_stack.clear()

    def undo(self):
        if not self._undo_stack:
            self.set_status("Nothing to undo.")
            return
        self._redo_stack.append(self._snapshot())
        self._restore_snapshot(self._undo_stack.pop())
        self.set_status(f"Undo  (undo stack: {len(self._undo_stack)} left)")

    def redo(self):
        if not self._redo_stack:
            self.set_status("Nothing to redo.")
            return
        self._undo_stack.append(self._snapshot())
        self._restore_snapshot(self._redo_stack.pop())
        self.set_status(f"Redo  (redo stack: {len(self._redo_stack)} left)")

    # -----------------------------
    # Theme
    # -----------------------------

    def toggle_theme(self):
        self.dark_mode = not self.dark_mode
        sv_ttk.set_theme("dark" if self.dark_mode else "light")
        self.configure_numeric_entry_styles()
        self.theme_toggle_btn.config(
            text="Switch to Light Mode" if self.dark_mode else "Switch to Dark Mode"
        )
        self.save_config()

    # -----------------------------
    # Database UI events
    # -----------------------------

    def on_database_change(self, event=None):
        db_name = self.database_var.get().strip()
        if db_name not in self.database_by_name:
            self.selection_dropdown["values"] = []
            self.db_info_label.config(text="Database info: invalid database")
            self._refresh_axis_inspector(None)
            self._rebuild_override_table(None, "", None)
            return

        db = self.database_by_name[db_name]
        self.db_info_label.config(text=(
            f"Database info: {db.name} | base_scenario={db.base_scenario} | "
            f"dimensions={len(db.varying_parameters)} | total runs={database_total_runs(db)} | "
            f"description={db.description or '-'}"
        ))
        self._refresh_axis_inspector(db)

        # Pre-select the scenario specified in databases.yaml (user can override)
        if db.base_scenario in self.valid_scenarios:
            self.db_scenario_var.set(db.base_scenario)
        elif self.valid_scenarios:
            self.db_scenario_var.set(self.valid_scenarios[0])

        self._rebuild_override_table(db, self.db_scenario_var.get(), None)

        selection_names = sorted(s.name for s in self.selections if s.database == db_name)
        self.selection_dropdown["values"] = selection_names
        if selection_names:
            self.selection_var.set(selection_names[0])
            self.on_selection_change()
        else:
            self.selection_var.set("")
            self.sel_info_label.config(text="Selection info: no selections available for this database")
            self._run_count_var.set("")

    def _refresh_axis_inspector(self, db: Optional[DatabaseDef]):
        """Repopulate the axis Treeview for the given database (or clear it)."""
        if not hasattr(self, "_axis_tree"):
            return
        for row in self._axis_tree.get_children():
            self._axis_tree.delete(row)
        if db is None:
            self._axis_tree.configure(height=4)
            return
        for axis in db.varying_parameters:
            if axis.nvals > 1:
                step_str = f"{(axis.max - axis.min) / (axis.nvals - 1):.6g}"
            else:
                step_str = "N/A"
            self._axis_tree.insert("", "end", values=(
                axis.name,
                f"{axis.min:.6g}",
                f"{axis.max:.6g}",
                str(axis.nvals),
                step_str,
            ))
        # Auto-size: show all rows up to a max of 7, no blank space
        n = len(db.varying_parameters)
        self._axis_tree.configure(height=max(1, min(n, 7)))

    def on_selection_change(self, event=None):
        sel_name = self.selection_var.get().strip()
        if not sel_name or sel_name not in self.selection_by_name:
            self.sel_info_label.config(text="Selection info: -")
            self._run_count_var.set("")
            if hasattr(self, "_run_count_badge"):
                self._run_count_badge.grid_remove()
            return

        sel = self.selection_by_name[sel_name]
        if sel.database not in self.database_by_name:
            self.sel_info_label.config(text=f"Selection info: database {sel.database!r} not found")
            self._run_count_var.set("")
            if hasattr(self, "_run_count_badge"):
                self._run_count_badge.grid_remove()
            return

        db = self.database_by_name[sel.database]
        try:
            count = selection_count(db, sel, self.scenarios, self.params)
            self.sel_info_label.config(text=(
                f"Selection info: {sel.name} | method={sel.method} | selected runs={count} | "
                f"description={sel.description or '-'}"
            ))
            self._run_count_var.set(f"  {count:,} runs  ")
            if hasattr(self, "_run_count_badge"):
                self._run_count_badge.config(background="#1a7abf")
                self._run_count_badge.grid()
        except Exception as e:
            self.sel_info_label.config(text=f"Selection info: invalid selection ({e})")
            self._run_count_var.set("  ⚠ error  ")
            if hasattr(self, "_run_count_badge"):
                self._run_count_badge.config(background="#c0392b")
                self._run_count_badge.grid()

        # Rebuild override table so axis colouring reflects the new selection
        self._rebuild_override_table(db, self.db_scenario_var.get(), sel)

    def _on_db_scenario_change(self, event=None):
        """Called when the user changes the scenario override dropdown in database mode."""
        scenario_name = self.db_scenario_var.get().strip()
        db_name = self.database_var.get().strip()
        db = self.database_by_name.get(db_name)
        sel_name = self.selection_var.get().strip()
        sel = self.selection_by_name.get(sel_name)
        self._rebuild_override_table(db, scenario_name, sel)
        self.update_database_scenario_summary()
        self.set_status(f"Database base scenario changed to: {scenario_name}")

    # -----------------------------
    # Sweeps parsing / preview
    # -----------------------------

    def parse_sweeps(self) -> Tuple[List[Tuple[int, np.ndarray, str]], List[str]]:
        errors: List[str] = []
        sweeps: List[Tuple[int, np.ndarray, str]] = []
        seen_params: List[str] = []

        for si, sw in enumerate(self.sweeps, start=1):
            pname = sw.param.get().strip()
            if pname == "":
                continue
            if pname not in self.labels:
                errors.append(f"Sweep #{si}: invalid parameter {pname!r}"); continue
            pdef = self.param_by_name[pname]
            if pdef.kind not in ("int", "float", "enum"):
                errors.append(f"Sweep #{si} ({pname}): only int/float/enum parameters can be swept (got {pdef.kind})")
                continue
            if pname in seen_params:
                errors.append(f"Sweep #{si}: parameter {pname!r} is duplicated across sweeps"); continue
            seen_params.append(pname)

            pts_txt, init_txt, fin_txt = sw.points.get().strip(), sw.init.get().strip(), sw.final.get().strip()
            if "" in (pts_txt, init_txt, fin_txt):
                errors.append(f"Sweep #{si} ({pname}): points/init/final must all be filled"); continue

            try:
                n = int(float(pts_txt))
                if n < 1: raise ValueError()
            except Exception:
                errors.append(f"Sweep #{si} ({pname}): Points must be an integer ≥ 1"); continue

            try:
                vals = np.linspace(float(init_txt), float(fin_txt), n)
            except Exception:
                errors.append(f"Sweep #{si} ({pname}): Init/Final must be numeric"); continue

            if pdef.kind == "int":
                bad = [v for v in vals if abs(v - int(round(v))) > 1e-12]
                if bad:
                    errors.append(
                        f"Sweep #{si} ({pname}): integer sweep generated non-integer values: "
                        + ", ".join(f"{x:.8g}" for x in bad[:5])
                    )
                    continue

            elif pdef.kind == "enum":
                allowed = {int(ch.value) for ch in (pdef.choices or [])}
                bad = []
                for v in vals:
                    try:
                        xi = format_sweep_enum_value(v)
                    except Exception:
                        bad.append(v); continue
                    if xi not in allowed:
                        bad.append(v)
                if bad:
                    errors.append(
                        f"Sweep #{si} ({pname}): generated invalid enum value(s): "
                        + ", ".join(f"{x:.8g}" for x in bad[:5])
                        + f". Allowed values: {', '.join(str(x) for x in sorted(allowed))}"
                    )
                    continue

            sweeps.append((self.labels.index(pname), vals, pname))

        return sweeps, errors


    # -----------------------------
    # Export helpers
    # -----------------------------

    def get_sims_dir(self) -> Path:
        d = Path(__file__).resolve().parent / "Sims"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def write_dat_file(self, file_path: Path, rows: List[List[str]], comment: str = ""):
        with open(file_path, "w", encoding="utf-8") as f:
            if comment:
                f.write(f"# {comment}\n")
            f.write(f"{len(rows)}\n{len(self.labels)}\n")
            f.write(" ".join(self.labels) + "\n")
            for row in rows:
                f.write(" ".join(row) + "\n")

    def get_current_manual_formatted_values(self) -> List[str]:
        value_errors, formatted_values = [], [""] * len(self.param_widgets)
        for i, pw in enumerate(self.param_widgets):
            v, err = pw.get_export_value()
            if err:
                value_errors.append(err)
            else:
                formatted_values[i] = str(v)
        if value_errors:
            msg = "Please fix these value issues:\n\n" + "\n".join(value_errors[:50])
            if len(value_errors) > 50:
                msg += f"\n\n… and {len(value_errors) - 50} more"
            raise ValueError(msg)
        return formatted_values

    def get_scenario_value_map(self, scenario_name: str) -> Dict[str, Any]:
        if scenario_name not in self.scenarios:
            raise ValueError(f"Scenario {scenario_name!r} not found")
        values = self.scenarios[scenario_name]
        if isinstance(values, dict):
            return {p.name: values.get(p.name, p.default if p.default is not None else "") for p in self.params}
        if isinstance(values, list):
            if len(values) != len(self.labels):
                raise ValueError(
                    f"Scenario {scenario_name!r} has {len(values)} values, expected {len(self.labels)}"
                )
            return {name: values[i] for i, name in enumerate(self.labels)}
        raise ValueError(f"Scenario {scenario_name!r} has unsupported format")

    def build_database_run_row(self, db: DatabaseDef, idx_tuple: Tuple[int, ...]) -> List[str]:
        # Use the scenario chosen in the GUI dropdown, falling back to the YAML-specified one
        chosen_scenario = self.db_scenario_var.get().strip() if hasattr(self, "db_scenario_var") else ""
        if not chosen_scenario or chosen_scenario not in self.scenarios:
            chosen_scenario = db.base_scenario
        value_map = self.get_scenario_value_map(chosen_scenario)

        # Apply per-parameter overrides from the override table (non-varying params only)
        for pname, pw in self._db_override_widgets.items():
            v, err = pw.get_export_value()
            if err:
                raise ValueError(f"Override for parameter {pname!r}: {err}")
            # Convert back to the raw type for value_map
            pdef = self.param_by_name.get(pname)
            if pdef and pdef.kind == "enum":
                # v is the integer string; store as int
                value_map[pname] = int(v)
            elif pdef and pdef.kind == "int":
                value_map[pname] = int(v)
            else:
                value_map[pname] = float(v) if v else v

        # Then overlay the varying-axis values for this specific run
        for axis, idx in zip(database_axes_in_order(db), idx_tuple):
            raw = float(build_axis_values(axis)[int(idx)])
            pdef = self.param_by_name[axis.name]
            rng_err = validate_range(raw, pdef.min, pdef.max)
            if rng_err:
                raise ValueError(
                    f"Database {db.name!r}: generated value for {axis.name!r} is out of allowed range: {rng_err}"
                )
            value_map[axis.name] = raw

        row: List[str] = []
        for p in self.params:
            raw = value_map.get(p.name, p.default if p.default is not None else "")
            if raw == "" or raw is None:
                raise ValueError(f"Base scenario {chosen_scenario!r}: parameter {p.name!r} is empty")
            if p.kind == "enum":
                try:
                    xi = int(round(float(raw)))
                except Exception:
                    raise ValueError(f"Parameter {p.name!r}: enum value {raw!r} is not numeric")
                allowed = {int(ch.value) for ch in (p.choices or [])}
                if xi not in allowed:
                    raise ValueError(
                        f"Parameter {p.name!r}: enum value {xi} not in allowed set {sorted(allowed)}"
                    )
                row.append(str(xi))
            else:
                row.append(format_value(str(raw), p.kind))
        return row

    def get_current_database_and_selection(self) -> Tuple[DatabaseDef, SelectionDef]:
        db_name = self.database_var.get().strip()
        sel_name = self.selection_var.get().strip()
        if db_name not in self.database_by_name:
            raise ValueError("Please choose a valid database.")
        if sel_name not in self.selection_by_name:
            raise ValueError("Please choose a valid selection.")
        db = self.database_by_name[db_name]
        sel = self.selection_by_name[sel_name]
        if sel.database != db.name:
            raise ValueError(f"Selection {sel.name!r} does not belong to database {db.name!r}")
        return db, sel

    # -----------------------------
    # Export
    # -----------------------------

    def export_manual_mode(self, file_path: Path):
        import datetime
        sweeps, sweep_errors = self.parse_sweeps()
        if sweep_errors:
            raise ValueError("\n".join(sweep_errors))

        formatted_values = self.get_current_manual_formatted_values()
        total_rows = int(np.prod([len(vals) for _, vals, _ in sweeps])) if sweeps else 1
        if total_rows > MAX_EXPORT_ROWS:
            raise ValueError(
                f"Manual export would generate {total_rows} rows, which exceeds the current limit of {MAX_EXPORT_ROWS}."
            )

        rows: List[List[str]] = []
        if not sweeps:
            rows.append(formatted_values)
        else:
            def recurse(level: int, current: List[str]):
                if level == len(sweeps):
                    rows.append(current); return
                idx, vals, pname = sweeps[level]
                pdef = self.param_by_name[pname]
                for v in vals:
                    nxt = current.copy()
                    if pdef.kind == "enum":
                        xi = format_sweep_enum_value(float(v))
                        allowed = {int(ch.value) for ch in (pdef.choices or [])}
                        if xi not in allowed:
                            raise ValueError(
                                f"Sweep for {pname} generated invalid enum value {xi}. Allowed: {sorted(allowed)}"
                            )
                        nxt[idx] = str(xi)
                    else:
                        nxt[idx] = format_value(str(v), pdef.kind)
                    recurse(level + 1, nxt)
            recurse(0, formatted_values)

        scenario = self.scenario_var.get().strip() or "—"
        date_str = datetime.date.today().isoformat()
        if sweeps:
            sweep_desc = ", ".join(
                f"{pname}({vals[0]:.6g}→{vals[-1]:.6g}, {len(vals)}pts)"
                for _, vals, pname in sweeps
            )
        else:
            sweep_desc = "—"
        comment = (
            f"mode=manual | scenario={scenario} | db=— | selection=— | "
            f"sweeps={sweep_desc} | runs={total_rows} | exported={date_str}"
        )

        self.write_dat_file(file_path, rows, comment)

    def export_database_mode(self, file_path: Path):
        import datetime
        db, sel = self.get_current_database_and_selection()
        count = selection_count(db, sel, self.scenarios, self.params)
        if count > MAX_EXPORT_ROWS:
            raise ValueError(
                f"Selection {sel.name!r} would generate {count} rows, which exceeds the current limit of {MAX_EXPORT_ROWS}.\n"
                f"Please use a smaller selection."
            )
        rows = [
            self.build_database_run_row(db, tuple(int(x) for x in idx_tuple))
            for idx_tuple in selection_index_iterator(db, sel, self.scenarios, self.params)
        ]
        chosen_scenario = self.db_scenario_var.get().strip() or db.base_scenario
        date_str = datetime.date.today().isoformat()
        comment = (
            f"mode=database | scenario={chosen_scenario} | db={db.name} | selection={sel.name} | "
            f"sweeps=— | runs={count} | exported={date_str}"
        )
        self.write_dat_file(file_path, rows, comment)

    def export_current_mode(self):
        is_db = self.mode_var.get().strip() == "database"

        if is_db:
            sim_no_txt = self.db_sim_no_var.get().strip()
            if not sim_no_txt.isdigit():
                messagebox.showerror("Invalid Input", "Simulation number must be an integer.")
                return
            file_path = self.get_sims_dir() / f"DB_{int(sim_no_txt):02d}.dat"
        else:
            sim_no_txt = self.manual_sim_no_var.get().strip()
            if not sim_no_txt.isdigit():
                messagebox.showerror("Invalid Input", "Simulation number must be an integer.")
                return
            file_path = self.get_sims_dir() / f"Sim_{int(sim_no_txt):02d}.dat"

        if file_path.exists() and not messagebox.askyesno("Overwrite?", f"{file_path.name} exists. Overwrite?"):
            return
        try:
            if is_db:
                self.export_database_mode(file_path)
            else:
                self.export_manual_mode(file_path)
            messagebox.showinfo("Success", f"Parameters saved to:\n{file_path}")
            self.set_status(f"Exported → {file_path.name}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to write file:\n{e}")
            self.set_status(f"Export failed: {e}")

    # -----------------------------
    # Load
    # -----------------------------

    def load_parameters(self):
        file_path = filedialog.askopenfilename(filetypes=[("DAT Files", "*.dat")])
        if not file_path:
            return
        file_path = Path(file_path)
        try:
            raw_lines = [ln.strip() for ln in file_path.read_text(encoding="utf-8").splitlines() if ln.strip()]

            # Extract and display any leading comment lines, then strip them for parsing
            comment_lines = [ln for ln in raw_lines if ln.startswith("#")]
            lines = [ln for ln in raw_lines if not ln.startswith("#")]

            if len(lines) < 3:
                messagebox.showerror("Error", "File format not recognized."); return

            try:
                rows, cols = int(lines[0]), int(lines[1])
            except Exception:
                messagebox.showerror("Error", "Header lines corrupted or missing."); return

            if cols != len(self.labels):
                messagebox.showerror("Mismatch", f"File has {cols} parameters, expected {len(self.labels)}."); return

            third = lines[2].split()
            has_header = (len(third) == cols) and all(tok in self.labels for tok in third)

            if has_header:
                header_names = third
                if len(lines) < 4:
                    messagebox.showerror("Error", "File has a header but no data rows."); return
                data = lines[3].split()
            else:
                header_names, data = self.labels[:], third

            if len(data) != cols:
                messagebox.showerror("Mismatch", "Parameter count mismatch in first data row."); return

            name_to_val = dict(zip(header_names, data))
            self._push_undo()
            for pw in self.param_widgets:
                pw.set_value(name_to_val.get(pw.pdef.name, ""))

            info = f"Loaded parameters from {file_path}\n(total rows in file: {rows})"
            if comment_lines:
                info += f"\n\n{chr(10).join(comment_lines)}"
            messagebox.showinfo("Loaded", info)
            self.set_status(f"Loaded parameters from {file_path.name}  ({rows} row(s) in file)")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {e}")


if __name__ == "__main__":
    root = tk.Tk()
    app = SimulationGUI(root)
    root.mainloop()