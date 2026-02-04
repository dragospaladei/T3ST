"""
Define_Simulations.py (refactored)

Key improvements implemented:
- Robust config directory discovery + friendly startup errors for missing YAMLs
- Stronger scenario handling:
  * accepts YAML scenario keys that are ints (0,1,...) by normalizing to str
  * validates scenario dict keys against parameters.yaml (catches typos/extras)
- Export formatting no longer destroys small numbers:
  * floats written with a magnitude-preserving format (8 significant digits)
  * ints / flags written as integers
- Export validates all parameters and all sweep fields first; shows one aggregated error dialog
- Sweep UI built in a loop (no copy/paste for sweep 1/2/3)
- Prevents duplicate sweep parameters (Param #1/#2/#3 must be distinct)
- Preview sweep now previews ALL configured sweeps (up to 3 plots) + shows total run count
- Parameter search box to quickly filter the long parameter list (grid rows hide/show)
- Uses pathlib consistently for paths

Dependencies: tkinter, numpy, matplotlib, sv_ttk, pyyaml
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

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

UI_TOTAL_COLUMNS = 10          # overall grid columns used in the main scrollable frame
UI_PARAM_COLS_PER_ROW = 6      # number of parameter widgets per row
UI_ENTRY_WIDTH = 8             # width of entry widgets (characters)
UI_COMBO_WIDTH = 5             # width of combobox widgets (characters)
TOOLTIP_FONT_FAMILY = "tahoma"
TOOLTIP_FONT_SIZE = 11
TOOLTIP_WRAP_PX = 420


# -----------------------------
# Config discovery / YAML loading
# -----------------------------

def discover_config_dir() -> Path:
    """
    Prefer ./config next to this script; fallback to script directory.
    This makes the GUI work whether YAML files are stored in ./config/ or alongside the script.
    """
    script_dir = Path(__file__).resolve().parent
    cfg = script_dir / "config"
    return cfg if cfg.is_dir() else script_dir


CONFIG_DIR = discover_config_dir()


def load_yaml(path: Path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def require_yaml_files(config_dir: Path, filenames: List[str]) -> Dict[str, Path]:
    """
    Returns mapping filename -> resolved Path. Raises with a friendly message if missing.
    """
    missing = []
    resolved: Dict[str, Path] = {}
    for name in filenames:
        p = config_dir / name
        if not p.exists():
            missing.append(str(p))
        else:
            resolved[name] = p
    if missing:
        msg = "Missing required config file(s):\n\n" + "\n".join(missing)
        raise FileNotFoundError(msg)
    return resolved


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

        # --- wrapping / font controls ---
        wrap_px = 420  # adjust: 300-600 is typical
        label = tk.Label(
            tw,
            text=self.text,
            justify="left",
            anchor="w",
            wraplength=wrap_px,          # <-- makes it wrap
            background="#ffffe0",
            relief="solid",
            borderwidth=1,
            font=(TOOLTIP_FONT_FAMILY, TOOLTIP_FONT_SIZE, "normal") # <-- bigger font example
        )
        label.pack(ipadx=6, ipady=4)

        # Keep tooltip on-screen
        tw.update_idletasks()
        w, h = tw.winfo_reqwidth(), tw.winfo_reqheight()
        scr_w, scr_h = tw.winfo_screenwidth(), tw.winfo_screenheight()
        x = min(x, scr_w - w - 10)
        y = min(y, scr_h - h - 10)
        tw.wm_geometry(f"+{x}+{y}")

    def hide(self, event=None):
        self._cancel_scheduled()
        if self.tip:
            self.tip.destroy()
            self.tip = None


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
    ui: str = ""  # e.g. "checkbox" for 0/1 ints
    choices: Optional[List[Choice]] = None


def _as_float_or_none(x) -> Optional[float]:
    if x is None or x == "":
        return None
    try:
        return float(x)
    except Exception:
        return None


def parse_parameters_yaml(obj) -> List[ParamDef]:
    """Parse parameters.yaml.

    Supports:
      - legacy format: YAML list of parameter-name strings
      - enriched v2: dict with key "parameters" as list of dicts
    """
    if obj is None:
        return []

    # Legacy: ["t0", "tmax", ...]
    if isinstance(obj, list) and all(isinstance(x, str) for x in obj):
        return [ParamDef(name=s, kind="float") for s in obj]

    # Enriched: {schema_version: 2, parameters: [ {name: ..., kind: ...}, ...]}
    if isinstance(obj, dict) and isinstance(obj.get("parameters"), list):
        out: List[ParamDef] = []
        for i, p in enumerate(obj["parameters"]):
            if not isinstance(p, dict) or "name" not in p:
                raise ValueError(f"parameters.yaml: item #{i} must be a dict with a 'name' field")
            name = str(p["name"]).strip()
            kind = str(p.get("kind", "float")).strip()
            group = str(p.get("group", "Ungrouped") or "Ungrouped")
            units = str(p.get("units", "") or "")
            help_txt = str(p.get("help", "") or "")
            default = p.get("default", None)
            ui = str(p.get("ui", "") or "")
            pmin = _as_float_or_none(p.get("min", None))
            pmax = _as_float_or_none(p.get("max", None))

            choices_raw = p.get("choices", None)
            choices: Optional[List[Choice]] = None
            if kind == "enum":
                if not isinstance(choices_raw, list) or not choices_raw:
                    raise ValueError(f"parameters.yaml: enum param {name!r} must provide non-empty 'choices' list")
                choices = []
                for ch in choices_raw:
                    if not isinstance(ch, dict) or "value" not in ch or "label" not in ch:
                        raise ValueError(
                            f"parameters.yaml: enum param {name!r} choices must be list of {{value, label}} dicts"
                        )
                    choices.append(Choice(value=int(ch["value"]), label=str(ch["label"])) )

            out.append(
                ParamDef(
                    name=name,
                    kind=kind,
                    group=group,
                    units=units,
                    help=help_txt,
                    default=default,
                    min=pmin,
                    max=pmax,
                    ui=ui,
                    choices=choices,
                )
            )
        return out

    raise ValueError("parameters.yaml: unsupported format (expected legacy list or enriched dict with 'parameters')")


def is_float_like(text: str) -> bool:
    if text.strip() == "":
        return False
    try:
        float(text)
        return True
    except Exception:
        return False


def format_value(text: str, kind: str) -> str:
    """
    Convert from GUI string to output string.
    - float: magnitude-preserving "general" format (keeps scientific notation when needed)
    - int: integer format (no decimals)
    """
    if kind == "int":
        # Accept floats like "1.0" and cast to int if they are integer-valued
        x = float(text)
        if not np.isfinite(x):
            raise ValueError("not finite")
        xi = int(round(x))
        if abs(x - xi) > 1e-12:
            raise ValueError(f"expected integer-like value, got {text!r}")
        return str(xi)

    # float
    x = float(text)
    if not np.isfinite(x):
        raise ValueError("not finite")
    # 8 significant digits is a good compromise; preserves tiny values via scientific notation
    return f"{x:.8g}"


def validate_range(x: float, pmin: Optional[float], pmax: Optional[float]) -> Optional[str]:
    if pmin is not None and x < pmin:
        return f"must be ≥ {pmin}"
    if pmax is not None and x > pmax:
        return f"must be ≤ {pmax}"
    return None


@dataclass
class ParamWidget:
    pdef: ParamDef
    label: ttk.Label
    widget: tk.Widget
    var: Optional[tk.Variable] = None  # StringVar/IntVar for non-Entry controls
    # For enums: map label -> value and value -> label
    enum_label_to_value: Optional[Dict[str, int]] = None
    enum_value_to_label: Optional[Dict[int, str]] = None

    def set_value(self, value: Any):
        """Set widget from a raw value (from scenario / file)."""
        if value is None:
            value = ""

        if self.pdef.kind == "enum":
            # Accept int codes, numeric strings, or labels
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

        # Entry-like
        if isinstance(self.widget, ttk.Entry):
            self.widget.config(state=tk.NORMAL)
            self.widget.delete(0, tk.END)
            self.widget.insert(0, str(value))

    def get_export_value(self) -> Tuple[Optional[str], Optional[str]]:
        """Return (formatted_value, error_message)."""
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
                # range-check uses the numeric value (before formatting)
                x = float(raw)
                if not np.isfinite(x):
                    return None, f"Parameter {name}: not finite"
                rng_err = validate_range(x, self.pdef.min, self.pdef.max)
                if rng_err:
                    return None, f"Parameter {name}: {rng_err}"
                return format_value(raw, self.pdef.kind), None
            except Exception as e:
                return None, f"Parameter {name}: {e}"

        # string (rare)
        return raw, None


# -----------------------------
# Sweep definition
# -----------------------------

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


# -----------------------------
# Main GUI
# -----------------------------

class SimulationGUI:
    CONFIG_FILE = "gui_config.json"
    REQUIRED_YAMLS = ["parameters.yaml", "scenarios.yaml"]

    def __init__(self, root: tk.Tk):
        self.root = root
        self.dark_mode = False
        self._frame_window_id = None

        # Load app config (theme)
        self.load_config()

        # Window geometry
        self.root.title("T3ST GUI")
        screen_w = self.root.winfo_screenwidth()
        screen_h = self.root.winfo_screenheight()
        w = int(screen_w * 0.8)
        h = int(screen_h * 0.8)
        x = (screen_w - w) // 2
        y = (screen_h - h) // 2
        self.root.geometry(f"{w}x{h}+{x}+{y}")

        # Theme
        sv_ttk.set_theme("dark" if self.dark_mode else "light")

        # Load YAMLs (with friendly errors)
        try:
            paths = require_yaml_files(CONFIG_DIR, self.REQUIRED_YAMLS)
            params_raw = load_yaml(paths["parameters.yaml"]) or {}
            scenarios_raw = load_yaml(paths["scenarios.yaml"]) or {}
        except Exception as e:
            messagebox.showerror("Startup error", str(e))
            raise

        # Parameters (supports legacy list or enriched schema)
        try:
            self.params: List[ParamDef] = parse_parameters_yaml(params_raw)
        except Exception as e:
            messagebox.showerror("Config error", f"Invalid parameters.yaml:\n\n{e}")
            raise

        if not self.params:
            messagebox.showerror("Config error", "parameters.yaml contains no parameters.")
            raise ValueError("Empty parameters.yaml")

        self.labels: List[str] = [p.name for p in self.params]
        self.param_by_name: Dict[str, ParamDef] = {p.name: p for p in self.params}

        # Normalize scenario keys to str (handles YAML keys like 0: {...})
        self.scenarios: Dict[str, object] = {str(k): v for k, v in (scenarios_raw or {}).items()}
        self.valid_scenarios = list(self.scenarios.keys())
        if not self.valid_scenarios:
            self.valid_scenarios = ["(none)"]

        # Validate scenario dict keys (catch typos / extras)
        self._scenario_warnings: List[str] = []
        self.validate_scenarios()

        # State
        self.param_widgets: List[ParamWidget] = []
        self.param_rows: List[Tuple[ttk.Label, tk.Widget]] = []  # for filtering
        self.sweeps: List[SweepWidgets] = []

        # UI setup
        self.setup_fonts()
        self.vcmd = (self.root.register(self.validate_numeric), "%P", "%W")
        self.create_widgets()

        # Apply scenario after UI exists
        self.root.after(50, self.on_scenario_change)

        # Show scenario warnings once at startup (non-fatal)
        if self._scenario_warnings:
            messagebox.showwarning(
                "Scenario warnings",
                "Some scenarios have issues:\n\n" + "\n".join(self._scenario_warnings[:30]) +
                ("\n\n… (more)" if len(self._scenario_warnings) > 30 else "")
            )

    # -----------------------------
    # Config (theme)
    # -----------------------------

    def load_config(self):
        cfg_path = Path(self.CONFIG_FILE)
        if cfg_path.exists():
            try:
                with open(cfg_path, "r", encoding="utf-8") as f:
                    cfg = json.load(f)
                self.dark_mode = bool(cfg.get("dark_mode", False))
            except Exception:
                self.dark_mode = False

    def save_config(self):
        with open(self.CONFIG_FILE, "w", encoding="utf-8") as f:
            json.dump({"dark_mode": self.dark_mode}, f)

    # -----------------------------
    # Startup validation
    # -----------------------------

    def validate_scenarios(self):
        labels_set = set(self.labels)
        for scen_name, values in self.scenarios.items():
            # New format: dict keyed by parameter name
            if isinstance(values, dict):
                extras = [k for k in values.keys() if k not in labels_set]
                if extras:
                    self._scenario_warnings.append(
                        f"Scenario {scen_name!r}: unknown parameter key(s): {extras}"
                    )
            # Old format: list
            elif isinstance(values, list):
                if len(values) != len(self.labels):
                    self._scenario_warnings.append(
                        f"Scenario {scen_name!r}: list has {len(values)} values, expected {len(self.labels)}"
                    )
            else:
                self._scenario_warnings.append(
                    f"Scenario {scen_name!r}: unsupported type {type(values).__name__}"
                )

    # -----------------------------
    # UI
    # -----------------------------

    def setup_fonts(self):
        default_font = tkFont.nametofont("TkDefaultFont")
        default_font.configure(size=10, weight="bold")

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

        style = ttk.Style()
        style.map("TButton", background=[("active", "#4b9ce2")], foreground=[("active", "white")])
        style.configure("Section.TLabel", font=("Helvetica", 12, "bold"))
        style.configure("Param.TLabel", font=("Consolas", 10, "bold"))
        style.configure("GroupHeader.TLabel", font=("Helvetica", 11, "bold"))

        self.build_layout()

    def _on_mousewheel(self, event):
        delta = -1 if event.delta > 0 else 1 if event.delta < 0 else 0
        if delta:
            self.canvas.yview_scroll(delta, "units")

    def validate_numeric(self, value_if_allowed: str, widget_name: str):
        """Allow floats incl. scientific notation; color red when invalid."""
        widget = self.root.nametowidget(widget_name)

        # Allow empty / intermediate states while typing
        if value_if_allowed in ("", "-", "+", ".", "-.", "+."):
            widget.configure(foreground="black")
            return True

        try:
            float(value_if_allowed)
            widget.configure(foreground="black")
            return True
        except ValueError:
            widget.configure(foreground="red")
            self.root.bell()
            return True  # allow continued editing

    def build_layout(self):
        row = 0

        # Banner
        banner_frame = ttk.Frame(self.scrollable_frame)
        banner_frame.grid(row=row, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="ew", pady=(5, 10))
        banner_label = ttk.Label(
            banner_frame,
            text="📈  T3ST GUI — v2.0 ",
            font=("Helvetica", 18, "bold"),
            anchor="center",
        )
        banner_label.pack(fill="x", pady=10)
        row += 1

        # Simulation setup
        ttk.Label(self.scrollable_frame, text="Simulation Setup", style="Section.TLabel").grid(
            row=row, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="w", padx=2, pady=(10, 0)
        )
        row += 1

        ttk.Label(self.scrollable_frame, text="Simulation no.:").grid(row=row, column=0, sticky="e")
        self.sim_no_var = tk.StringVar()
        self.sim_no_entry = ttk.Spinbox(self.scrollable_frame, from_=1, to=100, textvariable=self.sim_no_var, width=10)
        self.sim_no_entry.grid(row=row, column=1, padx=2, pady=5)

        ttk.Label(self.scrollable_frame, text="Scenario:").grid(row=row, column=2, sticky="e")
        self.scenario_var = tk.StringVar()
        self.scenario_dropdown = ttk.Combobox(
            self.scrollable_frame, textvariable=self.scenario_var,
            values=self.valid_scenarios, state="readonly", width=18
        )
        self.scenario_dropdown.grid(row=row, column=3, padx=2, pady=5)
        self.scenario_dropdown.set(self.valid_scenarios[0])
        self.scenario_dropdown.bind("<<ComboboxSelected>>", self.on_scenario_change)

        self.theme_toggle_btn = ttk.Button(
            self.scrollable_frame,
            text="Switch to Light Mode" if self.dark_mode else "Switch to Dark Mode",
            command=self.toggle_theme
        )
        self.theme_toggle_btn.grid(row=row, column=5)

        self.quit_btn_top = ttk.Button(self.scrollable_frame, text="Quit", command=self.root.quit)
        self.quit_btn_top.grid(row=row, column=6, padx=(55, 10))
        row += 1

        # Parameter filter
        ttk.Label(self.scrollable_frame, text="Parameter Set", style="Section.TLabel").grid(
            row=row, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="w", padx=2, pady=(20, 0)
        )
        row += 1

        ttk.Label(self.scrollable_frame, text="Filter:").grid(row=row, column=0, sticky="e")
        self.filter_var = tk.StringVar()
        filter_entry = ttk.Entry(self.scrollable_frame, textvariable=self.filter_var, width=24)
        filter_entry.grid(row=row, column=1, sticky="w", padx=5, pady=(0, 8))
        ToolTip(filter_entry, "Type to filter parameter names (substring match).")
        self.filter_var.trace_add("write", lambda *_: self.apply_filter())
        row += 1

        # Parameters (grouped)
        self.params_frame = ttk.Frame(self.scrollable_frame)
        self.params_frame.grid(row=row, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="ew")
        row += 1

        def tooltip_for(p: ParamDef) -> str:
            parts: List[str] = []
            if p.help:
                parts.append(p.help)
            if p.units:
                parts.append(f"Units: {p.units}")
            if p.default is not None and p.default != "":
                parts.append(f"Default: {p.default}")
            if p.kind == "enum" and p.choices:
                # keep compact
                opts = ", ".join([f"{c.value}={c.label}" for c in p.choices[:6]])
                if len(p.choices) > 6:
                    opts += ", …"
                parts.append(f"Choices: {opts}")
            return "\n".join(parts) if parts else f"Enter value for {p.name}"

        # Build grouped parameter layout in a *single* grid so columns align across groups.
        group_order: List[str] = []
        group_to_params: Dict[str, List[ParamDef]] = {}
        for p in self.params:
            g = p.group or "Ungrouped"
            if g not in group_to_params:
                group_to_params[g] = []
                group_order.append(g)
            group_to_params[g].append(p)

        # Configure columns (label, widget) pairs.
        total_param_cols = UI_PARAM_COLS_PER_ROW * 2
        for col in range(total_param_cols):
            # keep compact; no stretching
            self.params_frame.columnconfigure(col, weight=0)

        cur_row = 0
        for g in group_order:
            # Group header
            hdr = ttk.Label(self.params_frame, text=g, style="GroupHeader.TLabel")
            hdr.grid(row=cur_row, column=0, columnspan=total_param_cols, sticky="w", padx=4, pady=(8, 2))
            cur_row += 1

            num_columns = UI_PARAM_COLS_PER_ROW  # params per row
            params_in_group = group_to_params[g]
            for i, p in enumerate(params_in_group):
                r = cur_row + (i // num_columns)
                c = (i % num_columns) * 2

                label_text = p.name if not p.units else f"{p.name} [{p.units}]"
                lbl = ttk.Label(self.params_frame, text=label_text, style="Param.TLabel")
                lbl.grid(row=r, column=c, padx=(2, 2), pady=4, sticky="e")
                ToolTip(lbl, tooltip_for(p))

                # Widget
                if p.kind == "enum":
                    values = [ch.label for ch in (p.choices or [])]
                    cb = ttk.Combobox(self.params_frame, values=values, width=UI_COMBO_WIDTH, state="readonly")
                    cb.grid(row=r, column=c + 1, padx=(2, 6), pady=4, sticky="w")
                    enum_label_to_value = {ch.label: int(ch.value) for ch in (p.choices or [])}
                    enum_value_to_label = {int(ch.value): ch.label for ch in (p.choices or [])}
                    pw = ParamWidget(
                        pdef=p,
                        label=lbl,
                        widget=cb,
                        var=None,
                        enum_label_to_value=enum_label_to_value,
                        enum_value_to_label=enum_value_to_label,
                    )

                elif p.ui == "checkbox" and p.kind == "int":
                    v = tk.IntVar(value=0)
                    chk = ttk.Checkbutton(self.params_frame, variable=v)
                    chk.grid(row=r, column=c + 1, padx=(2, 6), pady=4, sticky="w")
                    pw = ParamWidget(pdef=p, label=lbl, widget=chk, var=v)

                else:
                    # Entry-like
                    validate = "key" if p.kind in ("int", "float") else "none"
                    vcmd = self.vcmd if validate == "key" else None
                    ent = ttk.Entry(self.params_frame, width=UI_ENTRY_WIDTH, validate=validate, validatecommand=vcmd)
                    ent.grid(row=r, column=c + 1, padx=(2, 6), pady=4, sticky="w")
                    pw = ParamWidget(pdef=p, label=lbl, widget=ent)

                # Apply defaults (if any)
                if p.default is not None and p.default != "":
                    pw.set_value(p.default)

                self.param_widgets.append(pw)
                self.param_rows.append((lbl, pw.widget))

            # advance to next group's start row
            rows_used = (len(params_in_group) + num_columns - 1) // num_columns
            cur_row += rows_used

        ttk.Separator(self.scrollable_frame, orient="horizontal").grid(
            row=row, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="ew", pady=10
        )
        row += 1

        # Sweep Configuration
        ttk.Label(self.scrollable_frame, text="Sweep Configuration", style="Section.TLabel").grid(
            row=row, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="w", padx=2, pady=(10, 0)
        )
        row += 1

        self.sweeps = []
        for si in range(1, 4):
            ttk.Label(self.scrollable_frame, text=f"Param #{si}:").grid(row=row, column=0, sticky="e")
            cb = ttk.Combobox(self.scrollable_frame, values=self.labels, width=12)
            cb.grid(row=row, column=1)

            ttk.Label(self.scrollable_frame, text="Points:").grid(row=row, column=2, sticky="e")
            pts = ttk.Entry(self.scrollable_frame, width=10)
            pts.grid(row=row, column=3)

            ttk.Label(self.scrollable_frame, text="Init:").grid(row=row, column=4, sticky="e")
            init = ttk.Entry(self.scrollable_frame, width=10)
            init.grid(row=row, column=5)

            ttk.Label(self.scrollable_frame, text="Final:").grid(row=row, column=6, sticky="e")
            fin = ttk.Entry(self.scrollable_frame, width=10)
            fin.grid(row=row, column=7)

            sw = SweepWidgets(
                label=ttk.Label(self.scrollable_frame, text=f"Sweep {si}"),
                param=cb, points=pts, init=init, final=fin
            )
            sw.attach_tooltips()
            self.sweeps.append(sw)
            row += 1

        # Actions
        ttk.Label(self.scrollable_frame, text="Actions", style="Section.TLabel").grid(
            row=row, column=0, columnspan=UI_TOTAL_COLUMNS, sticky="w", padx=2, pady=(20, 0)
        )
        row += 1

        ttk.Button(self.scrollable_frame, text="Export Parameters", command=self.export_parameters).grid(row=row, column=0, pady=20)
        ttk.Button(self.scrollable_frame, text="Load Parameters", command=self.load_parameters).grid(row=row, column=1, pady=20)
        ttk.Button(self.scrollable_frame, text="Preview Sweep(s)", command=self.preview_sweeps).grid(row=row, column=2, pady=20)
        ttk.Button(self.scrollable_frame, text="Clear All Fields", command=self.clear_all_fields).grid(row=row, column=3, pady=20)
        ttk.Button(self.scrollable_frame, text="Reset to Scenario", command=self.reset_to_scenario_defaults).grid(row=row, column=4, pady=20)
        ttk.Button(self.scrollable_frame, text="Quit", command=self.root.quit).grid(row=row, column=5, pady=20)

    def apply_filter(self):
        q = self.filter_var.get().strip().lower()
        for (lbl, entry) in self.param_rows:
            name = lbl.cget("text").lower()
            if (not q) or (q in name):
                lbl.grid()
                entry.grid()
            else:
                lbl.grid_remove()
                entry.grid_remove()

    # -----------------------------
    # Scenario / parameters
    # -----------------------------

    def clear_all_fields(self):
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

        values = self.scenarios[scenario]

        # New scenarios.yaml format: dict keyed by parameter name
        if isinstance(values, dict):
            for pw in self.param_widgets:
                name = pw.pdef.name
                if name in values:
                    v = values.get(name)
                else:
                    # fall back to declared default (if any)
                    v = pw.pdef.default if pw.pdef.default is not None else ""
                pw.set_value(v)
            return

        # Backward compatibility: old format (list)
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

        messagebox.showerror(
            "Scenario error",
            f"Scenario '{scenario}' has unsupported format: {type(values)}"
        )

    def reset_to_scenario_defaults(self):
        self.on_scenario_change()

    # -----------------------------
    # Theme
    # -----------------------------

    def toggle_theme(self):
        self.dark_mode = not self.dark_mode
        sv_ttk.set_theme("dark" if self.dark_mode else "light")
        self.theme_toggle_btn.config(
            text="Switch to Light Mode" if self.dark_mode else "Switch to Dark Mode"
        )
        self.save_config()

    # -----------------------------
    # Sweeps parsing / preview
    # -----------------------------

    def parse_sweeps(self) -> Tuple[List[Tuple[int, np.ndarray, str]], List[str]]:
        """
        Returns:
          sweeps: list of (param_index, values, param_name)
          errors: list of error strings
        """
        errors: List[str] = []
        sweeps: List[Tuple[int, np.ndarray, str]] = []

        seen_params: List[str] = []

        for si, sw in enumerate(self.sweeps, start=1):
            pname = sw.param.get().strip()
            # empty means unused sweep row
            if pname == "":
                continue

            if pname not in self.labels:
                errors.append(f"Sweep #{si}: invalid parameter {pname!r}")
                continue

            pkind = self.param_by_name[pname].kind
            if pkind not in ("int", "float"):
                errors.append(f"Sweep #{si} ({pname}): only int/float parameters can be swept (got {pkind})")
                continue

            if pname in seen_params:
                errors.append(f"Sweep #{si}: parameter {pname!r} is duplicated across sweeps")
                continue
            seen_params.append(pname)

            pts_txt = sw.points.get().strip()
            init_txt = sw.init.get().strip()
            fin_txt = sw.final.get().strip()

            if pts_txt == "" or init_txt == "" or fin_txt == "":
                errors.append(f"Sweep #{si} ({pname}): points/init/final must all be filled")
                continue

            try:
                n = int(float(pts_txt))
                if n < 1:
                    raise ValueError("points must be >= 1")
            except Exception:
                errors.append(f"Sweep #{si} ({pname}): Points must be an integer ≥ 1")
                continue

            try:
                i = float(init_txt)
                f = float(fin_txt)
                vals = np.linspace(i, f, n)
            except Exception:
                errors.append(f"Sweep #{si} ({pname}): Init/Final must be numeric")
                continue

            sweeps.append((self.labels.index(pname), vals, pname))

        return sweeps, errors

    def preview_sweeps(self):
        sweeps, errors = self.parse_sweeps()
        if errors:
            messagebox.showerror("Invalid Sweep(s)", "\n".join(errors))
            return

        if not sweeps:
            messagebox.showinfo("Preview", "No sweeps configured.")
            return

        total_rows = int(np.prod([len(v) for _, v, _ in sweeps]))
        messagebox.showinfo("Sweep Summary", f"Configured sweeps: {len(sweeps)}\nTotal runs: {total_rows}")

        # Show one figure per sweep (simple and clear)
        for _, vals, pname in sweeps:
            plt.figure(figsize=(6, 4))
            plt.plot(range(len(vals)), vals, marker="o")
            plt.title(f"Sweep of '{pname}'")
            plt.xlabel("Index")
            plt.ylabel(pname)
            plt.grid(True)
            plt.tight_layout()
        plt.show()

    # -----------------------------
    # Export / load
    # -----------------------------

    def export_parameters(self):
        # Validate simulation number
        sim_no_txt = self.sim_no_entry.get().strip()
        if not sim_no_txt.isdigit():
            messagebox.showerror("Invalid Input", "Simulation number must be an integer.")
            return
        sim_no = int(sim_no_txt)

        # Validate sweeps
        sweeps, sweep_errors = self.parse_sweeps()
        if sweep_errors:
            messagebox.showerror("Invalid Sweep(s)", "\n".join(sweep_errors))
            return

        # Validate all parameter widgets (collect all errors)
        value_errors: List[str] = []
        formatted_values: List[str] = [""] * len(self.param_widgets)

        for i, pw in enumerate(self.param_widgets):
            v, err = pw.get_export_value()
            if err:
                value_errors.append(err)
            else:
                formatted_values[i] = str(v)

        if value_errors:
            # show up to 50 to avoid giant popups
            msg = "Please fix these value issues:\n\n" + "\n".join(value_errors[:50])
            if len(value_errors) > 50:
                msg += f"\n\n… and {len(value_errors) - 50} more"
            messagebox.showerror("Invalid Parameter Values", msg)
            return

        # Build output directory
        script_dir = Path(__file__).resolve().parent
        sims_dir = script_dir / "Sims"
        sims_dir.mkdir(parents=True, exist_ok=True)
        file_path = sims_dir / f"Sim_{sim_no:02d}.dat"

        if file_path.exists():
            if not messagebox.askyesno("Overwrite?", f"{file_path.name} exists. Overwrite?"):
                return

        # Write file
        try:
            total_rows = int(np.prod([len(vals) for _, vals, _ in sweeps])) if sweeps else 1

            with open(file_path, "w", encoding="utf-8") as f:
                f.write(f"{total_rows}\n{len(self.labels)}\n")
                # Header row: parameter names define the columns
                f.write(" ".join(self.labels) + "\n")

                def write_row(row_vals_formatted: List[str]):
                    f.write(" ".join(row_vals_formatted) + "\n")

                if not sweeps:
                    write_row(formatted_values)
                else:
                    # Recursively write cartesian product of sweeps
                    def recurse(level: int, current: List[str]):
                        if level == len(sweeps):
                            write_row(current)
                            return
                        idx, vals, pname = sweeps[level]
                        kind = self.param_by_name[pname].kind
                        for v in vals:
                            nxt = current.copy()
                            nxt[idx] = format_value(str(v), kind)
                            recurse(level + 1, nxt)

                    recurse(0, formatted_values)

            messagebox.showinfo("Success", f"Parameters saved to:\n{file_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to write file:\n{e}")

    def load_parameters(self):
        file_path = filedialog.askopenfilename(filetypes=[("DAT Files", "*.dat")])
        if not file_path:
            return
        file_path = Path(file_path)

        try:
            lines = [ln.strip() for ln in file_path.read_text(encoding="utf-8").splitlines() if ln.strip()]
            if len(lines) < 3:
                messagebox.showerror("Error", "File format not recognized.")
                return

            try:
                rows = int(lines[0])
                cols = int(lines[1])
            except Exception:
                messagebox.showerror("Error", "Header lines corrupted or missing.")
                return

            if cols != len(self.labels):
                messagebox.showerror("Mismatch", f"File has {cols} parameters, expected {len(self.labels)}.")
                return

            third = lines[2].split()
            has_header = (len(third) == cols) and all(tok in self.labels for tok in third)

            if has_header:
                header_names = third
                if len(lines) < 4:
                    messagebox.showerror("Error", "File has a header but no data rows.")
                    return
                data = lines[3].split()
            else:
                header_names = self.labels[:]  # assume canonical order
                data = third

            if len(data) != cols:
                messagebox.showerror("Mismatch", "Parameter count mismatch in first data row.")
                return

            name_to_val = dict(zip(header_names, data))
            for pw in self.param_widgets:
                val = name_to_val.get(pw.pdef.name, "")
                pw.set_value(val)

            messagebox.showinfo("Loaded", f"Loaded parameters from {file_path}\n(total rows in file: {rows})")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {e}")


if __name__ == "__main__":
    root = tk.Tk()
    app = SimulationGUI(root)
    root.mainloop()
