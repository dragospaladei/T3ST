import os
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from scipy.ndimage import gaussian_filter1d
from tkinter import filedialog, messagebox, ttk


SCAN_PARAMETER_OFFSETS = [19, 15, 32, 34, 26, 27, 42, 44, 46, 47, 48, 50, 53, 63, 66]


# ----------------------------
# Low-level IO / utilities
# ----------------------------

def read_table_flat(path: str) -> np.ndarray:
    with open(path, "r") as f:
        tokens = f.read().split()
    return np.array(tokens, dtype=object)



def safe_name(name: str) -> str:
    return name.replace("_", "")



def as_int(x: Any) -> int:
    return int(round(float(x)))



def read_bin_f64(path: str) -> np.ndarray:
    return np.fromfile(path, dtype=np.float64)



def reshape_or_raise(arr: np.ndarray, shape: Tuple[int, ...], path: str) -> np.ndarray:
    expected = int(np.prod(shape))
    if arr.size != expected:
        raise ValueError(
            f"Size mismatch for {path}: got {arr.size} floats, expected {expected} for shape {shape}"
        )
    return arr.reshape(shape)



def aver(noact: np.ndarray, dtt: np.ndarray, Np: float, pas: float) -> np.ndarray:
    weighted = noact[..., None] * dtt
    ts = weighted.mean(axis=0).mean(axis=0) / Np
    return gaussian_filter1d(ts, sigma=pas, mode="nearest")



def speed1(dtt: np.ndarray, tau: np.ndarray, rhoi: float, R0: float, vth: float) -> np.ndarray:
    num = rhoi * (np.roll(dtt, -1) - np.roll(dtt, +1))
    denom = (R0 * (2.0 * tau[1])) / vth
    return np.column_stack([tau, num / denom])



def difuz1(dtt2: np.ndarray, dtt1: np.ndarray, tau: np.ndarray, rhoi: float, R0: float, vth: float) -> np.ndarray:
    base = dtt2 - dtt1**2
    num = (rhoi**2) * (np.roll(base, -1) - np.roll(base, +1))
    denom = (R0 * (4.0 * tau[1])) / vth
    return np.column_stack([tau, num / denom])



def speed2(dtt: np.ndarray, tau: np.ndarray, rhoi: float, R0: float, vth: float) -> np.ndarray:
    num = rhoi * (dtt - dtt[0])
    denom = (R0 * (tau + 1e-9)) / vth
    return np.column_stack([tau, num / denom])



def difuz2(dtt2: np.ndarray, dtt1: np.ndarray, tau: np.ndarray, rhoi: float, R0: float, vth: float) -> np.ndarray:
    base0 = dtt2[0] - dtt1[0] ** 2
    num = (rhoi**2) * ((dtt2 - dtt1**2) - base0)
    denom = (R0 * (2.0 * tau + 1e-9)) / vth
    return np.column_stack([tau, num / denom])


# ----------------------------
# Data model
# ----------------------------
@dataclass
class RunData:
    pars: Dict[str, Any]
    raw_names: List[str]
    raw_values: List[Any]
    noact: np.ndarray
    Q1: np.ndarray
    VT: np.ndarray
    V0: np.ndarray
    LT0: np.ndarray
    Q2: np.ndarray
    LTT: np.ndarray
    L0T: np.ndarray
    L00: np.ndarray
    X: Optional[np.ndarray] = None
    Y: Optional[np.ndarray] = None
    Z: Optional[np.ndarray] = None
    Vp: Optional[np.ndarray] = None
    Mu: Optional[np.ndarray] = None
    H: Optional[np.ndarray] = None
    q1: Optional[np.ndarray] = None
    q2: Optional[np.ndarray] = None
    q3: Optional[np.ndarray] = None
    Pc: Optional[np.ndarray] = None
    ck1: Optional[np.ndarray] = None
    ck2: Optional[np.ndarray] = None
    ck3: Optional[np.ndarray] = None
    Xallin: Optional[np.ndarray] = None
    Yallin: Optional[np.ndarray] = None
    Xallout: Optional[np.ndarray] = None
    Yallout: Optional[np.ndarray] = None
    Hallout: Optional[np.ndarray] = None
    q1allout: Optional[np.ndarray] = None
    q2allout: Optional[np.ndarray] = None



def load_run(sim_folder: str, run_index: int, garbage: bool = True) -> RunData:
    run_dir = os.path.join(sim_folder, f"Run_{run_index:05d}")

    flat = read_table_flat(os.path.join(run_dir, "file_01.dat"))
    n = flat.size
    half = n // 2
    raw_names = [str(x) for x in flat[:half]]
    values_raw = flat[half:half + half]

    raw_values: List[Any] = []
    for v in values_raw:
        try:
            raw_values.append(float(v))
        except Exception:
            raw_values.append(v)

    safe_names = [safe_name(s) for s in raw_names]
    pars: Dict[str, Any] = dict(zip(safe_names, raw_values))

    Nreal = as_int(pars["Nreal"])
    Nt = as_int(pars["Nt"])
    Nloop = as_int(pars["Nloop"])
    ntraj = as_int(pars["ntraj"])

    noact = reshape_or_raise(
        read_bin_f64(os.path.join(run_dir, "file_24.bin")),
        (Nreal, Nloop),
        os.path.join(run_dir, "file_24.bin"),
    )

    dar11 = reshape_or_raise(
        read_bin_f64(os.path.join(run_dir, "file_11.bin")),
        (Nreal, Nloop, Nt + 1, 4),
        os.path.join(run_dir, "file_11.bin"),
    )
    Q1, VT, V0, LT0 = dar11[..., 0], dar11[..., 1], dar11[..., 2], dar11[..., 3]

    dar12 = reshape_or_raise(
        read_bin_f64(os.path.join(run_dir, "file_12.bin")),
        (Nreal, Nloop, Nt + 1, 4),
        os.path.join(run_dir, "file_12.bin"),
    )
    Q2, LTT, L0T, L00 = dar12[..., 0], dar12[..., 1], dar12[..., 2], dar12[..., 3]

    rd = RunData(
        pars=pars,
        raw_names=raw_names,
        raw_values=raw_values,
        noact=noact,
        Q1=Q1,
        VT=VT,
        V0=V0,
        LT0=LT0,
        Q2=Q2,
        LTT=LTT,
        L0T=L0T,
        L00=L00,
    )

    if garbage:
        shp4 = (Nreal, Nloop, ntraj, Nt + 1)
        rd.X = reshape_or_raise(read_bin_f64(os.path.join(run_dir, "file_02.bin")), shp4, "file_02.bin")
        rd.Y = reshape_or_raise(read_bin_f64(os.path.join(run_dir, "file_03.bin")), shp4, "file_03.bin")
        rd.Z = reshape_or_raise(read_bin_f64(os.path.join(run_dir, "file_04.bin")), shp4, "file_04.bin")
        rd.Vp = reshape_or_raise(read_bin_f64(os.path.join(run_dir, "file_05.bin")), shp4, "file_05.bin")
        rd.Mu = reshape_or_raise(read_bin_f64(os.path.join(run_dir, "file_06.bin")), shp4, "file_06.bin")
        rd.H = reshape_or_raise(read_bin_f64(os.path.join(run_dir, "file_07.bin")), shp4, "file_07.bin")
        rd.q1 = reshape_or_raise(read_bin_f64(os.path.join(run_dir, "file_08.bin")), shp4, "file_08.bin")
        rd.q2 = reshape_or_raise(read_bin_f64(os.path.join(run_dir, "file_09.bin")), shp4, "file_09.bin")
        rd.q3 = reshape_or_raise(read_bin_f64(os.path.join(run_dir, "file_10.bin")), shp4, "file_10.bin")
        rd.Pc = reshape_or_raise(read_bin_f64(os.path.join(run_dir, "file_25.bin")), shp4, "file_25.bin")
        rd.ck1 = reshape_or_raise(read_bin_f64(os.path.join(run_dir, "file_26.bin")), shp4, "file_26.bin")
        rd.ck2 = reshape_or_raise(read_bin_f64(os.path.join(run_dir, "file_27.bin")), shp4, "file_27.bin")
        rd.ck3 = reshape_or_raise(read_bin_f64(os.path.join(run_dir, "file_28.bin")), shp4, "file_28.bin")

        rd.Xallin = read_bin_f64(os.path.join(run_dir, "file_13.bin"))
        rd.Yallin = read_bin_f64(os.path.join(run_dir, "file_14.bin"))
        rd.Xallout = read_bin_f64(os.path.join(run_dir, "file_15.bin"))
        rd.Yallout = read_bin_f64(os.path.join(run_dir, "file_16.bin"))
        rd.Hallout = read_bin_f64(os.path.join(run_dir, "file_17.bin"))
        rd.q1allout = read_bin_f64(os.path.join(run_dir, "file_18.bin"))
        rd.q2allout = read_bin_f64(os.path.join(run_dir, "file_19.bin"))

    return rd


class SimulationDataset:
    def __init__(self, sim_folder: str, run_start: int, run_end: int, garbage: bool = True):
        self.sim_folder = sim_folder
        self.run_start = run_start
        self.run_end = run_end
        self.garbage = garbage
        self.simulation_number = self._infer_simulation_number(sim_folder)

        self.runs: List[RunData] = []
        self._computed = False

        self.t_list: List[np.ndarray] = []
        self.ScriptQ1: List[np.ndarray] = []
        self.ScriptQ2: List[np.ndarray] = []
        self.Pi1: List[np.ndarray] = []
        self.Pi2: List[np.ndarray] = []
        self.ScriptH1: List[np.ndarray] = []
        self.ScriptH2: List[np.ndarray] = []
        self.ScriptVt1: List[np.ndarray] = []
        self.ScriptVt2: List[np.ndarray] = []
        self.DQ1: List[np.ndarray] = []
        self.DQ2: List[np.ndarray] = []
        self.VQ1: List[np.ndarray] = []
        self.VQ2: List[np.ndarray] = []
        self.DQas: Optional[np.ndarray] = None
        self.VQas: Optional[np.ndarray] = None
        self.DQas2: Optional[np.ndarray] = None
        self.VQas2: Optional[np.ndarray] = None

    @staticmethod
    def _infer_simulation_number(sim_folder: str) -> Optional[int]:
        base = os.path.basename(os.path.normpath(sim_folder))
        if base.lower().startswith("sim_"):
            try:
                return int(base.split("_")[-1])
            except ValueError:
                return None
        return None

    def load(self):
        self.runs = [load_run(self.sim_folder, i, garbage=self.garbage) for i in range(self.run_start, self.run_end + 1)]
        self._computed = False

    def get_scan_parameter(self) -> Tuple[Optional[str], Optional[np.ndarray]]:
        if not self.runs or self.simulation_number is None:
            return None, None
        if self.simulation_number < 1 or self.simulation_number > len(SCAN_PARAMETER_OFFSETS):
            return None, None
        offset = SCAN_PARAMETER_OFFSETS[self.simulation_number - 1]
        first = self.runs[0]
        if offset < 1 or offset > len(first.raw_names):
            return None, None
        name = first.raw_names[offset - 1]
        values = np.array([float(r.raw_values[offset - 1]) for r in self.runs], dtype=float)
        return name, values

    def compute_derived(self, pas: float = 1.0, frac: float = 0.7, frac2: float = 0.95):
        if not self.runs:
            raise RuntimeError("No runs loaded")

        tmax = np.array([r.pars["tmax"] for r in self.runs], dtype=float)
        Nt = np.array([as_int(r.pars["Nt"]) for r in self.runs], dtype=int)
        dt = tmax / Nt
        self.t_list = [dt[k] * np.arange(Nt[k] + 1, dtype=float) for k in range(len(self.runs))]

        self.ScriptQ1 = [aver(r.noact, r.Q1, r.pars["Np"], pas) for r in self.runs]
        self.ScriptQ2 = [aver(r.noact, r.Q2, r.pars["Np"], pas) for r in self.runs]
        self.Pi1 = [aver(r.noact, r.VT, r.pars["Np"], pas) for r in self.runs]
        self.Pi2 = [aver(r.noact, r.LTT, r.pars["Np"], pas) for r in self.runs]
        self.ScriptH1 = [aver(r.noact, r.V0, r.pars["Np"], pas) for r in self.runs]
        self.ScriptH2 = [aver(r.noact, r.L0T, r.pars["Np"], pas) for r in self.runs]
        self.ScriptVt1 = [aver(r.noact, r.LT0, r.pars["Np"], pas) for r in self.runs]
        self.ScriptVt2 = [aver(r.noact, r.L00, r.pars["Np"], pas) for r in self.runs]

        self.VQ1, self.DQ1, self.VQ2, self.DQ2 = [], [], [], []
        for k, r in enumerate(self.runs):
            tau = self.t_list[k]
            rhoi = float(r.pars["rhoi"])
            R0 = float(r.pars["R0"])
            vth = float(r.pars["vth"])

            vq1 = speed1(self.ScriptQ1[k], tau, rhoi, R0, vth)
            dq1 = difuz1(self.ScriptQ2[k], self.ScriptQ1[k], tau, rhoi, R0, vth)
            vq2 = speed2(self.ScriptQ1[k], tau, rhoi, R0, vth)
            dq2 = difuz2(self.ScriptQ2[k], self.ScriptQ1[k], tau, rhoi, R0, vth)

            dq1[0, 1] = 0.0
            dq1[Nt[k], 1] = dq1[Nt[k] - 1, 1]
            vq1[0, 1] = vq1[1, 1]

            self.VQ1.append(vq1)
            self.DQ1.append(dq1)
            self.VQ2.append(vq2)
            self.DQ2.append(dq2)

        self.VQas = np.array([self.VQ1[k][int(frac * Nt[k]): int(Nt[k] - 1), 1].mean() for k in range(len(self.runs))])
        self.DQas = np.array([self.DQ1[k][int(frac * Nt[k]): int(Nt[k] - 1), 1].mean() for k in range(len(self.runs))])
        self.VQas2 = np.array([self.VQ2[k][int(frac2 * Nt[k]): int(Nt[k] - 1), 1].mean() for k in range(len(self.runs))])
        self.DQas2 = np.array([self.DQ2[k][int(frac2 * Nt[k]): int(Nt[k] - 1), 1].mean() for k in range(len(self.runs))])
        self._computed = True


# ----------------------------
# Plot recipes
# ----------------------------
def _new_figure(nrows=1, ncols=1, figsize=(7.5, 5.0)):
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    if isinstance(axes, np.ndarray):
        axes_list = axes.ravel().tolist()
    else:
        axes_list = [axes]
    return fig, axes_list



def plot_dq1_all(ds: SimulationDataset, params: Dict[str, Any]):
    fig, (ax,) = _new_figure()
    for k, arr in enumerate(ds.DQ1):
        ax.plot(arr[:, 0], arr[:, 1], label=f"run {ds.run_start + k}")
    ax.set_xlabel("t")
    ax.set_ylabel("DQ1")
    ax.set_title("DQ1 (all runs)")
    ax.legend()
    fig.tight_layout()
    return fig



def plot_vq1_all(ds: SimulationDataset, params: Dict[str, Any]):
    fig, (ax,) = _new_figure()
    for k, arr in enumerate(ds.VQ1):
        ax.plot(arr[:, 0], arr[:, 1], label=f"run {ds.run_start + k}")
    ax.set_xlabel("t [R0/vth]")
    ax.set_ylabel("V [m^2/s]")
    ax.set_title("VQ1 (all runs)")
    ax.legend()
    fig.tight_layout()
    return fig



def plot_scan_vqas(ds: SimulationDataset, params: Dict[str, Any]):
    name, values = ds.get_scan_parameter()
    if name is None or values is None or ds.VQas is None:
        raise RuntimeError("Could not infer the varied scan parameter for this simulation folder.")
    fig, (ax,) = _new_figure()
    ax.plot(values, ds.VQas, marker="o")
    ax.set_xlabel(name)
    ax.set_ylabel("V [m^2/s]")
    ax.set_title("VQas vs scan parameter")
    fig.tight_layout()
    return fig



def plot_scan_dqas(ds: SimulationDataset, params: Dict[str, Any]):
    name, values = ds.get_scan_parameter()
    if name is None or values is None or ds.DQas is None:
        raise RuntimeError("Could not infer the varied scan parameter for this simulation folder.")
    fig, (ax,) = _new_figure()
    ax.plot(values, ds.DQas, marker="o")
    ax.set_xlabel(name)
    ax.set_ylabel("D [m^2/s]")
    ax.set_title("DQas vs scan parameter")
    fig.tight_layout()
    return fig



def plot_xy_in_out(ds: SimulationDataset, params: Dict[str, Any]):
    as_idx = int(params.get("as", 1)) - 1
    nmax = int(params.get("nmax", 10000))
    r = ds.runs[as_idx]
    if r.Xallout is None or r.Yallout is None or r.Xallin is None or r.Yallin is None:
        raise RuntimeError("Xallin/out not loaded (set garbage=True)")
    fig, (ax,) = _new_figure(figsize=(7, 6))
    ax.scatter(r.Xallout[:nmax], r.Yallout[:nmax], s=4, label="out")
    ax.scatter(r.Xallin[:nmax], r.Yallin[:nmax], s=4, label="in")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title(f"X-Y in/out (run as={as_idx + 1}, nmax={nmax})")
    ax.legend()
    fig.tight_layout()
    return fig



def plot_q_variance_and_mean(ds: SimulationDataset, params: Dict[str, Any]):
    as_idx = int(params.get("as", 1)) - 1
    fig, axes = _new_figure(nrows=1, ncols=2, figsize=(12, 4.5))
    ax1, ax2 = axes
    q1 = ds.ScriptQ1[as_idx]
    q2 = ds.ScriptQ2[as_idx]
    ax1.plot(q2 - q1**2)
    ax1.set_title("ScriptQ2 - ScriptQ1^2")
    ax2.plot(q1)
    ax2.set_title("ScriptQ1")
    fig.tight_layout()
    return fig



def plot_dqas(ds: SimulationDataset, params: Dict[str, Any]):
    if ds.DQas is None:
        raise RuntimeError("Derived quantities not computed")
    fig, (ax,) = _new_figure()
    x = np.arange(ds.run_start, ds.run_end + 1)
    ax.plot(x, ds.DQas, marker="o")
    ax.set_xlabel("run index")
    ax.set_ylabel("DQas")
    ax.set_title("DQas vs run")
    fig.tight_layout()
    return fig



def plot_phase_space_triplet(ds: SimulationDataset, params: Dict[str, Any]):
    os_idx = int(params.get("os", 1)) - 1
    whc = int(params.get("whc", 1)) - 1
    loop = int(params.get("loop", 1)) - 1
    traj = int(params.get("traj", 1)) - 1

    r = ds.runs[os_idx]
    for name in ("q1", "q2", "q3", "X", "Y"):
        if getattr(r, name) is None:
            raise RuntimeError(f"{name} not loaded (set garbage=True)")

    q1 = r.q1[:, loop, traj, whc]
    q2 = r.q2[:, loop, traj, whc]
    q3 = r.q3[:, loop, traj, whc]
    X = r.X[:, loop, traj, whc]
    Y = r.Y[:, loop, traj, whc]

    fig, axes = _new_figure(nrows=1, ncols=3, figsize=(14, 4.5))
    ax1, ax2, ax3 = axes
    ax1.plot(q1, q2)
    ax1.set_xlabel("x (q1)")
    ax1.set_ylabel("y (q2)")
    ax2.plot(X, Y)
    ax2.set_xlabel("R (X)")
    ax2.set_ylabel("Z (Y)")
    ax3.plot(q2, q3)
    ax3.set_xlabel("z (q2)")
    ax3.set_ylabel("y (q3)")
    fig.tight_layout()
    return fig



def plot_compare_two_runs_q123(ds: SimulationDataset, params: Dict[str, Any]):
    as1 = int(params.get("as1", 1)) - 1
    as2 = int(params.get("as2", 2)) - 1
    loop = int(params.get("loop", 1)) - 1
    traj = int(params.get("traj", 1)) - 1
    pidx = int(params.get("particle", 1)) - 1

    r1 = ds.runs[as1]
    r2 = ds.runs[as2]
    for r in (r1, r2):
        for name in ("q1", "q2", "q3"):
            if getattr(r, name) is None:
                raise RuntimeError(f"{name} not loaded (set garbage=True)")

    fig, axes = _new_figure(nrows=1, ncols=3, figsize=(14, 4.5))
    ax1, ax2, ax3 = axes
    ax1.plot(r1.q1[pidx, loop, traj, :], label=f"as1={as1 + 1}")
    ax1.plot(r2.q1[pidx, loop, traj, :], label=f"as2={as2 + 1}")
    ax2.plot(r1.q2[pidx, loop, traj, :], label=f"as1={as1 + 1}")
    ax2.plot(-r2.q2[pidx, loop, traj, :], label=f"-as2={as2 + 1}")
    ax3.plot(r1.q3[pidx, loop, traj, :], label=f"as1={as1 + 1}")
    ax3.plot(-r2.q3[pidx, loop, traj, :], label=f"-as2={as2 + 1}")
    for ax, title in zip(axes, ["q1 time series", "q2 time series", "q3 time series"]):
        ax.set_title(title)
        ax.legend()
    fig.tight_layout()
    return fig



def plot_h_and_pc_changes(ds: SimulationDataset, params: Dict[str, Any]):
    as_idx = int(params.get("as", 1)) - 1
    loop = int(params.get("loop", 1)) - 1
    traj = int(params.get("traj", 1)) - 1
    ncurves = int(params.get("ncurves", 5))
    scale_h = float(params.get("scale_h", 1e2))
    scale_pc = float(params.get("scale_pc", 1e4))

    r = ds.runs[as_idx]
    if r.H is None or r.Pc is None:
        raise RuntimeError("H/Pc not loaded (set garbage=True)")

    pmax = min(ncurves, r.H.shape[0])
    fig, axes = _new_figure(nrows=1, ncols=3, figsize=(16, 4.5))
    a1, a2, a3 = axes
    for o in range(pmax):
        h = r.H[o, loop, traj, :]
        pc = r.Pc[o, loop, traj, :]
        a1.plot(h, label=f"p{o + 1}")
        a2.plot((h / (h[0] + 1e-300) - 1.0) * scale_h, label=f"p{o + 1}")
        a3.plot((pc / (pc[0] + 1e-300) - 1.0) * scale_pc, label=f"p{o + 1}")
    a1.set_title("H")
    a2.set_title(f"(H/H0-1)*{scale_h:g}")
    a3.set_title(f"(Pc/Pc0-1)*{scale_pc:g}")
    for ax in axes:
        ax.legend()
    fig.tight_layout()
    return fig


PLOT_RECIPES: Dict[str, Tuple[Callable[[SimulationDataset, Dict[str, Any]], plt.Figure], bool]] = {
    "VQas vs scan parameter": (plot_scan_vqas, True),
    "DQas vs scan parameter": (plot_scan_dqas, True),
    "VQ1 all runs": (plot_vq1_all, True),
    "DQ1 all runs": (plot_dq1_all, True),
    "X-Y in/out scatter": (plot_xy_in_out, False),
    "Q variance & mean": (plot_q_variance_and_mean, True),
    "DQas vs run": (plot_dqas, True),
    "Phase/trajectory triplet": (plot_phase_space_triplet, False),
    "Compare two runs q1/q2/q3": (plot_compare_two_runs_q123, False),
    "H & Pc time series + relative change": (plot_h_and_pc_changes, False),
}


# ----------------------------
# Tkinter GUI
# ----------------------------
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Simulation Plotter")
        self.geometry("1250x760")
        self.ds: Optional[SimulationDataset] = None
        self.current_fig: Optional[plt.Figure] = None
        self._build_ui()

    def _build_ui(self):
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        left = ttk.Frame(self, padding=10)
        left.grid(row=0, column=0, sticky="ns")

        right = ttk.Frame(self, padding=5)
        right.grid(row=0, column=1, sticky="nsew")
        right.rowconfigure(0, weight=1)
        right.columnconfigure(0, weight=1)

        ttk.Label(left, text="Data folder (Sim_XX)", font=("Segoe UI", 10, "bold")).grid(row=0, column=0, sticky="w")
        self.var_sim_folder = tk.StringVar(value="")
        ttk.Entry(left, textvariable=self.var_sim_folder, width=42).grid(row=1, column=0, sticky="we", pady=(2, 4))
        ttk.Button(left, text="Browse...", command=self._browse_sim_folder).grid(row=2, column=0, sticky="we", pady=(0, 10))

        frm_runs = ttk.LabelFrame(left, text="Run range", padding=8)
        frm_runs.grid(row=3, column=0, sticky="we", pady=(0, 10))
        self.var_run_start = tk.IntVar(value=1)
        self.var_run_end = tk.IntVar(value=20)
        ttk.Label(frm_runs, text="Start").grid(row=0, column=0, sticky="w")
        ttk.Entry(frm_runs, textvariable=self.var_run_start, width=8).grid(row=0, column=1, sticky="w", padx=(6, 0))
        ttk.Label(frm_runs, text="End").grid(row=1, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(frm_runs, textvariable=self.var_run_end, width=8).grid(row=1, column=1, sticky="w", padx=(6, 0), pady=(6, 0))

        frm_opts = ttk.LabelFrame(left, text="Load/compute options", padding=8)
        frm_opts.grid(row=4, column=0, sticky="we", pady=(0, 10))
        self.var_garbage = tk.BooleanVar(value=True)
        ttk.Checkbutton(frm_opts, text="Load big arrays (garbage=True)", variable=self.var_garbage).grid(row=0, column=0, sticky="w")
        self.var_pas = tk.DoubleVar(value=1.0)
        ttk.Label(frm_opts, text="Smoothing pas").grid(row=1, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(frm_opts, textvariable=self.var_pas, width=10).grid(row=2, column=0, sticky="w")

        ttk.Button(left, text="Load data", command=self._load_data).grid(row=5, column=0, sticky="we", pady=(0, 8))
        ttk.Separator(left).grid(row=6, column=0, sticky="we", pady=10)

        ttk.Label(left, text="Plot recipe", font=("Segoe UI", 10, "bold")).grid(row=7, column=0, sticky="w")
        self.var_recipe = tk.StringVar(value=list(PLOT_RECIPES.keys())[0])
        ttk.Combobox(left, textvariable=self.var_recipe, values=list(PLOT_RECIPES.keys()), state="readonly", width=42).grid(row=8, column=0, sticky="we", pady=(2, 10))

        frm_p = ttk.LabelFrame(left, text="Plot parameters", padding=8)
        frm_p.grid(row=9, column=0, sticky="we")

        def add_param(row, label, var, width=10):
            ttk.Label(frm_p, text=label).grid(row=row, column=0, sticky="w")
            ttk.Entry(frm_p, textvariable=var, width=width).grid(row=row, column=1, sticky="w", padx=(6, 0), pady=2)

        self.var_as = tk.IntVar(value=1)
        self.var_os = tk.IntVar(value=1)
        self.var_nmax = tk.IntVar(value=10000)
        self.var_whc = tk.IntVar(value=5)
        self.var_loop = tk.IntVar(value=1)
        self.var_traj = tk.IntVar(value=1)
        self.var_as1 = tk.IntVar(value=1)
        self.var_as2 = tk.IntVar(value=2)
        self.var_particle = tk.IntVar(value=1)
        self.var_ncurves = tk.IntVar(value=5)

        add_param(0, "as (run idx)", self.var_as)
        add_param(1, "os (run idx)", self.var_os)
        add_param(2, "nmax", self.var_nmax)
        add_param(3, "whc (time idx)", self.var_whc)
        add_param(4, "loop", self.var_loop)
        add_param(5, "traj", self.var_traj)
        add_param(6, "as1", self.var_as1)
        add_param(7, "as2", self.var_as2)
        add_param(8, "particle", self.var_particle)
        add_param(9, "ncurves", self.var_ncurves)

        ttk.Button(left, text="Plot", command=self._plot).grid(row=10, column=0, sticky="we", pady=(10, 6))
        ttk.Button(left, text="Save plot...", command=self._save_plot).grid(row=11, column=0, sticky="we")

        plot_frame = ttk.Frame(right)
        plot_frame.grid(row=0, column=0, sticky="nsew")
        plot_frame.rowconfigure(1, weight=1)
        plot_frame.columnconfigure(0, weight=1)

        self.toolbar_frame = ttk.Frame(plot_frame)
        self.toolbar_frame.grid(row=0, column=0, sticky="ew")
        self.canvas_frame = ttk.Frame(plot_frame)
        self.canvas_frame.grid(row=1, column=0, sticky="nsew")

        self.fig = plt.Figure(figsize=(7.5, 5.5), dpi=100)
        ax = self.fig.add_subplot(111)
        ax.text(0.5, 0.5, "Load data, then Plot", ha="center", va="center")
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbar_frame)
        self.toolbar.update()

    def _browse_sim_folder(self):
        folder = filedialog.askdirectory(title="Select Sim_XX folder (contains Run_00001, ...)")
        if folder:
            self.var_sim_folder.set(folder)

    def _load_data(self):
        sim_folder = self.var_sim_folder.get().strip()
        if not sim_folder or not os.path.isdir(sim_folder):
            messagebox.showerror("Error", "Please select a valid Sim_XX folder.")
            return
        run_start = int(self.var_run_start.get())
        run_end = int(self.var_run_end.get())
        if run_end < run_start:
            messagebox.showerror("Error", "Run end must be >= run start.")
            return
        try:
            ds = SimulationDataset(sim_folder, run_start, run_end, garbage=bool(self.var_garbage.get()))
            ds.load()
            ds.compute_derived(pas=float(self.var_pas.get()))
            self.ds = ds
            sim_num_text = f" (simulation {ds.simulation_number})" if ds.simulation_number is not None else ""
            messagebox.showinfo("Loaded", f"Loaded runs {run_start}..{run_end}{sim_num_text} from:\n{sim_folder}")
        except Exception as e:
            messagebox.showerror("Load failed", str(e))
            self.ds = None

    def _gather_plot_params(self) -> Dict[str, Any]:
        return {
            "as": self.var_as.get(),
            "os": self.var_os.get(),
            "nmax": self.var_nmax.get(),
            "whc": self.var_whc.get(),
            "loop": self.var_loop.get(),
            "traj": self.var_traj.get(),
            "as1": self.var_as1.get(),
            "as2": self.var_as2.get(),
            "particle": self.var_particle.get(),
            "ncurves": self.var_ncurves.get(),
        }

    def _plot(self):
        if self.ds is None:
            messagebox.showerror("Error", "Load data first.")
            return
        recipe = self.var_recipe.get()
        plot_fn, needs_derived = PLOT_RECIPES[recipe]
        if needs_derived:
            try:
                self.ds.compute_derived(pas=float(self.var_pas.get()))
            except Exception as e:
                messagebox.showerror("Compute failed", str(e))
                return
        try:
            fig = plot_fn(self.ds, self._gather_plot_params())
            self._draw_figure(fig)
            self.current_fig = fig
        except Exception as e:
            messagebox.showerror("Plot failed", str(e))

    def _draw_figure(self, fig: plt.Figure):
        for w in self.canvas_frame.winfo_children():
            w.destroy()
        for w in self.toolbar_frame.winfo_children():
            w.destroy()
        self.canvas = FigureCanvasTkAgg(fig, master=self.canvas_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbar_frame)
        self.toolbar.update()
        self.canvas.draw()

    def _save_plot(self):
        if self.current_fig is None:
            messagebox.showerror("Error", "Nothing to save yet. Make a plot first.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            self.current_fig.savefig(path, dpi=200, bbox_inches="tight")
            messagebox.showinfo("Saved", f"Saved:\n{path}")
        except Exception as e:
            messagebox.showerror("Save failed", str(e))


def main():
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()

