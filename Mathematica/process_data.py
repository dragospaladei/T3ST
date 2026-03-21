import os
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d


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


@dataclass
class SimulationResult:
    simulation: int
    sim_folder: str
    runs: List[RunData]
    t_list: List[np.ndarray]
    ScriptQ1: List[np.ndarray]
    ScriptQ2: List[np.ndarray]
    Pi1: List[np.ndarray]
    Pi2: List[np.ndarray]
    ScriptH1: List[np.ndarray]
    ScriptH2: List[np.ndarray]
    ScriptVt1: List[np.ndarray]
    ScriptVt2: List[np.ndarray]
    VQ1: List[np.ndarray]
    DQ1: List[np.ndarray]
    VQ2: List[np.ndarray]
    DQ2: List[np.ndarray]
    VQas: np.ndarray
    DQas: np.ndarray
    VQas2: np.ndarray
    DQas2: np.ndarray
    shear: np.ndarray
    nub0: np.ndarray
    nub: np.ndarray
    Gamma_gyro: np.ndarray
    scan_parameter_name: Optional[str]
    scan_parameter_values: Optional[np.ndarray]


# ----------------------------
# Loading and processing
# ----------------------------

def load_run(sim_folder: str, run_index: int, garbage: bool = True) -> RunData:
    run_dir = os.path.join(sim_folder, f"Run_{run_index:05d}")
    print("Loading:", run_dir)

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
    Np = as_int(pars["Np"])
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



def get_scan_parameter(simulation: int, runs: List[RunData]) -> Tuple[Optional[str], Optional[np.ndarray]]:
    if not runs:
        return None, None
    if simulation < 1 or simulation > len(SCAN_PARAMETER_OFFSETS):
        return None, None

    offset = SCAN_PARAMETER_OFFSETS[simulation - 1]  # Mathematica is 1-based
    first = runs[0]
    if offset < 1 or offset > len(first.raw_names):
        return None, None

    name = first.raw_names[offset - 1]
    values = []
    for run in runs:
        if offset <= len(run.raw_values):
            values.append(float(run.raw_values[offset - 1]))
        else:
            values.append(np.nan)
    return name, np.array(values, dtype=float)



def process_simulation(sim_folder: str, simulation: int, run_start: int, run_end: int, garbage: bool = True, pas: float = 1.0) -> SimulationResult:
    runs = [load_run(sim_folder, i, garbage=garbage) for i in range(run_start, run_end + 1)]
    nsim2 = len(runs)

    tmax = np.array([r.pars["tmax"] for r in runs], dtype=float)
    Nt = np.array([as_int(r.pars["Nt"]) for r in runs], dtype=int)
    dt = tmax / Nt
    t_list = [dt[k] * np.arange(Nt[k] + 1, dtype=float) for k in range(nsim2)]

    ScriptQ1 = [aver(r.noact, r.Q1, r.pars["Np"], pas) for r in runs]
    ScriptQ2 = [aver(r.noact, r.Q2, r.pars["Np"], pas) for r in runs]
    Pi1 = [aver(r.noact, r.VT, r.pars["Np"], pas) for r in runs]
    Pi2 = [aver(r.noact, r.LTT, r.pars["Np"], pas) for r in runs]
    ScriptH1 = [aver(r.noact, r.V0, r.pars["Np"], pas) for r in runs]
    ScriptH2 = [aver(r.noact, r.L0T, r.pars["Np"], pas) for r in runs]
    ScriptVt1 = [aver(r.noact, r.LT0, r.pars["Np"], pas) for r in runs]
    ScriptVt2 = [aver(r.noact, r.L00, r.pars["Np"], pas) for r in runs]

    VQ1: List[np.ndarray] = []
    DQ1: List[np.ndarray] = []
    VQ2: List[np.ndarray] = []
    DQ2: List[np.ndarray] = []

    for k, r in enumerate(runs):
        tau = t_list[k]
        rhoi = float(r.pars["rhoi"])
        R0 = float(r.pars["R0"])
        vth = float(r.pars["vth"])

        vq1 = speed1(ScriptQ1[k], tau, rhoi, R0, vth)
        dq1 = difuz1(ScriptQ2[k], ScriptQ1[k], tau, rhoi, R0, vth)
        vq2 = speed2(ScriptQ1[k], tau, rhoi, R0, vth)
        dq2 = difuz2(ScriptQ2[k], ScriptQ1[k], tau, rhoi, R0, vth)

        dq1[0, 1] = 0.0
        dq1[Nt[k], 1] = dq1[Nt[k] - 1, 1]
        vq1[0, 1] = vq1[1, 1]

        VQ1.append(vq1)
        DQ1.append(dq1)
        VQ2.append(vq2)
        DQ2.append(dq2)

    r00 = np.array([r.pars["r00"] for r in runs], dtype=float)
    s1 = np.array([r.pars["s1"] for r in runs], dtype=float)
    s3 = np.array([r.pars["s3"] for r in runs], dtype=float)
    a0 = np.array([r.pars["a0"] for r in runs], dtype=float)
    shear = (2.0 * r00**2 * s3) / (a0**2 * s1 + r00**2 * s3)

    eps0 = 8.85e-12
    logL = 17.0
    e = 1.6e-19
    mi = 1.6e-27

    ndens = np.array([r.pars["ndens"] for r in runs], dtype=float)
    Zeff = np.array([r.pars["Zeff"] for r in runs], dtype=float)
    Zw = np.array([r.pars["Zw"] for r in runs], dtype=float)
    Aw = np.array([r.pars["Aw"] for r in runs], dtype=float)
    Aeff = np.array([r.pars["Aeff"] for r in runs], dtype=float)
    vth = np.array([r.pars["vth"] for r in runs], dtype=float)
    R0 = np.array([r.pars["R0"] for r in runs], dtype=float)
    q00 = np.array([r.pars["q00"] for r in runs], dtype=float)
    rhoi = np.array([r.pars["rhoi"] for r in runs], dtype=float)

    nub0 = (
        (4.0 * logL * ndens * ((Zeff * Zw * e**2) / (eps0 * mi * Aw))**2 * (np.sqrt(Aeff / 2.0) ** 3))
        / ((8.0 * np.pi) * vth**3)
    )
    nub = (
        (4.0 * logL * ndens * ((Zeff * Zw * e**2) / (eps0 * mi * Aw))**2 * (np.sqrt(Aeff / 2.0) ** 3))
        * ((r00 * R0) * q00 * np.sqrt(2.0))
        / (((8.0 * np.pi) * vth**3) * (vth * np.sqrt(r00)))
    )

    frac = 0.7
    frac2 = 0.95
    VQas = np.array([VQ1[k][int(frac * Nt[k]): int(Nt[k] - 1), 1].mean() for k in range(nsim2)])
    DQas = np.array([DQ1[k][int(frac * Nt[k]): int(Nt[k] - 1), 1].mean() for k in range(nsim2)])
    VQas2 = np.array([VQ2[k][int(frac2 * Nt[k]): int(Nt[k] - 1), 1].mean() for k in range(nsim2)])
    DQas2 = np.array([DQ2[k][int(frac2 * Nt[k]): int(Nt[k] - 1), 1].mean() for k in range(nsim2)])
    Gamma_gyro = (rhoi / (a0 * R0)) ** 2 * vth * Aw

    scan_parameter_name, scan_parameter_values = get_scan_parameter(simulation, runs)

    return SimulationResult(
        simulation=simulation,
        sim_folder=sim_folder,
        runs=runs,
        t_list=t_list,
        ScriptQ1=ScriptQ1,
        ScriptQ2=ScriptQ2,
        Pi1=Pi1,
        Pi2=Pi2,
        ScriptH1=ScriptH1,
        ScriptH2=ScriptH2,
        ScriptVt1=ScriptVt1,
        ScriptVt2=ScriptVt2,
        VQ1=VQ1,
        DQ1=DQ1,
        VQ2=VQ2,
        DQ2=DQ2,
        VQas=VQas,
        DQas=DQas,
        VQas2=VQas2,
        DQas2=DQas2,
        shear=shear,
        nub0=nub0,
        nub=nub,
        Gamma_gyro=Gamma_gyro,
        scan_parameter_name=scan_parameter_name,
        scan_parameter_values=scan_parameter_values,
    )


# ----------------------------
# Plot/export helpers
# ----------------------------

def plot_multi_curves(curves: List[np.ndarray], title: str, xlabel: str, ylabel: str) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    for k, arr in enumerate(curves):
        arr = np.asarray(arr)
        ax.plot(arr[:, 0], arr[:, 1], label=f"run {k + 1}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    return fig



def plot_scan(values: np.ndarray, y: np.ndarray, xlabel: str, ylabel: str, title: str) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.plot(values, y, marker="o")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    return fig



def plot_in_out(run: RunData, as_index: int, nmax: int = 10_000) -> plt.Figure:
    if run.Xallout is None or run.Yallout is None or run.Xallin is None or run.Yallin is None:
        raise RuntimeError("Xallin/Xallout data not loaded. Set garbage=True.")
    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    ax.scatter(run.Xallout[:nmax], run.Yallout[:nmax], s=4, label="out")
    ax.scatter(run.Xallin[:nmax], run.Yallin[:nmax], s=4, label="in")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title(f"In/Out scatter for run {as_index}")
    ax.legend()
    fig.tight_layout()
    return fig



def export_default_figures(result: SimulationResult, nmax: int = 10_000) -> None:
    folder = result.sim_folder

    if result.scan_parameter_name is not None and result.scan_parameter_values is not None:
        fig = plot_scan(
            result.scan_parameter_values,
            result.VQas,
            result.scan_parameter_name,
            "V [m^2/s]",
            f"Simulation {result.simulation}: VQas scan",
        )
        fig.savefig(os.path.join(folder, "fig_1_a.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)

        fig = plot_scan(
            result.scan_parameter_values,
            result.DQas,
            result.scan_parameter_name,
            "D [m^2/s]",
            f"Simulation {result.simulation}: DQas scan",
        )
        fig.savefig(os.path.join(folder, "fig_1_b.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)

    fig = plot_multi_curves(result.VQ1, "VQ1", "t [R0/vth]", "V [m^2/s]")
    fig.savefig(os.path.join(folder, "fig_2_a.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig = plot_multi_curves(result.DQ1, "DQ1", "t [R0/vth]", "D [m^2/s]")
    fig.savefig(os.path.join(folder, "fig_2_b.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    for as_index, run in enumerate(result.runs, start=1):
        if run.Xallin is None:
            continue
        fig = plot_in_out(run, as_index=as_index, nmax=nmax)
        fig.savefig(os.path.join(folder, f"fig_3_{as_index}.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)


# ----------------------------
# Main script
# ----------------------------

def main() -> None:
    sims = (1, 15)
    run_range = (1, 20)
    garbage = True
    pas = 1.0
    nmax = 10_000

    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(script_dir)
    data_root = os.path.join(repo_root, "data")

    all_results: List[SimulationResult] = []
    for simulation in range(sims[0], sims[1] + 1):
        sim_folder = os.path.join(data_root, f"Sim_{simulation:02d}")
        if not os.path.isdir(sim_folder):
            print(f"Skipping missing folder: {sim_folder}")
            continue

        print(f"\n=== Processing simulation {simulation:02d} ===")
        result = process_simulation(
            sim_folder=sim_folder,
            simulation=simulation,
            run_start=run_range[0],
            run_end=run_range[1],
            garbage=garbage,
            pas=pas,
        )
        export_default_figures(result, nmax=nmax)
        all_results.append(result)

        print("shear:", result.shear)
        print("nub0 :", result.nub0)
        print("nub  :", result.nub)
        print("VQas :", result.VQas)
        print("DQas :", result.DQas)
        print("VQas2:", result.VQas2)
        print("DQas2:", result.DQas2)
        print("Gamma_gyro:", result.Gamma_gyro)

    print(f"\nProcessed {len(all_results)} simulations.")


if __name__ == "__main__":
    main()

