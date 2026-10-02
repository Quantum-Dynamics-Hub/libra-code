"""Shared numerical work for the manuscript LDR/DVR examples.

Entrypoints and settings live in each model directory. Paths are based on this
file, so the examples can run from a source checkout without PYTHONPATH.
"""

import argparse
import copy
import json
from pathlib import Path
import sys

import torch

SRC = Path(__file__).resolve().parents[4] / "src"
sys.path.insert(0, str(SRC))
from libra_py.dynamics.ldr_torch.compute import ldr_solver
from libra_py.dynamics.exact_torch.compute import exact_tdse_solver_multistate


def potential(q, settings):
    """Float64 SAC, DAC, and FLV potentials; q has shape (ndof, *grid)."""
    x = q[0]
    h = torch.zeros((*x.shape, 2, 2), dtype=torch.float64, device=x.device)
    model = settings["model"]
    if model == "tully1":
        h[..., 0, 0] = 0.01 * x.sign() * (1 - torch.exp(-1.6 * x.abs()))
        h[..., 1, 1] = -h[..., 0, 0]
        coupling = 0.005 * torch.exp(-x**2)
    elif model == "tully2":
        h[..., 1, 1] = 0.05 - 0.10 * torch.exp(-0.28 * x**2)
        coupling = 0.015 * torch.exp(-0.06 * x**2)
    elif model == "flv":
        y = q[1]
        h[..., 0, 0] = 0.01 * (x - 4)**2 + 0.05 * y**2
        h[..., 1, 1] = 0.01 * (x - 3)**2 + 0.05 * y**2 + 0.01
        coupling = settings["gamma"] * y * torch.exp(-3 * (x - 3)**2 - 1.5 * y**2)
    else:
        raise ValueError(f"Unknown model: {model}")
    h[..., 0, 1] = h[..., 1, 0] = coupling
    return h.to(torch.cdouble)


def grid_axes(settings, method, device):
    axes = []
    for (lo, hi), dx in zip(settings[f"{method}_bounds"], settings[f"{method}_spacing"]):
        intervals = round((hi - lo) / dx)
        if intervals < 1 or abs(intervals * dx - (hi - lo)) > 1e-8:
            raise ValueError("Grid spacing must divide each specified interval")
        axes.append(torch.linspace(lo, hi, intervals + 1, dtype=torch.float64, device=device))
    return axes


def coordinates(axes):
    return torch.stack(torch.meshgrid(*axes, indexing="ij"))


def align_1d(frames):
    """Align adjacent real-model eigenvectors, retaining eigenvectors as columns."""
    overlaps = (frames[:-1].conj() * frames[1:]).sum(dim=-2)
    factors = torch.where(overlaps.abs() > 1e-14,
                          overlaps.conj() / overlaps.abs().clamp_min(1e-14),
                          torch.ones_like(overlaps))
    phases = torch.cat([torch.ones_like(factors[:1]), factors.cumprod(dim=0)])
    return frames * phases[:, None, :]


def electronic_data(q, settings):
    energies, frames = torch.linalg.eigh(potential(q, settings))
    if q.shape[0] == 1:
        frames = align_1d(frames)
    return energies, frames


def reference_state(settings, device):
    q0 = torch.tensor(settings["q0"], dtype=torch.float64, device=device)[:, None]
    _, frames = torch.linalg.eigh(potential(q0, settings))
    return frames[0, :, settings["istate"]]


def gaussian(q, settings):
    """Normalized nuclear envelope with the specified probability std. deviations."""
    shape = (len(settings["q0"]),) + (1,) * (q.ndim - 1)
    q0 = q.new_tensor(settings["q0"]).reshape(shape)
    p0 = q.new_tensor(settings["p0"]).reshape(shape)
    sigma = q.new_tensor(settings["sigma"]).reshape(shape)
    alpha0 = 1 / (4 * sigma**2)
    prefactor = torch.prod((2 * alpha0 / torch.pi)**0.25)
    return prefactor * torch.exp((-alpha0 * (q - q0)**2 + 1j * p0 * (q - q0)).sum(dim=0))


class DVRReference(exact_tdse_solver_multistate):
    """Stream a double-precision SOFT reference using exact_torch's operators.

    The example supplies its own grids and recording loop: use fftfreq for
    correct odd/even FFT ordering, include both t=0 and the final time, and
    store selected densities instead of every representation at every time.
    The core electronic transforms and potential/kinetic operators are reused.
    """

    def __init__(self, settings, device):
        axes = grid_axes(settings, "dvr", device)
        super().__init__({
            "grid_size": [len(a) for a in axes], "mass": settings["mass"],
            "dt": settings["dt"], "Nstates": 2, "device": device,
            "potential_fn": potential, "potential_fn_params": settings,
            "psi0_fn": gaussian, "psi0_fn_params": settings,
            "representation": "diabatic", "method": "split-operator",
        })
        self.settings = settings
        self.axes = axes
        self.initialize_grids()
        self.initialize_operators()
        if self.ndim == 1:
            self.eigvecs = align_1d(self.eigvecs)
        # Identical physical preparation to LDR: chi0(q) phi_istate(q0).
        self.psi_r_dia = gaussian(self.Q, settings)[..., None] * reference_state(settings, device)
        self.psi_r_dia /= (self.psi_r_dia.abs().square().sum() * self.dV).sqrt()
        self.update_adi_r()
        self.transform_r2k(0)

    def initialize_grids(self):
        self.grid_size = torch.tensor([len(a) for a in self.axes])
        self.ngrid = int(self.grid_size.prod())
        self.Q = coordinates(self.axes)
        self.dq = torch.stack([a[1] - a[0] for a in self.axes])
        self.dV = self.dq.prod()
        self.mass = self.Q.new_tensor(self.settings["mass"])
        k_axes = [2 * torch.pi * torch.fft.fftfreq(len(a), d=float(a[1] - a[0]),
                                                 dtype=torch.float64, device=a.device)
                  for a in self.axes]
        self.K = coordinates(k_axes)
        shape = tuple(len(a) for a in self.axes) + (2,)
        self.psi_r_dia = torch.zeros(shape, dtype=torch.cdouble, device=self.Q.device)
        self.psi_r_adi = torch.zeros_like(self.psi_r_dia)
        self.psi_k_dia = torch.zeros_like(self.psi_r_dia)
        self.psi_k_adi = torch.zeros_like(self.psi_r_dia)

    def step(self):
        self.psi_r_dia = torch.einsum("...ij,...j->...i", self.expV_half, self.psi_r_dia)
        self.transform_r2k(0)
        self.psi_k_dia *= self.expT[..., None]
        self.transform_k2r(0)
        self.psi_r_dia = torch.einsum("...ij,...j->...i", self.expV_half, self.psi_r_dia)


def make_ldr(settings, device):
    axes = grid_axes(settings, "ldr", device)
    q = coordinates(axes).flatten(start_dim=1)
    energies, frames = electronic_data(q, settings)
    n = q.shape[1]
    se = torch.einsum("ndi,mdj->injm", frames.conj(), frames).reshape(2*n, 2*n)
    b_el = torch.einsum("ndi,d->in", frames.conj(), reference_state(settings, device))
    sigma = q.new_tensor(settings["sigma"])
    mass = q.new_tensor(settings["mass"])
    alpha = q.new_tensor([1 / (2 * float(a[1] - a[0])**2) for a in axes])
    solver = ldr_solver({
        "device": device, "qgrid": q.T, "alpha": alpha,
        "q0": settings["q0"], "p0": settings["p0"], "mass": mass,
        "k": 1 / (4 * mass * sigma**4), "nstates": 2, "istate": settings["istate"],
        "E": energies.T, "s_elec": se, "elec_ampl": b_el, "dt": settings["dt"],
    })
    solver.buildSH()
    solver.compute_propagator()
    solver.initialize_C()
    solver.C_curr = solver.C0.clone()
    return solver, axes, frames


def ldr_density(solver, axes, frames, evaluation_axes):
    """Reconstruct total nuclear density, including electronic interference."""
    c_dia = torch.einsum("ndi,in->nd", frames, solver.C_curr.reshape(2, -1))
    kernels = [(2 * solver.alpha[f] / torch.pi)**0.25 *
               torch.exp(-solver.alpha[f] * (a[:, None] - axes[f][None, :])**2)
               for f, a in enumerate(evaluation_axes)]
    kernels = [k.to(torch.cdouble) for k in kernels]
    if len(axes) == 1:
        psi = kernels[0] @ c_dia
    else:
        c_dia = c_dia.reshape(len(axes[0]), len(axes[1]), 2)
        psi = torch.einsum("ax,xyi,by->abi", kernels[0], c_dia, kernels[1])
    norm = torch.vdot(solver.C_curr, solver.S @ solver.C_curr).real
    return psi.abs().square().sum(dim=-1) / norm


def record_ldr(solver):
    c = solver.C_curr
    sc = solver.S @ c
    norm = torch.vdot(c, sc).real
    z = (solver.S_half @ c).reshape(2, -1)
    rho = z @ z.conj().T / norm
    total = torch.vdot(c, solver.H @ c).real / norm
    # Exact expectation of the implemented symmetrized endpoint potential.
    pe = torch.vdot(solver.E.reshape(-1) * c, sc).real / norm
    pos = torch.stack([torch.vdot(solver.qgrid[:, f].repeat(2) * c, sc).real / norm
                       for f in range(solver.ndof)])
    return rho, total - pe, pe, total, pos, norm


def record_dvr(solver):
    solver.update_adi_r()
    solver.transform_r2k(0)
    psi = solver.psi_r_dia
    density = psi.abs().square().sum(dim=-1)
    norm = density.sum() * solver.dV
    a = solver.psi_r_adi.reshape(-1, 2)
    rho = a.T @ a.conj() * solver.dV / norm
    ke = (solver.psi_k_dia.abs().square() * solver.T[..., None]).sum()
    ke = ke * solver.dV / solver.ngrid / norm
    pe = torch.einsum("...i,...ij,...j->", psi.conj(), solver.V, psi).real * solver.dV / norm
    pos = torch.stack([(solver.Q[f] * density).sum() * solver.dV / norm
                       for f in range(solver.ndim)])
    return rho, ke, pe, ke + pe, pos, norm


def simulate(settings, method, output_dir, device):
    dt = settings["dt"]
    steps = round(settings["duration"] / dt)
    snapshot_steps = {round(t / dt) for t in settings["snapshots"]}
    if abs(steps * dt - settings["duration"]) > 1e-8 or any(
            abs(round(t / dt) * dt - t) > 1e-8 for t in settings["snapshots"]):
        raise ValueError("Duration and snapshot times must be multiples of dt")
    if any(s < 0 or s > steps for s in snapshot_steps):
        raise ValueError("Snapshot times must fall within the simulation")
    evaluation_axes = grid_axes(settings, "dvr", device)
    print(f"{settings['title']}: {method.upper()}, {steps} steps, dt={dt}", flush=True)
    if method == "ldr":
        solver, axes, frames = make_ldr(settings, device)
    else:
        solver = DVRReference(settings, device)
    names = ["rho", "kinetic_energy", "potential_energy", "total_energy", "position", "norm"]
    series = {name: [] for name in names}
    times, densities, density_times = [], [], []
    for step in range(steps + 1):
        if step % settings["save_every"] == 0 or step == steps:
            values = record_ldr(solver) if method == "ldr" else record_dvr(solver)
            for name, value in zip(names, values):
                series[name].append(value.detach().cpu())
            times.append(step * dt)
        if step in snapshot_steps:
            if method == "ldr":
                density = ldr_density(solver, axes, frames, evaluation_axes)
            else:
                density = solver.psi_r_dia.abs().square().sum(dim=-1)
                density = density / (density.sum() * solver.dV)
            densities.append(density.detach().cpu())
            density_times.append(step * dt)
        if step % max(1, steps // 10) == 0:
            print(f"  {method}: t={step * dt:g}", flush=True)
        if step < steps:
            if method == "ldr":
                solver.C_curr = solver.U @ solver.C_curr
            else:
                solver.step()
    result = {name: torch.stack(value) for name, value in series.items()}
    result.update(time=torch.tensor(times, dtype=torch.float64),
                  density=torch.stack(densities),
                  density_time=torch.tensor(density_times, dtype=torch.float64),
                  axes=[a.cpu() for a in evaluation_axes], settings=settings, method=method)
    output_dir.mkdir(parents=True, exist_ok=True)
    torch.save(result, output_dir / f"{method}.pt")
    (output_dir / f"{method}_settings.json").write_text(json.dumps(settings, indent=2) + "\n")
    print(f"Saved {output_dir / (method + '.pt')}", flush=True)
    return result


def quick_settings(settings):
    settings = copy.deepcopy(settings)
    settings["title"] += " (quick execution check)"
    settings.update(duration=20.0, save_every=5, snapshots=[0.0, 10.0, 20.0])
    if settings["model"] == "flv":
        settings["ldr_spacing"] = [0.3, 0.3]
        settings["dvr_spacing"] = [0.1, 0.1]
    else:
        # Retain momentum resolution for p0=30; a coarser FFT grid aliases
        # this wavepacket. Reduce the box instead for the short smoke run.
        settings["ldr_bounds"] = settings["dvr_bounds"] = [[-12.0, 4.0]]
        settings["ldr_spacing"] = settings["dvr_spacing"] = [0.05]
    return settings


def run_cases(cases, output_dir):
    parser = argparse.ArgumentParser(description="Run the local manuscript LDR/DVR examples")
    parser.add_argument("--quick", action="store_true", help="Small execution check in output/quick")
    parser.add_argument("--method", choices=["both", "ldr", "dvr"], default="both")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("threads must be positive")
    torch.set_num_threads(args.threads)
    output_dir = output_dir / "quick" if args.quick else output_dir
    for name, settings in cases.items():
        config = quick_settings(settings) if args.quick else copy.deepcopy(settings)
        folder = output_dir / name if name else output_dir
        methods = ["ldr", "dvr"] if args.method == "both" else [args.method]
        for method in methods:
            simulate(config, method, folder, torch.device(args.device))
