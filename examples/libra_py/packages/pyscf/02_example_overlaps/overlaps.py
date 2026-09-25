"""Compare the available PySCF overlap strategies on two small geometries.

Run with an installed Libra/PySCF environment. --source-tree optionally selects
this checkout's libra_py sources while retaining its installed compiled modules.
"""
import argparse
import json
from pathlib import Path

import numpy as np

METHODS = ("casscf_hf", "casscf", "cisd", "tda", "tddft")
SYSTEMS = {"lih": dict(charge=0, multiplicity=1, nelecas=(1, 1)),
           "lih_plus": dict(charge=1, multiplicity=2, nelecas=(1, 0))}


def make_strategy(method, system):
    from libra_py.packages.pyscf.implementations import CASSCF, CISD, TDDFT
    common = dict(basis="sto-3g", unit="Bohr", charge=system["charge"])
    if method.startswith("casscf"):
        return CASSCF(norbcas=2, nelecas=system["nelecas"], nroots=2,
                      spin_multiplicity=system["multiplicity"],
                      overlap_orbitals="hf" if method == "casscf_hf" else "casscf", **common)
    if method == "cisd":
        return CISD(nroots=2, spin_multiplicity=system["multiplicity"], **common)
    return TDDFT(atom_labels=("Li", "H"), nexc=1, spin=system["multiplicity"]-1,
                 xc="pbe0", grid_level=1, use_tda=method == "tda", **common)


def check_convergence(state):
    for name in ("mf", "mc", "myci", "td"):
        solver = getattr(state, name, None)
        if solver is not None and not np.all(np.atleast_1d(getattr(solver, "converged", True))):
            raise RuntimeError(f"Unconverged {name}; overlap comparison stopped")


def orbital_metric_error(state):
    # Report full-MO orthonormality, not just CI coefficient normalization.
    mo = state.mc.mo_coeff if getattr(state, "mc", None) is not None else state.mf.mo_coeff
    blocks = [mo] if np.asarray(mo).ndim == 2 else mo
    return max(float(np.linalg.norm(c.T.conj() @ state.mf.get_ovlp() @ c - np.eye(c.shape[1])))
               for c in blocks)


def calculate(method, system, displacement, sequential):
    from libra_py.packages.pyscf.interfaces import ES_Request, MolecularGeometry
    previous_geom = MolecularGeometry(("Li", "H"), np.array([[0., 0., 0.], [0., 0., 3.0]]))
    current_geom = MolecularGeometry(("Li", "H"), np.array([[0., 0., 0.], [0., 0., 3.0+displacement]]))
    request = ES_Request(n_singlets=2, gradient_state=None, nacv=False, time_overlap=True)
    previous_strategy = make_strategy(method, system)
    previous_result = previous_strategy.compute_result(previous_geom, request)
    check_convergence(previous_strategy.get_state())
    # Preserve a separate snapshot before any sequential update.
    previous = previous_strategy.copy().get_state()
    current_strategy = previous_strategy if sequential else make_strategy(method, system)
    current_result = current_strategy.compute_result(current_geom, request)
    current = current_strategy.get_state()
    check_convergence(current)
    overlap = current_strategy.compute_time_overlap(current, previous)
    # In a sequential run this also checks the public request/result route.
    if sequential:
        np.testing.assert_allclose(overlap, current_result.time_overlap, atol=1e-10)
    self_overlap = current_strategy.compute_time_overlap(current, current)
    if not np.all(np.isfinite(overlap)):
        raise ValueError("Nonfinite overlap")
    return dict(method=method, charge=system["charge"], multiplicity=system["multiplicity"],
                basis="sto-3g", distances_bohr=[3.0, 3.0+displacement],
                initialization="sequential" if sequential else "independent",
                energies_previous_hartree=np.asarray(previous_result.H_el).tolist(),
                energies_current_hartree=np.asarray(current_result.H_el).tolist(),
                overlap=overlap.tolist(), self_overlap_current=self_overlap.tolist(),
                overlap_singular_values=np.linalg.svd(overlap, compute_uv=False).tolist(),
                overlap_orthogonality_defect=float(np.linalg.norm(overlap.T @ overlap-np.eye(2))),
                self_overlap_identity_error=float(np.linalg.norm(self_overlap-np.eye(2))),
                previous_mo_metric_error=orbital_metric_error(previous),
                current_mo_metric_error=orbital_metric_error(current))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--system", choices=("all", *SYSTEMS), default="all")
    p.add_argument("--method", choices=("all", *METHODS), default="all")
    p.add_argument("--displacement", type=float, default=0.1, help="Li-H displacement in Bohr")
    p.add_argument("--sequential", action="store_true", help="Reuse the strategy between geometries, as in dynamics")
    p.add_argument("--source-tree", action="store_true", help="Use this checkout's Python sources")
    p.add_argument("--output", type=Path, default=Path("overlap_results"), help="Fresh result directory")
    args = p.parse_args()
    if not np.isfinite(args.displacement) or 3.0+args.displacement <= 0:
        p.error("The displaced bond length must be finite and positive")
    if args.source_tree:
        import libra_py
        repo = Path(__file__).resolve().parents[5]
        libra_py.__path__.insert(0, str(repo / "src/libra_py"))
    from pyscf import lib
    import pyscf
    lib.num_threads(1)
    args.output.mkdir(parents=True, exist_ok=False)
    systems = SYSTEMS if args.system == "all" else {args.system: SYSTEMS[args.system]}
    methods = METHODS if args.method == "all" else (args.method,)
    results = {}
    for name, system in systems.items():
        for method in methods:
            key = f"{name}_{method}"
            row = calculate(method, system, args.displacement, args.sequential)
            row["pyscf_version"] = pyscf.__version__
            results[key] = row
            (args.output / f"{key}.json").write_text(json.dumps(row, indent=2)+"\n")
            np.savetxt(args.output / f"{key}_overlap.csv", row["overlap"], delimiter=",")
            print(f"\n{key}: <previous i | current j> (column signs aligned)")
            print(np.asarray(row["overlap"]))
            print(f"MO metric error: {row['current_mo_metric_error']:.3e}; "
                  f"self-overlap error: {row['self_overlap_identity_error']:.3e}")
    (args.output / "summary.json").write_text(json.dumps(results, indent=2)+"\n")


if __name__ == "__main__":
    main()
