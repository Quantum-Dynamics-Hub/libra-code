"""
NBRA Phase-1 example for LiF CASSCF(2e,5o), 2 singlet roots.

Phase 1:
    1. Define one configured ES strategy prototype.
    2. Define a predetermined Li-F bond-stretch trajectory.
    3. Call strategy_compute_adi() for every geometry.
    4. Collect:
         - adiabatic energies
         - neighboring-frame time overlaps
         - Hvib matrices
    5. Write a human-readable CSV table.
    6. Write an Hvib-only NPZ file for Phase 2.

strategy_compute_adi() internally owns the lifetime of the
previous/current ES strategy snapshots.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

from liblibra_core import MATRIX

from libra_py import units


# ============================================================================
# LOCAL REPOSITORY IMPORT
# ============================================================================

def _prepend_repo_root() -> None:
    file_path = Path(__file__).resolve()

    for parent in file_path.parents:
        if (parent / "src" / "libra_py" / "__init__.py").is_file():
            repo_root = str(parent / "src")

            if repo_root not in sys.path:
                sys.path.insert(0, repo_root)

            build_root = str(parent / "_build_venv" / "src")
            if (parent / "_build_venv" / "src").exists() and build_root not in sys.path:
                sys.path.insert(0, build_root)

            return


if __name__ == "__main__" and __package__ is None:
    _prepend_repo_root()


# ============================================================================
# IMPORTS
# ============================================================================

from libra_py.packages.pyscf.implementations.casscf import CASSCF
from libra_py.packages.pyscf.methods import strategy_compute_adi


# ============================================================================
# CHEMICAL SYSTEM
# ============================================================================

NSTATES = 2

ATOM_LABELS = [
    "Li",
    "F",
]

BASIS = {
    "Li": "sto-3g",
    "F": "6-311+g*",
}

CAS_LIST = [
    4,
    7,
    11,
    14,
    17,
]

NORBCAS = 5
NELECAS = 2


# ============================================================================
# PREDETERMINED NUCLEAR TRAJECTORY
# ============================================================================

R_MIN = 7.0       # Bohr
R_MAX = 11.0      # Bohr

# Spatial separation between neighboring ES calculations.
DR = 0.1        # Bohr

# Physical time separating neighboring frames.
DT_FS = 0.2       # fs

AU_PER_FS = 41.3413745758
DT_AU = DT_FS * AU_PER_FS


bond_grid = np.arange(
    R_MIN,
    R_MAX + 1e-12,
    DR,
)

# ============================================================================
# CONFIGURED ES STRATEGY PROTOTYPE
# ============================================================================
#
# The input file chooses the ES method and configures it.
#
# This object is a prototype. strategy_compute_adi() does not use this
# object itself as a trajectory snapshot. Instead, it creates independent
# snapshots using:
#
#     current = es_strategy.clone()
#
# This allows strategy_compute_adi() to maintain:
#
#     previous = ES(R[n-1])
#     current  = ES(R[n])
#
# for computing neighboring-frame time overlaps.
# ============================================================================

es_obj = CASSCF(
    norbcas=NORBCAS,
    nelecas=NELECAS,
    nroots=NSTATES,
    basis=BASIS,
    unit="Bohr",
    charge=0,
    cas_list=CAS_LIST,
)

# ============================================================================
# MODEL PARAMETERS
# ============================================================================
#
# User-owned configuration:
#
#     es_strategy
#     atom_labels
#     nstates
#     requested ES properties
#
# Runtime state such as "_es_previous" is created and managed internally
# by strategy_compute_adi().
# ============================================================================

model_params = {

    "atom_labels": ATOM_LABELS,

    "nstates": NSTATES,

    # Configured ES strategy prototype
    "es_strategy": es_obj,

    # No gradients needed for NBRA Phase 1
    "gradient_state": None,

    "hessian_state": None,

    # Use neighboring-frame time overlaps
    "nacv": False,

    "time_overlap": True,

    "H_soc": False,

    # Used for time-overlap -> Hvib conversion
    "dt": DT_AU,
}


# ============================================================================
# BUILD ONE LIBRA q FRAME
# ============================================================================

def make_q(r_lif_bohr: float) -> MATRIX:

    coords_bohr = [
        0.0, 0.0, 0.0,
        0.0, 0.0, r_lif_bohr,
    ]

    q = MATRIX(
        len(coords_bohr),
        1,
    )

    for i, value_bohr in enumerate(coords_bohr):

        q.set(
            i,
            0,
            value_bohr * units.Angst,
        )

    return q


# ============================================================================
# LIBRA CMATRIX -> NUMPY
# ============================================================================

def cmatrix_to_numpy(
    matrix,
    nstates: int,
) -> np.ndarray:

    return np.array(
        [
            [
                matrix.get(i, j)
                for j in range(nstates)
            ]
            for i in range(nstates)
        ],
        dtype=np.complex128,
    )


# ============================================================================
# STORAGE
# ============================================================================

frames_fs = []
energies = []
overlaps = []
hvibs = []


# ============================================================================
# NBRA PHASE-1 LOOP
# ============================================================================
#
# The high-level logic is intentionally simple:
#
#     geometry
#         ↓
#     strategy_compute_adi()
#         ↓
#     energy + time overlap + Hvib
#         ↓
#     save
#
# The driver does not know about:
#
#     ES_Request
#     ES_Result
#     clone()
#     previous strategy
#     current strategy
#
# ============================================================================


def complex_parts(z):
 
     if z is None:
         return [
             "",
             "",
         ]
 
     return [
         float(np.real(z)),
         float(np.imag(z)),
     ]

def overlap_parts(z):
    
    if z is None:
        return [
            "",
        ]
    
    return [
        float(np.real(z)),
    ]

for istep, r_lif in enumerate(bond_grid):

    time_fs = istep * DT_FS

    frames_fs.append(
        time_fs
    )


    # ------------------------------------------------------------------------
    # 1. Build the geometry frame
    # ------------------------------------------------------------------------

    q = make_q(
        float(r_lif)
    )


    # ------------------------------------------------------------------------
    # 2. Run ES calculation through the generic adapter
    # ------------------------------------------------------------------------

    result = strategy_compute_adi(
        q,
        model_params,
        [0],
    )


    # ------------------------------------------------------------------------
    # 3. Extract energy of every adiabatic root
    # ------------------------------------------------------------------------

    frame_energies = np.array(
        [
            result.ham_adi.get(
                i,
                i,
            ).real

            for i in range(NSTATES)
        ],
        dtype=np.float64,
    )


    energies.append(
        frame_energies
    )


    # ------------------------------------------------------------------------
    # 4. First frame has no previous-frame overlap
    # ------------------------------------------------------------------------

    if istep == 0:

        overlaps.append(
            None
        )

        hvibs.append(
            None
        )

        continue


    # ------------------------------------------------------------------------
    # 5. Extract S(R[n-1], R[n])
    # ------------------------------------------------------------------------

    frame_overlap = cmatrix_to_numpy(
        result.time_overlap_adi,
        NSTATES,
    )


    overlaps.append(
        frame_overlap
    )


    # ------------------------------------------------------------------------
    # 6. Extract Hvib
    #
    # Hvib was already constructed by strategy_compute_adi().
    # Do not reconstruct it here.
    # ------------------------------------------------------------------------

    frame_hvib = cmatrix_to_numpy(
        result.hvib_adi,
        NSTATES,
    )


    hvibs.append(
        frame_hvib
    )


# ============================================================================
# WRITE HUMAN-READABLE PHASE-1 TABLE
# ============================================================================

output_dir = Path(__file__).resolve().parents[3]
output_dir.mkdir(parents=True, exist_ok=True)

csv_file = output_dir / "lif_nbra_phase1_table.csv"


with open(
    csv_file,
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.writer(f)


    writer.writerow([

        "frame",
        "time_fs",
        "R_bohr",

        "E0_au",
        "E1_au",

        "S00_re",
        "S00_im",

        "S01_re",
        "S01_im",

        "S10_re",
        "S10_im",

        "S11_re",
        "S11_im",

        "Hvib00_re",
        "Hvib00_im",

        "Hvib01_re",
        "Hvib01_im",

        "Hvib10_re",
        "Hvib10_im",

        "Hvib11_re",
        "Hvib11_im",
    ])


    for iframe, r_lif in enumerate(bond_grid):

        e = energies[iframe]
        s = overlaps[iframe]
        h = hvibs[iframe]


        row = [

            iframe,

            float(
                frames_fs[iframe]
            ),

            float(
                r_lif
            ),

            float(
                e[0]
            ),

            float(
                e[1]
            ),
        ]


        # Time overlap

        if s is None:

            row += [
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
            ]

        else:

            row += (
                overlap_parts(s[0, 0])
                + overlap_parts(s[0, 1])
                + overlap_parts(s[1, 0])
                + overlap_parts(s[1, 1])
            )


        # Hvib

        if h is None:

            row += [
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
            ]

        else:

            row += (
                complex_parts(h[0, 0])
                + complex_parts(h[0, 1])
                + complex_parts(h[1, 0])
                + complex_parts(h[1, 1])
            )


        writer.writerow(
            row
        )


# ============================================================================
# WRITE Hvib-ONLY DATASET FOR NBRA PHASE 2
# ============================================================================

Hvib_table = np.array(
    hvibs[1:],
    dtype=np.complex128,
)


np.savez(

    output_dir / "lif_hvib_only.npz",

    Hvib_table=Hvib_table,

    dt_fs=np.array(
        DT_FS
    ),

    dt_au=np.array(
        DT_AU
    ),

    bond_grid=bond_grid,

    time_fs=np.array(
        frames_fs
    ),
)


# ============================================================================
# SUMMARY
# ============================================================================

print(
    "Wrote lif_nbra_phase1_table.csv"
)

print(
    "Wrote lif_hvib_only.npz"
)

print(
    f"Number of geometry frames = {len(bond_grid)}"
)

print(
    f"Number of neighboring-frame overlaps = {len(overlaps) - 1}"
)

print(
    f"Number of Hvib matrices = {len(Hvib_table)}"
)

print(
    f"Coordinate spacing DR = {DR} Bohr"
)

print(
    f"Frame time spacing DT_FS = {DT_FS} fs"
)

print(
    f"Effective prescribed velocity = "
    f"{DR / DT_FS:.6f} Bohr/fs"
)