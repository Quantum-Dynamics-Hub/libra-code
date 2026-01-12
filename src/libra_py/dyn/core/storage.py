"""
 rank-generalizable
 backend-agnostic

"""

class TensorStorage:
    """
    Shared numerical storage for all TBFs.
    Lives on CPU/GPU (NumPy, PyTorch, JAX).
    """

    def __init__(self, ndof, nstates, device="cpu"):
        self.ndof = ndof
        self.nstates = nstates

        # Core nuclear DOFs
        self.R = None        # (N_tbf, ndof)
        self.P = None        # (N_tbf, ndof)

        # Electronic structure
        self.C = None        # (N_tbf, nstates)
        self.H_adi = None    # (N_tbf, nstates, nstates)
        self.H_dia = None    # (N_tbf, nstates, nstates)

        # Forces, energies, etc.
        self.forces = None   # (N_tbf, ndof)
        self.energies = None # (N_tbf, nstates)

