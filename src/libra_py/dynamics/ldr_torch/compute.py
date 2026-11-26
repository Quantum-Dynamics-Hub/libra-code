# *********************************************************************************
# * Copyright (C) 2025 Daeho Han and Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# ***********************************************************************************
"""
.. module:: compute
   :platform: Unix, Windows
   :synopsis: This module implements functions for doing local diabatic representation (LDR) dynamics with PyTorch
       List of classes:
           * ldr_tdse_solver_multistate

.. moduleauthor:: Daeho Han and Alexey V. Akimov

"""

__author__ = "Daeho Han, Alexey V. Akimov"
__copyright__ = "Copyright 2025 Alexey V. Akimov"
__credits__ = ["Daeho Han", "Alexey V. Akimov"]
__license__ = "GNU-3"
__version__ = "1.0"
__maintainer__ = "Alexey V. Akimov"
__email__ = "alexvakimov@gmail.com"
__url__ = "https://github.com/Quantum-Dynamics-Hub/libra-code"



import torch

from . import electronic


class ldr_tdse_solver_multistate:
    """LDR dynamics driven by a diabatic model function, as in exact_torch.

    Supply ``potential_fn(q, potential_fn_params)`` returning shape
    (npoints, nstates, nstates) for q of shape (ndof, npoints). Libra prepares
    reference energies, eigenvectors, electronic overlaps, and initial
    electronic overlaps in a fixed orthonormal diabatic model basis.

    ``qgrid`` has shape (ngrids, ndof). Alternatively, ``grid_size``, ``q_min``
    and ``q_max`` define a regular grid with both endpoints included.
    ``alpha`` is a positive scalar, dimension vector, or (ngrids, ndof) array.
    ``Nstates`` and ``initial_state_index`` follow the exact_torch names.

    The initial nuclear state is the analytic Gaussian specified by q0, p0,
    k and mass. With potential_fn, ``initial_representation="diabatic"``
    (the default) selects a unit vector in the model's diabatic basis;
    ``"adiabatic"`` selects the adiabatic eigenvector at q0. Both use
    ``initial_state_index`` and hold the electronic vector fixed across the
    Gaussian. This is not a locally adiabatic state at every coordinate.
    ``initial_electronic_state`` overrides this selection with a normalized
    vector supplied in the model's diabatic basis, regardless of
    ``initial_representation``.

    Schemes: diagonal, BAT0 (alias symmetrized, the default), BAT1, midpoint.
    BAT1 computes energy gradients at the reference geometries from
    ``potential_fn`` using Torch autograd in a separate helper. The model
    must preserve the computational graph from coordinates to potential.
    Midpoint always evaluates potential_fn at Gaussian-weighted midpoints
    through a separate helper; ``midpoint_batch_size`` defaults to 256.

    Existing calls with precomputed E, s_elec and elec_ampl remain supported
    when potential_fn is omitted. These arrays use state-major order.

    ``psi_r_dia`` and ``psi_r_adi`` are ordinary (ngrids, nstates) arrays
    evaluated on qgrid. Enable their (T, ngrids, nstates) histories with
    "psi_r_dia" and "psi_r_adi" in properties_to_save. ``rho_dia_all`` and
    ``rho_adi_all`` are initialized to None as placeholders for a future
    electronic density-matrix calculation. Spatial densities are abs(psi)**2.
    Without potential_fn, only native-center adiabatic amplitudes exist.
    """

    def __init__(self, params):
        params = dict(params)
        if "nstates" in params or "istate" in params:
            raise ValueError("Use Nstates and initial_state_index as parameter names")
        self.prefix = params.get("prefix", "ldr-solution")
        self.device = torch.device(params.get("device", "cuda" if torch.cuda.is_available() else "cpu"))
        self.hbar = 1.0
        self.hamiltonian_scheme = params.get("hamiltonian_scheme", "symmetrized")
        self.dE = params.get("dE")
        self.potential_fn = params.get("potential_fn")
        self.potential_fn_params = params.get("potential_fn_params", {})
        self.initial_representation = params.get("initial_representation", "diabatic")
        if self.potential_fn is None and "initial_representation" in params:
            raise ValueError("initial_representation requires potential_fn; "
                             "use elec_ampl with precomputed electronic data")
        if "potential_gradient_fn" in params:
            raise ValueError("BAT1 uses autograd with potential_fn; remove potential_gradient_fn")
        self.eigvecs = None
        self.midpoint_batch_size = params.get("midpoint_batch_size", 256)
        self.q0 = torch.as_tensor(params.get("q0", [0.0]), dtype=torch.float64, device=self.device)
        self.p0 = torch.as_tensor(params.get("p0", [0.0]), dtype=torch.float64, device=self.device)
        self.k = torch.as_tensor(params.get("k", [0.001]), dtype=torch.float64, device=self.device)
        self.mass = torch.as_tensor(params.get("mass", [2000.0]), dtype=torch.float64, device=self.device)
        self.alpha = torch.as_tensor(params.get("alpha", [18.0]), dtype=torch.float64, device=self.device)
        if "grid_size" in params:
            if "qgrid" in params:
                raise ValueError("Supply qgrid or grid_size, not both")
            sizes, lower, upper = params["grid_size"], params.get("q_min"), params.get("q_max")
            if lower is None or upper is None or not (len(sizes) == len(lower) == len(upper)):
                raise ValueError("grid_size, q_min and q_max must have the same length")
            if not sizes or any(not isinstance(n, int) or n < 1 for n in sizes):
                raise ValueError("grid_size must contain positive integers")
            if any(hi <= lo for lo, hi in zip(lower, upper)):
                raise ValueError("q_max must exceed q_min in every dimension")
            axes = [torch.linspace(lo, hi, n, dtype=torch.float64, device=self.device)
                    for lo, hi, n in zip(lower, upper, sizes)]
            params["qgrid"] = torch.stack(torch.meshgrid(*axes, indexing="ij"), dim=-1).reshape(-1, len(sizes))
        self.qgrid = torch.as_tensor(
            params.get("qgrid", [[-10 + i*0.1] for i in range(201)]),
            dtype=torch.float64, device=self.device)
        if self.qgrid.ndim != 2 or min(self.qgrid.shape) < 1 or not torch.isfinite(self.qgrid).all():
            raise ValueError("qgrid must be a finite, nonempty (ngrids, ndof) array")
        self.ngrids = len(self.qgrid) # N
        self.ndof = self.qgrid.shape[1] 
        if not torch.all(torch.isfinite(self.alpha) & (self.alpha > 0)):
            raise ValueError("Gaussian widths must be positive and finite")
        if self.alpha.ndim == 0 or (
                self.alpha.ndim == 1 and self.alpha.numel() in (1, self.ndof)):
            self._alpha_grid = self.alpha.expand(self.ngrids, self.ndof)
        elif self.alpha.shape == (self.ngrids, self.ndof):
            self._alpha_grid = self.alpha
        else:
            raise ValueError("alpha must be scalar, have shape (ndof,), "
                             "or have shape (ngrids, ndof)")
        self.nstates = params.get("Nstates", 2)
        self.istate = params.get("initial_state_index", 0)
        if not 0 <= self.istate < self.nstates:
            raise ValueError("initial_state_index must index one of Nstates electronic states")
        self.Nstates, self.initial_state_index = self.nstates, self.istate
        if self.potential_fn is not None:
            array_keys = {"E", "s_elec", "elec_ampl", "dE", "midpoint_potential", "midpoint_evaluator"}
            if array_keys.intersection(params):
                raise ValueError("potential_fn generates electronic data; remove precomputed electronic inputs")
            with torch.no_grad():
                E, self.eigvecs, se, ampl = electronic.prepare_reference_data(
                    self.qgrid, self.q0, self.nstates, self.istate,
                    self.potential_fn, self.potential_fn_params, params.get("initial_electronic_state"),
                    initial_representation=self.initial_representation)
            params.update(E=E, s_elec=se, elec_ampl=ampl)
            if self.hamiltonian_scheme == "BAT1":
                # Continue with the eigensolver-selected reference basis at degeneracies.
                # The resulting Hellmann-Feynman gradients depend on that basis.
                # if torch.any((E[1:] - E[:-1]).abs() < 1e-12):
                #     raise ValueError("BAT1 energy gradients are not unique at an electronic degeneracy")
                self.dE = electronic.compute_energy_gradients(
                    self.qgrid, self.eigvecs, self.potential_fn, self.potential_fn_params)
        amplitudes = torch.as_tensor(
            params.get("elec_ampl", torch.ones(self.ngrids)),
            dtype=torch.cdouble, device=self.device,
        )
        if amplitudes.shape == (self.ngrids,):
            self.elec_ampl = torch.zeros(
                self.nstates, self.ngrids, dtype=torch.cdouble, device=self.device
            )
            self.elec_ampl[self.istate] = amplitudes
        elif amplitudes.shape == (self.nstates, self.ngrids):
            self.elec_ampl = amplitudes
        else:
            raise ValueError("elec_ampl must have shape (ngrids,) or (nstates, ngrids)")

        self.save_every_n_steps = params.get("save_every_n_steps", 1)
        self.properties_to_save = params.get("properties_to_save", ["time", "population_right"])
        self.dt = params.get("dt", 0.01)
        self.nsteps = params.get("nsteps", 500)
        self.ndim = self.nstates * self.ngrids

        energies = torch.as_tensor(
            params.get("E", torch.zeros(self.nstates, self.ngrids)),
            dtype=torch.cdouble, device=self.device,
        )
        if energies.shape != (self.nstates, self.ngrids):
            raise ValueError("E must have shape (nstates, ngrids)")
        if energies.is_complex():
            if torch.any(energies.imag != 0):
                raise ValueError("E must contain real electronic energies")
            energies = energies.real
        self.E = energies.to(dtype=torch.float64)

        if "s_elec" in params:
            self.s_elec = torch.as_tensor(
                params["s_elec"], dtype=torch.cdouble, device=self.device
            )
        else:
            # The same electronic state overlaps with itself at EVERY center.
            self.s_elec = torch.kron(
                torch.eye(self.nstates, dtype=torch.cdouble, device=self.device),
                torch.ones(self.ngrids, self.ngrids, dtype=torch.cdouble, device=self.device),
            )
        if self.s_elec.shape != (self.ndim, self.ndim):
            raise ValueError("s_elec must have shape (nstates * ngrids, nstates * ngrids)")
        self.s_elec = self.s_elec.contiguous()

        # Computed with LDR methods
        self.C0 = torch.zeros(self.ndim, dtype=torch.cdouble, device=self.device)
        self.C_curr = torch.zeros(self.ndim, dtype=torch.cdouble, device=self.device)

        self.s_nucl = torch.eye(self.ngrids, dtype=torch.cdouble, device=self.device)
        self.t_nucl = torch.zeros(self.ngrids, self.ngrids, dtype=torch.cdouble, device=self.device)

        self.S, self.H = torch.zeros(self.ndim, self.ndim, dtype=torch.cdouble, device=self.device), torch.zeros(self.ndim, self.ndim, dtype=torch.cdouble, device=self.device)
        self.U = torch.zeros(self.ndim, self.ndim, dtype=torch.cdouble, device=self.device)
        self.S_half = torch.zeros(self.ndim, self.ndim, dtype=torch.cdouble, device=self.device)
        
        self.time = []
        self.kinetic_energy = []
        self.potential_energy = []
        self.total_energy = []
        self.average_pos = []
        self.population_right = []
        self.denmat = []
        self.norm = []
        self.C_save = []

        self.psi_batch_size = params.get("psi_batch_size", 256)
        if (isinstance(self.psi_batch_size, bool)
                or not isinstance(self.psi_batch_size, int) or self.psi_batch_size < 1):
            raise ValueError("psi_batch_size must be a positive integer")
        if self.eigvecs is None and "psi_r_dia" in self.properties_to_save:
            raise ValueError("Diabatic amplitudes require potential_fn")

        shape = (self.ngrids, self.nstates)
        nsnaps = len(range(0, self.nsteps, self.save_every_n_steps))
        self.psi_r_dia = (torch.zeros(shape, dtype=torch.cdouble, device=self.device)
                          if self.eigvecs is not None else None)
        self.psi_r_adi = torch.zeros(shape, dtype=torch.cdouble, device=self.device)
        self.psi_r_dia_all = (torch.zeros((nsnaps, *shape), dtype=torch.cdouble, device=self.device)
                              if "psi_r_dia" in self.properties_to_save else None)
        self.psi_r_adi_all = (torch.zeros((nsnaps, *shape), dtype=torch.cdouble, device=self.device)
                              if "psi_r_adi" in self.properties_to_save else None)
        self.rho_dia_all = None
        self.rho_adi_all = None

    def update_psi_r(self):
        """Reconstruct both (ngrids, nstates) amplitudes on qgrid.

        Coefficients are state-major: C_curr[i*N+n]. Eigenvector columns
        eigvecs[n,a,i] are <diabatic_a|adiabatic_i(q_n)>. Sum every center
        and reference state in the common diabatic frame, then project onto
        each local adiabatic frame. No renormalization is applied.

        Called during initialization, when saving histories, and at the end
        of propagation. Call explicitly after manually editing C_curr.
        The basis, model, and qgrid remain fixed during a run.
        """
        N, M = self.ngrids, self.nstates
        batch_size = self.psi_batch_size
        if isinstance(batch_size, bool) or not isinstance(batch_size, int) or batch_size < 1:
            raise ValueError("psi_batch_size must be a positive integer")
        coefficients = torch.as_tensor(self.C_curr, dtype=torch.cdouble, device=self.device)
        if coefficients.shape != (M*N,) or not torch.isfinite(coefficients).all():
            raise ValueError("coefficients must be a finite state-major vector of size M*N")
        c = coefficients.reshape(M, N)
        alpha = self._alpha_grid
        prefactor = (2 * alpha / torch.pi).pow(0.25).prod(dim=1)
        if self.eigvecs is None:
            # Precomputed overlaps suffice for the native adiabatic projection.
            # They do not specify a named common diabatic frame.
            if self.s_elec is None:
                raise ValueError("Electronic overlaps are required without reference frames")
            overlaps = self.s_elec.reshape(M, N, M, N)
            psi_dia = None
        else:
            if self.eigvecs.shape != (N, M, M):
                raise ValueError("reference_frames must have shape (N, M, M)")
            center_vectors = torch.einsum("nai,in->na", self.eigvecs, c)
            psi_dia = torch.empty(N, M, dtype=torch.cdouble, device=self.device)
        psi_adi = torch.empty(N, M, dtype=torch.cdouble, device=self.device)
        for start in range(0, N, batch_size):
            stop = min(start + batch_size, N)
            delta = self.qgrid[start:stop, None, :] - self.qgrid[None, :, :]
            # Normalized Gaussian g_n(q_p); temporary arrays are bounded by batch_size.
            kernels = prefactor[None, :] * torch.exp(-(alpha[None, :, :] * delta.square()).sum(-1))
            kernels = kernels.to(torch.cdouble)
            if psi_dia is None:
                psi_adi[start:stop] = torch.einsum(
                    "kpin,pn,in->pk", overlaps[:, start:stop, :, :], kernels, c)
            else:
                psi_dia[start:stop] = kernels @ center_vectors
                psi_adi[start:stop] = torch.einsum(
                    "pak,pa->pk", self.eigvecs[start:stop].conj(), psi_dia[start:stop])
        self.psi_r_dia, self.psi_r_adi = psi_dia, psi_adi
        return self.psi_r_dia, self.psi_r_adi

    def chi_overlap(self):
        """
        Compute nuclear overlap matrix s_nucl[i, j] for the mesh qmesh 
        from the Gaussian basis, g(x; q) = \exp(-\alpha * (x-q)**2).
        """
        delta = self.qgrid[:, None, :] - self.qgrid[None, :, :]    # (N, N, D)
        if self.alpha.ndim < 2:
            exponent = -0.5 * torch.sum(self.alpha * delta**2, dim=2)
            self.s_nucl = torch.exp(exponent)
        else:
            alpha_i = self._alpha_grid[:, None, :]
            alpha_j = self._alpha_grid[None, :, :]
            alpha_sum = alpha_i + alpha_j
            beta = alpha_i * alpha_j / alpha_sum
            prefactor = torch.sqrt(2 * torch.sqrt(alpha_i * alpha_j) / alpha_sum)
            self.s_nucl = prefactor.prod(dim=2) * torch.exp(
                -torch.sum(beta * delta**2, dim=2))

    def chi_kinetic(self):
        """
        Compute nuclear kinetic energy matrix t_nucl[i,j] = <g(x; qgrid[i]) | T | g(x; qgrid[j])>,
        with T = \sum_{\nu} -0.5* m_ν^{-1} \partial^{2}/\partial x_{\nu}^2.
        """
        delta = self.qgrid[:, None, :] - self.qgrid[None, :, :]               # (N, N, D)
        if self.alpha.ndim < 2:
            tau = self.alpha / (2.0 * self.mass) * (1.0 - self.alpha * delta**2)
        else:
            alpha_i = self._alpha_grid[:, None, :]
            alpha_j = self._alpha_grid[None, :, :]
            beta = alpha_i * alpha_j / (alpha_i + alpha_j)
            tau = beta / self.mass * (1.0 - 2.0 * beta * delta**2)
        tau_sum = torch.sum(tau, dim=2)                                       # (N, N)
    
        self.t_nucl = self.s_nucl * tau_sum                                   # (N, N)

    def build_compound_overlap(self):
        """
        Build the compound nuclear-electronic overlap matrix self.S (ndim, ndim)
        """
        N, s, ndim = self.ngrids, self.nstates, self.ndim
    
        # Reshape s_elec[a, b] -> (i, n, j, m) with:
        #   a = i * N + n
        #   b = j * N + m
        s_elec_4d = self.s_elec.view(s, N, s, N) # (i, n, j, m)
    
        s_nucl_4d = self.s_nucl[None, :, None, :] # (1, n, 1, m)
    
        S_4d = s_elec_4d * s_nucl_4d
    
        # Reshape back to (ndim, ndim) with compound indices
        self.S = S_4d.reshape(ndim, ndim)

    def build_compound_hamiltonian(self):
        """Build H = T + V using the selected potential approximation.

        Both nuclear integrals must be built before calling this method.
        Composite indices are state-major: (i,n) -> i*ngrids+n.
        """
        N, s, ndim = self.ngrids, self.nstates, self.ndim
        scheme = self.hamiltonian_scheme
        s_elec_4d = self.s_elec.view(s, N, s, N)
        T_4d = self.t_nucl[None, :, None, :]
        S_4d = self.s_nucl[None, :, None, :]

        if scheme == 'diagonal':
            # Only V is diagonal; retain all kinetic couplings.
            self.H = (s_elec_4d * T_4d).reshape(ndim, ndim)
            self.H.diagonal().add_(self.E.reshape(ndim))
            return
        if scheme == 'midpoint':
            self.H = self._midpoint_hamiltonian(s_elec_4d)
            return
        if scheme == 'as_is':
            # Legacy, generally non-Hermitian form, retained for comparison.
            E_pair = self.E[None, None, :, :]
        elif scheme in ['symmetrized', 'BAT0']:
            E_pair = 0.5 * (self.E[:, :, None, None] + self.E[None, None, :, :])
        elif scheme == 'BAT1':
            E_pair = 0.5 * (self.E[:, :, None, None] + self.E[None, None, :, :])
            gradients = self._reference_energy_gradients()
            # Accumulate one coordinate at a time without allocating an
            # (nstates, N, nstates, N, ndof) intermediate.
            for a in range(self.ndof):
                an = self._alpha_grid[:, None, a]
                am = self._alpha_grid[None, :, a]
                delta = self.qgrid[None, :, a] - self.qgrid[:, None, a]
                # qbar-q_n and qbar-q_m, with center-dependent widths.
                from_bra = am / (an + am) * delta
                from_ket = -an / (an + am) * delta
                E_pair += 0.5 * (
                    gradients[:, :, a][:, :, None, None]
                    * from_bra[None, :, None, :]
                    + gradients[:, :, a][None, None, :, :]
                    * from_ket[None, :, None, :]
                )
        else:
            raise ValueError(f"Unknown Hamiltonian scheme: {scheme}")

        self.H = (s_elec_4d * (T_4d + E_pair * S_4d)).reshape(ndim, ndim)

    def _reference_energy_gradients(self):
        """Validate BAT1 gradients at q_n (energy gradients, not forces)."""
        if self.dE is None:
            raise ValueError("BAT1 requires dE at the reference geometries")
        gradients = torch.as_tensor(self.dE, dtype=torch.cdouble, device=self.device)
        if gradients.shape != (self.nstates, self.ngrids, self.ndof):
            raise ValueError("dE must have shape (nstates, ngrids, ndof)")
        if torch.any(gradients.imag != 0) or not torch.all(torch.isfinite(gradients)):
            raise ValueError("dE must contain finite real energy gradients")
        self.dE = gradients.real.to(dtype=torch.float64)
        return self.dE

    def _midpoint_hamiltonian(self, s_elec_4d):
        """Build the midpoint matrix using the same model as at the centers."""
        if self.potential_fn is None:
            raise ValueError("midpoint requires potential_fn")
        reduced = electronic.compute_midpoint_potential(
            self.qgrid, self._alpha_grid, self.eigvecs, self.potential_fn,
            self.potential_fn_params, self.midpoint_batch_size)
        H = s_elec_4d * self.t_nucl[None, :, None, :]
        H += reduced * self.s_nucl[None, :, None, :]
        return H.reshape(self.ndim, self.ndim)

    def compute_propagator(self):
        """
        Compute U = S^-1/2 exp(-i S^-1/2 H S^-1/2 dt / hbar) S^1/2.

        The overlap must be positive definite; linearly dependent basis
        functions must be removed before propagation.
  
        """
        S = self.S
        H = self.H
        dt = self.dt
    
        evals_S, evecs_S = torch.linalg.eigh(S)
        if not torch.all(torch.isfinite(evals_S) & (evals_S > 0)):
            raise ValueError("The compound overlap S must be positive definite")
        evecs_S = evecs_S.to(dtype=torch.cdouble)
        # eigh returns eigenvectors as COLUMNS of Q: S Q = Q diag(evals_S).
        # Thus S^p = Q diag(evals_S^p) Q^dagger, not Q^dagger diag(...) Q.
        # Broadcasting the eigenvalue vector below scales Q's columns.
        self.S_half = (evecs_S * evals_S.sqrt()) @ evecs_S.conj().T
        S_invhalf = (evecs_S * evals_S.rsqrt()) @ evecs_S.conj().T
    
        H_ortho = S_invhalf @ H @ S_invhalf
    
        evals_H, evecs_H = torch.linalg.eigh(H_ortho)
    
        exp_diag = torch.diag(torch.exp(-1j * evals_H * dt / self.hbar))
        U_ortho = evecs_H @ exp_diag @ evecs_H.conj().T
    
        self.U = S_invhalf @ U_ortho @ self.S_half


    def initialize_C(self):
        """
        Project a normalized Gaussian onto the compound basis and normalize.

        The target nuclear wavefunction is proportional to
        exp[-sum(alpha0 * (q-q0)^2) + i p0.(q-q0)/hbar], where
        alpha0 = sqrt(k * mass)/2. Electronic overlaps are supplied through
        ``elec_ampl`` (see the class docstring). Build S before calling this
        method. Solve S C0 = b with b_a = <psi_a|Psi0>; the overlaps themselves
        are NOT expansion coefficients in a nonorthogonal basis.
        """
        alpha0 = 0.5 * torch.sqrt(self.k * self.mass)
        if torch.any(self.alpha <= 0) or torch.any(alpha0 <= 0):
            raise ValueError("Gaussian widths must be positive")
        delta = self.qgrid - self.q0
        alpha = self._alpha_grid
        width = alpha + alpha0
        momentum = self.p0 / self.hbar
        # Analytic overlap of normalized Gaussians, in coordinates relative
        # to q0 to avoid cancellation between large absolute positions.
        prefactor = torch.prod((4 * alpha * alpha0 / width**2)**0.25, dim=1)
        exponent = torch.sum(
            -alpha * alpha0 / width * delta**2
            -momentum**2 / (4 * width)
            +1j * alpha / width * momentum * delta,
            dim=1,
        )
        nuclear_overlap = prefactor * torch.exp(exponent)
        b = (self.elec_ampl * nuclear_overlap[None, :]).reshape(self.ndim)
        coefficients = torch.linalg.solve(self.S, b)
        norm_squared = torch.vdot(coefficients, self.S @ coefficients).real
        if not torch.isfinite(norm_squared) or norm_squared <= 0:
            raise ValueError("The initial state must have a finite, nonzero projected norm")
        self.C0 = coefficients / norm_squared.sqrt()
        self.C_curr = self.C0.clone()
        self.update_psi_r()
    
    def propagate(self):
        """
        Propagate coefficient.
        """
        # Initialize first step with normalized initial wavefunction
        self.C_curr = self.C0.clone()

        print(F"step = 0")
        self.save_results(0)
        
        for step in range(1, self.nsteps):
            C_vec = self.C_curr.clone()
            self.C_curr = self.U @ C_vec

            if step % self.save_every_n_steps == 0:
                print(F"step = {step}")
                self.save_results(step)
        self.update_psi_r()

    def save_results(self, step):
        save_wavefunction = bool({"psi_r_dia", "psi_r_adi"}.intersection(
            self.properties_to_save))
        if "time" in self.properties_to_save or save_wavefunction:
            self.time.append(step*self.dt)
        if "norm" in self.properties_to_save:
            overlap = torch.matmul(self.S, self.C_curr)
            self.norm.append(torch.sqrt(torch.vdot(self.C_curr, overlap)))
        if "population_right" in self.properties_to_save:
            self.population_right.append(self.compute_populations())
        if "denmat" in self.properties_to_save:
            self.denmat.append(self.compute_denmat())
        if "kinetic_energy" in self.properties_to_save:
            self.kinetic_energy.append(self.compute_kinetic_energy())
        if "potential_energy" in self.properties_to_save:
            self.potential_energy.append(self.compute_potential_energy())
        if "total_energy" in self.properties_to_save:
            self.total_energy.append(self.compute_total_energy())
        if "average_pos" in self.properties_to_save:
            self.average_pos.append(self.compute_average_pos())
        if "C_save" in self.properties_to_save:
            self.C_save.append(self.C_curr)
        if save_wavefunction:
            self.update_psi_r()
            isnap = step // self.save_every_n_steps
            if self.psi_r_dia_all is not None:
                self.psi_r_dia_all[isnap] = self.psi_r_dia
            if self.psi_r_adi_all is not None:
                self.psi_r_adi_all[isnap] = self.psi_r_adi

    def compute_populations(self):
        """Return the diagonal of the Lowdin electronic density matrix.

        Populations use the same orthogonalized coefficients as coherences.
        They are nonnegative and sum to C_curr^+ S C_curr, which is one
        for a normalized wavefunction, without separate renormalization.
        """
        return self.compute_denmat().diagonal().real

    def compute_denmat(self):
        """
        Compute electronic density matrix for a single step using the orthogonalization.
        """
        N, s = self.ngrids, self.nstates
        C_vec = self.C_curr
    
        # Orthogonalize coefficients: C_ortho = S^{1/2} C
        C_ortho = self.S_half @ C_vec
      
        C_blocks = C_ortho.view(s, N)
    
        rho = C_blocks @ C_blocks.conj().T # (s, s)
    
        return rho

    def compute_kinetic_energy(self):
        """
        Compute nuclear kinetic energy as C^+ T C / C^+ S C for a single step.
        """
        N, s, ndim = self.ngrids, self.nstates, self.ndim
    
        # Rebuild compound kinetic matrix: T_4d * s_elec_4d
        s_elec_4d = self.s_elec.view(s, N, s, N)
        T_4d = self.t_nucl[None, :, None, :]
        T_compound = (s_elec_4d * T_4d).reshape(ndim, ndim)
    
        C_vec = self.C_curr
    
        numer = torch.vdot(C_vec, T_compound @ C_vec).real
        denom = torch.vdot(C_vec, self.S @ C_vec).real
    
        return numer / denom
    
    
    def compute_potential_energy(self):
        """Return <V> = <H> - <T> for the Hamiltonian actually built.

        This follows the selected scheme without storing another dense
        compound matrix or repeating expensive midpoint evaluations.
        """
        return self.compute_total_energy() - self.compute_kinetic_energy()
    
    
    def compute_total_energy(self):
        """
        Compute total energy as C^+ H C / C^+ S C for a single step.
        """
        C_vec = self.C_curr
    
        numer = torch.vdot(C_vec, self.H @ C_vec).real
        denom = torch.vdot(C_vec, self.S @ C_vec).real
    
        return numer / denom

    def compute_average_pos(self):
        """
        Compute average position as <q_i> = \sum_i C^+ Q C / C^+ S C for a single step.
        """
        N, s, ndim = self.ngrids, self.nstates, self.ndim
        
        C_vec = self.C_curr

        denom = torch.vdot(C_vec, self.S @ C_vec).real
        s_elec_4d = self.s_elec.view(s, N, s, N)
        
        avg_q = []
        for idof in range(self.ndof):
            alpha_i = self._alpha_grid[:, None, idof]
            alpha_j = self._alpha_grid[None, :, idof]
            q_med = (alpha_i * self.qgrid[:, None, idof]
                     + alpha_j * self.qgrid[None, :, idof]) / (alpha_i + alpha_j)
            q_nucl = self.s_nucl * q_med 
            Q_4d = q_nucl[None, :, None, :]
            Q_4d_compound = s_elec_4d * Q_4d
            Q_compound = Q_4d_compound.reshape(ndim, ndim)

            numer = torch.vdot(C_vec, Q_compound @ C_vec).real
            avg_q.append(numer / denom)

        return avg_q

    def save(self):
        self.update_psi_r()
        torch.save( {"q0":self.q0,
                     "p0":self.p0,
                     "k":self.k,
                     "mass":self.mass,
                     "alpha":self.alpha,
                     "qgrid":self.qgrid,
                     "nstates":self.nstates,
                     "istate":self.istate,
                     "s_nucl":self.s_nucl,
                     "t_nucl":self.t_nucl,
                     "E":self.E,
                     "s_elec":self.s_elec,
                     "eigvecs": self.eigvecs,
                     "psi_r_dia": self.psi_r_dia,
                     "psi_r_adi": self.psi_r_adi,
                     "psi_r_dia_all": self.psi_r_dia_all,
                     "psi_r_adi_all": self.psi_r_adi_all,
                     "rho_dia_all": self.rho_dia_all,
                     "rho_adi_all": self.rho_adi_all,
                     "S":self.S,
                     "H":self.H,
                     "U":self.U,
                     "C_save":self.C_save,
                     "save_every_n_steps":self.save_every_n_steps,
                     "hamiltonian_scheme": self.hamiltonian_scheme,
                     "dE": self.dE,
                     "dt":self.dt, "nsteps":self.nsteps,
                     "time":self.time,
                     "kinetic_energy":self.kinetic_energy,
                     "potential_energy":self.potential_energy,
                     "total_energy":self.total_energy,
                     "average_pos":self.average_pos,
                     "population_right":self.population_right,
                     "denmat":self.denmat,
                     "norm":self.norm
                    }, F"{self.prefix}.pt" )

    def buildSH(self):
        self.chi_overlap()
        self.chi_kinetic()
        self.build_compound_overlap()
        self.build_compound_hamiltonian()
    
    def initialize_operators(self):
        """Build the LDR matrices and propagator, then project the initial state."""
        self.buildSH()
        self.compute_propagator()
        self.initialize_C()

    def solve(self):
        print("Building LDR matrices and initializing the state")
        self.initialize_operators()
        print("Propagating Coefficients")
        self.propagate()
        self.save()


@torch.no_grad()
def wf_calc_nD(filename, axes, potential_fn=None, potential_fn_params=None,
               *, times=None, batch_size=256, normalize=True, device="cpu"):
    """Reconstruct saved LDR amplitudes and densities on a Cartesian display grid.

    Args:
        filename: Path to an LDR ``.pt`` file containing ``qgrid``, ``alpha``,
            ``eigvecs``, ``C_save`` and ``time``. Normalization also needs ``norm``.
        axes: One 1D coordinate array per nuclear dimension. These display
            coordinates are independent of the original Gaussian centers.
        potential_fn: Optional model with the usual LDR signature
            ``potential_fn(q, params)``, where q has shape (ndof, npoints).
            Supply the original model to obtain adiabatic components at the
            display points. Without it, only diabatic and total results exist.
        potential_fn_params: Parameters for the original model.
        times: Saved physical times to select, in the requested order; None
            selects every saved time. Times are not interpolated.
        batch_size: Number of display points evaluated together (default 256).
        normalize: Divide amplitudes by the saved amplitude norm (default True).
            False retains the raw coefficient normalization. The display box
            and individual electronic states are never separately normalized.
        device: Torch device for reconstruction (default "cpu").

    Returns:
        Dictionary with ``axes``, ``time``, ``psi_r_dia``, ``psi_r_adi``,
        ``rho_dia``, ``rho_adi`` and ``rho_total``. Amplitudes and state densities
        have shape (ntime, len(axis_0), ..., len(axis_Dminus1), nstates).
        ``rho_total`` omits the state axis; adiabatic entries are None when
        potential_fn is omitted. Here rho denotes a spatial density, not an
        electronic density matrix. Amplitudes are complex128, densities float64.

    The complex Gaussian contributions are summed before squaring, retaining
    interference. Every Gaussian is evaluated, without density interpolation
    or tail truncation. Projection uses fresh electronic states at each display
    point; individual adiabatic components are basis-dependent at an exact
    degeneracy. The saved reference eigenvectors must define the same complete
    diabatic model basis as potential_fn. No propagation or file output occurs.
    """
    if isinstance(batch_size, bool) or not isinstance(batch_size, int) or batch_size < 1:
        raise ValueError("batch_size must be a positive integer")
    device = torch.device(device)
    # Memory-map large saved matrices; only the reconstruction arrays are read.
    saved = torch.load(filename, map_location="cpu", weights_only=True, mmap=True)
    q = torch.as_tensor(saved["qgrid"], dtype=torch.float64, device=device)
    alpha = torch.as_tensor(saved["alpha"], dtype=torch.float64, device=device).expand_as(q)
    if saved.get("eigvecs") is None:
        raise ValueError("A separate display grid requires saved eigvecs in a common diabatic basis")
    frames = torch.as_tensor(saved["eigvecs"], dtype=torch.cdouble, device=device)
    n, ndof = q.shape
    nstates = frames.shape[-1]
    if frames.shape != (n, nstates, nstates):
        raise ValueError("eigvecs must have shape (ncenters, nstates, nstates)")
    if not torch.isfinite(q).all() or not torch.all(torch.isfinite(alpha) & (alpha > 0)):
        raise ValueError("Centers must be finite and Gaussian widths must be finite and positive")

    saved_times = torch.as_tensor(saved["time"], dtype=torch.float64, device="cpu")
    if saved_times.ndim != 1 or len(saved_times) == 0 or len(saved["C_save"]) != len(saved_times):
        raise ValueError("Save both time and C_save at the same snapshots")
    if times is None:
        indices = torch.arange(len(saved_times))
    else:
        requested = torch.as_tensor(times, dtype=torch.float64, device="cpu")
        if requested.ndim != 1 or requested.numel() == 0 or not torch.isfinite(requested).all():
            raise ValueError("times must be a nonempty sequence of finite saved times")
        indices = (saved_times[:, None] - requested).abs().argmin(dim=0)
        if not torch.allclose(saved_times[indices], requested, rtol=0, atol=1e-8):
            raise ValueError("A requested time is absent from the saved trajectory")
    coefficients = torch.stack([
        torch.as_tensor(saved["C_save"][i], dtype=torch.cdouble, device=device).reshape(nstates, n)
        for i in indices.tolist()
    ])
    if normalize:
        if len(saved.get("norm", [])) != len(saved_times):
            raise ValueError("normalize=True requires a saved norm for every coefficient snapshot")
        # LDR's norm is sqrt(C^dagger S C), not the probability norm itself.
        norms = torch.stack([
            torch.as_tensor(saved["norm"][i], dtype=torch.cdouble, device=device).reshape(())
            for i in indices.tolist()
        ]).abs()
        if not torch.all(torch.isfinite(norms) & (norms > 0)):
            raise ValueError("Saved amplitude norms must be finite and positive")
        coefficients = coefficients / norms[:, None, None]
    if not torch.isfinite(coefficients).all() or not torch.isfinite(frames).all():
        raise ValueError("Saved coefficients and electronic frames must be finite")

    axes = tuple(torch.as_tensor(a, dtype=torch.float64, device=device) for a in axes)
    if len(axes) != ndof or any(a.ndim != 1 or a.numel() == 0 or not torch.isfinite(a).all()
                                for a in axes):
        raise ValueError("axes must contain one nonempty finite 1D array per nuclear dimension")
    points = torch.stack(torch.meshgrid(*axes, indexing="ij"), dim=-1).reshape(-1, ndof)
    ntime = len(indices)
    prefactor = (2 * alpha / torch.pi).pow(0.25).prod(dim=1)
    # U_n[a,i] C[t,i,n] gives each center's contribution in the common frame.
    center_vectors = torch.einsum("nai,tin->nta", frames, coefficients).reshape(n, -1)
    psi_dia = torch.empty(ntime, len(points), nstates, dtype=torch.cdouble, device=device)
    psi_adi = torch.empty_like(psi_dia) if potential_fn is not None else None
    model_params = {} if potential_fn_params is None else potential_fn_params
    for start in range(0, len(points), batch_size):
        stop = min(start + batch_size, len(points))
        p = points[start:stop]
        delta = p[:, None, :] - q[None, :, :]
        gaussian = prefactor * torch.exp(-(alpha * delta.square()).sum(dim=-1))
        psi = (gaussian.to(torch.cdouble) @ center_vectors).reshape(len(p), ntime, nstates)
        psi_dia[:, start:stop] = psi.permute(1, 0, 2)
        if psi_adi is not None:
            _, local_frames = electronic.diagonalize_potential(
                p.T, potential_fn, model_params, nstates)
            chi = torch.einsum("pak,pta->ptk", local_frames.conj(), psi)
            psi_adi[:, start:stop] = chi.permute(1, 0, 2)

    shape = (ntime, *(len(a) for a in axes), nstates)
    psi_dia = psi_dia.reshape(shape)
    psi_adi = psi_adi.reshape(shape) if psi_adi is not None else None
    rho_dia = psi_dia.abs().square()
    rho_adi = psi_adi.abs().square() if psi_adi is not None else None
    return {"axes": axes, "time": saved_times[indices].to(device),
            "psi_r_dia": psi_dia, "psi_r_adi": psi_adi,
            "rho_dia": rho_dia, "rho_adi": rho_adi,
            "rho_total": rho_dia.sum(dim=-1)}
