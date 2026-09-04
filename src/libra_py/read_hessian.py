# *********************************************************************************
# * Copyright (C) 2026 Somesh Chandra
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
# Reads the Hessian from different packages
# Covert Hessian to numpy array as an input for wigner sampling
import numpy as np
import h5py

"""
read_hessian.py module provides a unified interface for reading Hessian matrices (second 
derivatives of energy with respect to nuclear coordinates) from three popular 
quantum chemistry packages:

  1. DFTB+ (semi-empirical tight-binding)
  2. OpenMolcas (ab initio multi-configuration methods)
  3. PySCF (Python-based ab initio package)

All functions return a symmetric (3N, 3N) numpy array in atomic units 
(Hartree/Bohr²), making them directly compatible with Libra molecular dynamics 
and normal mode analysis.

================================================================================
USAGE QUICK START
================================================================================

from read_hessian import read_openmolcas
from read_hessian import read_dftbplus
from read_hessian import read_pyscf

or 
from read_hessian import read_hessian
source = dftb+/ openmolcas/ pyscf

# DFTB+ (from hessian.out file)
H = read_hessian('dftb+', filename='hessian.out', natoms=None)

# OpenMolcas (from .slapaf.h5 HDF5 file)
H = read_hessian('openmolcas', filename='molecule.slapaf.h5')

# PySCF (from computed mf.Hessian().kernel() output)
from pyscf import gto, dft
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
mf = dft.RKS(mol)
mf.xc = 'b3lyp'
mf.kernel()
H = read_hessian('pyscf', hessian=mf.Hessian().kernel())

# All return the same format: (3N, 3N) symmetric matrix
print(H.shape)  # (6, 6) for H2
print(np.allclose(H, H.T))  # True (symmetric)


All functions return a numpy array of shape (3N, 3N) in Hartree/Bohr^2.
"""

def read_dftbplus(filename="hessian.out", natoms=None):
    """
    Reads DFTB+ hessian.out file.
    Automatically detects if it is in packed lower-triangular format or full matrix format.
    
    Args:
        filename: Path to hessian.out file
        natoms: Number of atoms (required for full matrix detection)
    
    Returns:
        A symmetric (3N, 3N) numpy array
    """
    # 1. Read all numbers
    numbers = []
    with open(filename, 'r') as f:
        for line in f:
            vals = line.strip().split()
            for v in vals:
                try:
                    numbers.append(float(v))
                except ValueError:
                    pass  # skip headers/text
    
    total = len(numbers)
    
    # 2. Determine shape based on total numbers
    if natoms is not None:
        ndof = 3 * natoms
        expected_packed = ndof * (ndof + 1) // 2
        expected_full = ndof * ndof
        
        if total == expected_packed:
            # Reconstruct from lower-triangular packed
            hessian = np.zeros((ndof, ndof))
            idx = 0
            for i in range(ndof):
                for j in range(i + 1):
                    hessian[i, j] = numbers[idx]
                    hessian[j, i] = numbers[idx]
                    idx += 1
            return hessian
            
        elif total == expected_full:
            # Reconstruct from full matrix
            hessian = np.array(numbers).reshape((ndof, ndof))
            # Ensure it is symmetric
            return (hessian + hessian.T) / 2.0
            
        else:
            raise ValueError(f"Found {total} numbers. For {natoms} atoms, expected {expected_packed} (packed) or {expected_full} (full).")
            
    else:
        raise ValueError("natoms must be provided to correctly parse the DFTB+ Hessian.")


def read_openmolcas(filename):
    """
    Reads Hessian from OpenMolcas .slapaf.h5 file.

    Supports either:
      - packed lower-triangular 1D format
      - full 2D square matrix format

    Returns:
        A symmetric (3N, 3N) numpy array
    """
    with h5py.File(filename, 'r') as f:
        if 'HESSIAN' in f:
            raw = np.array(f['HESSIAN'][:])
        elif 'Hessian' in f:
            raw = np.array(f['Hessian'][:])
        elif 'hessian' in f:
            raw = np.array(f['hessian'][:])
        else:
            available_keys = list(f.keys())
            raise KeyError(f"Hessian not found. Available keys: {available_keys}")

        if 'COORDINATES' in f:
            q_eq = np.array(f['COORDINATES'][0]).flatten()
            ndof = len(q_eq)
        else:
            raise KeyError("COORDINATES dataset not found; cannot infer Hessian dimension")

    # Case 1: already a full square matrix
    if raw.ndim == 2:
        hessian = raw

    # Case 2: packed lower-triangular 1D array
    elif raw.ndim == 1:
        expected = ndof * (ndof + 1) // 2
        if raw.size != expected:
            raise ValueError(
                f"Packed Hessian has {raw.size} elements, expected {expected} for ndof={ndof}"
            )

        hessian = np.zeros((ndof, ndof))
        idx = 0
        for i in range(ndof):
            for j in range(i + 1):
                hessian[i, j] = raw[idx]
                hessian[j, i] = raw[idx]
                idx += 1
    else:
        raise ValueError(f"Unsupported Hessian shape: {raw.shape}")

    # Ensure symmetry
    if not np.allclose(hessian, hessian.T, rtol=1e-8, atol=1e-10):
        hessian = 0.5 * (hessian + hessian.T)

    return hessian

def read_pyscf(hess_pyscf):
    """
    Converts PySCF Hessian from (N, N, 3, 3) to (3N, 3N).

    Args:
        hess_pyscf: Output from mf.Hessian().kernel()

    Returns:
        A symmetric (3N, 3N) numpy array
    """
    hess_pyscf = np.asarray(hess_pyscf)

    if hess_pyscf.ndim != 4:
        raise ValueError(
            f"PySCF Hessian must have 4 dimensions (N, N, 3, 3), got shape {hess_pyscf.shape}"
        )

    if hess_pyscf.shape[2:] != (3, 3):
        raise ValueError(
            f"Last two dimensions of PySCF Hessian must be (3, 3), got {hess_pyscf.shape[2:]}"
        )

    if hess_pyscf.shape[0] != hess_pyscf.shape[1]:
        raise ValueError(
            f"First two dimensions of PySCF Hessian must match, got {hess_pyscf.shape[:2]}"
        )

    N = hess_pyscf.shape[0]
    ndof = 3 * N
    hessian = hess_pyscf.transpose(0, 2, 1, 3).reshape(ndof, ndof)  # ← Changed

    if not np.allclose(hessian, hessian.T, rtol=1e-8, atol=1e-10):
        hessian = 0.5 * (hessian + hessian.T)

    return hessian

# someone want to import by this method
#from read_hessian import read_hessian
#source = dftb+
#source = openmolcas
#source =  pyscf
def read_hessian(source, **kwargs):
    """
    Unified interface for reading Hessians from different packages.
    
    Usage:
        # DFTB+
        H = read_hessian('dftb+', filename='hessian.out')
        
        # OpenMolcas  
        H = read_hessian('openmolcas', filename='molecule.slapaf.h5')
        
        # PySCF (requires already computed hessian)
        H = read_hessian('pyscf', hessian=mf.Hessian().kernel())
    
    Args:
        source: 'dftb+', 'openmolcas', or 'pyscf'
        **kwargs: Arguments specific to each source
        
    Returns:
        A (3N, 3N) numpy array
    """
    source = source.lower()
    
    if source in ['dftb+', 'dftb', 'dftbplus']:
        return read_dftbplus(kwargs.get('filename', 'hessian.out'), 
                            kwargs.get('natoms', None))
    
    elif source in ['openmolcas', 'molcas']:
        return read_openmolcas(kwargs.get('filename'))
    
    elif source in ['pyscf', 'py']:
        return read_pyscf(kwargs.get('hessian'))
    
    else:
        raise ValueError(f"Unknown source: {source}. Use 'dftb+', 'openmolcas', or 'pyscf'")
