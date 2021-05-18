# Computing the molecular overlaps using Libint and CP2K

Here are 3 python test files for computing both the atomic and molecular orbitals (MOs) overlap matrices. In computing the MOs overlap, we need to make sure that the 
diagonal elements of the MO overlap matrix are all `1`. 

## `test.py` 

This file tests the libint library in Libra for computing the atomic orbital overlap matrix.

## `test_molden.py`

This file computes the MO overlaps from a `molden` file format printed out by CP2K. 

## `test_molog.py`

This file computes the MO overlaps from a `MOLog` file printed out by CP2K.
