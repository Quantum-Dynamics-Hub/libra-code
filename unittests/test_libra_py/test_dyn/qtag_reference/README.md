# Frozen standalone QTAG oracle

This is the former `src/libra_py/dyn/methods/qtag` implementation, retained only
for regression comparisons while production QTAG now specializes DynamicsEngine.
Only imports of the shared GWP integrals were made absolute when relocating it.
It is not installed with libra_py and must not be imported by production code.
Test conftest adds this directory's parent to the test import path. Existing
formula/C++ tests continue to exercise it, and engine tests compare histories
and spatial densities against it. No second production QTAG driver is retained.
