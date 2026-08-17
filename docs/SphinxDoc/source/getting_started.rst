Getting started
===============

Libra's Python layer lives in ``src/libra_py`` and uses the compiled
``liblibra_core`` extension for performance-sensitive algorithms and data
structures. From a configured development checkout, verify both layers with:

.. code-block:: bash

   python -c "import libra_py, liblibra_core; print('Libra is available')"

Using the API reference
-----------------------

The :doc:`reference/python_api` is organized by Python package and module.
Each module page lists its public functions, classes, arguments, return values,
and source links when these are present in its docstrings.

The :doc:`reference/boost_python` covers symbols exported from the C++ core.
These objects are imported from ``liblibra_core`` on Linux and macOS. Some
historical Windows environments use the compatible ``cyglibra_core`` name.

Building these docs
-------------------

.. code-block:: bash

   cd docs/SphinxDoc
   python -m pip install -r requirements.txt
   make html

The local website is written to ``build/html/index.html``. The recursive
Python API stubs are regenerated automatically at the beginning of every
Sphinx build, so newly added modules appear without manual toctree edits.

