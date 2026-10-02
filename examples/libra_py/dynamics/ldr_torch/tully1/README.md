# Tully 1: simple avoided crossing

From this directory:

```bash
python run.py
python plot.py
```

The default run covers 1300 a.u. with 841 grid points, dt=1, and both LDR and DVR.
No environment variables or command-line parameters are required.
Edit `settings.py` to change the parameters or the local `OUTPUT_DIR`.
Results and PNG/PDF comparison figures are written under `output/`.

For a short execution check in separate `output/quick/` directories:

```bash
python run.py --quick
python plot.py --quick
```

See the [shared instructions](../README.md) for physical conventions, numerical methods, output fields, and optional separate LDR/DVR runs.
