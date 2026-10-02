# 2D FLV model

From this directory:

```bash
python run.py
python plot.py
```

The default run covers 3000 a.u. for both weak (gamma=0.01) and strong (gamma=0.08) coupling, with LDR and DVR.
No environment variables or command-line parameters are required.
Edit `settings.py` to change the parameters or the local `OUTPUT_DIR`.
Results and PNG/PDF comparison figures are written under `output/`.

For a short execution check in separate `output/quick/` directories:

```bash
python run.py --quick
python plot.py --quick
```

The manuscript LDR grid has 7381 Gaussian centers. Its dense compound matrices require tens of GB of working memory. Use quick mode to check setup on a smaller machine. The two coupling cases have separate `output/weak/` and `output/strong/` results and density figures; `output/comparison.png` and `.pdf` compare both cases side by side.

See the [shared instructions](../README.md) for physical conventions, numerical methods, output fields, and optional separate LDR/DVR runs.
