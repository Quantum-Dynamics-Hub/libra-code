from contextlib import redirect_stdout
from io import StringIO
import os
from pathlib import Path
import runpy
import shutil
import sys


ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "src"
EXAMPLES = ROOT / "examples" / "libra_py" / "dyn" / "propagation"
SAVING_EXAMPLES = EXAMPLES / "ehrenfest_tully1_saving"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def test_propagation_examples_run():
    for name in (
        "example_nuclear.py",
        "example_electronic.py",
        "example_coupled.py",
        "example_engine.py",
        "example_ehrenfest_tully1_ensemble.py",
    ):
        stream = StringIO()
        old_nsteps = os.environ.get("EHRENFEST_TULLY1_NSTEPS")
        if name == "example_ehrenfest_tully1_ensemble.py":
            os.environ["EHRENFEST_TULLY1_NSTEPS"] = "30"
        try:
            with redirect_stdout(stream):
                runpy.run_path(str(EXAMPLES / name), run_name="__main__")
            assert stream.getvalue()
        finally:
            if old_nsteps is None:
                os.environ.pop("EHRENFEST_TULLY1_NSTEPS", None)
            else:
                os.environ["EHRENFEST_TULLY1_NSTEPS"] = old_nsteps
            shutil.rmtree(ROOT / "ehrenfest_tully1_outputs", ignore_errors=True)


def test_ehrenfest_tully1_saving_examples_run():
    old_nsteps = os.environ.get("EHRENFEST_TULLY1_NSTEPS")
    old_path = list(sys.path)
    os.environ["EHRENFEST_TULLY1_NSTEPS"] = "5"
    sys.path.insert(0, str(SAVING_EXAMPLES))
    try:
        for name in (
            "01_hdf5_stride.py",
            "02_hdf5_less_frequent_flush.py",
            "03_hdf5_chunked_timeseries.py",
            "04_hdf5_fewer_fields.py",
            "05_npz_saver.py",
            "06_json_saver.py",
        ):
            stream = StringIO()
            with redirect_stdout(stream):
                runpy.run_path(str(SAVING_EXAMPLES / name), run_name="__main__")
            assert stream.getvalue()
    finally:
        if old_nsteps is None:
            os.environ.pop("EHRENFEST_TULLY1_NSTEPS", None)
        else:
            os.environ["EHRENFEST_TULLY1_NSTEPS"] = old_nsteps
        sys.path[:] = old_path
        for path in ROOT.glob("*_outputs"):
            if path.name.startswith(("01_hdf5_", "02_hdf5_", "03_hdf5_", "04_hdf5_", "05_npz_", "06_json_")):
                shutil.rmtree(path, ignore_errors=True)
        shutil.rmtree(ROOT / ".ehrenfest_tully1_saving_cache", ignore_errors=True)
