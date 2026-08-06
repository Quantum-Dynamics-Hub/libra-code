from contextlib import contextmanager, redirect_stdout
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
TEST_OUTPUT = Path(__file__).resolve().parent / "output"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


@contextmanager
def _example_workdir(name):
    """Run an example in a test-local directory and always remove artifacts."""

    run_dir = TEST_OUTPUT / name
    run_dir.mkdir(parents=True, exist_ok=True)
    previous_cwd = Path.cwd()
    os.chdir(run_dir)
    try:
        yield run_dir
    finally:
        os.chdir(previous_cwd)
        shutil.rmtree(run_dir, ignore_errors=True)
        if TEST_OUTPUT.exists() and not any(TEST_OUTPUT.iterdir()):
            TEST_OUTPUT.rmdir()


def test_propagation_examples_run():
    with _example_workdir("propagation_examples"):
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


def test_ehrenfest_tully1_saving_examples_run():
    with _example_workdir("saving_examples"):
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
