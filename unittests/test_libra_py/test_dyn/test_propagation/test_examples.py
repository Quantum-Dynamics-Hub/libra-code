from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import runpy
import sys


ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "src"
EXAMPLES = ROOT / "examples" / "libra_py" / "dyn" / "propagation"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def test_propagation_examples_run():
    for name in (
        "example_nuclear.py",
        "example_electronic.py",
        "example_coupled.py",
        "example_engine.py",
    ):
        stream = StringIO()
        with redirect_stdout(stream):
            runpy.run_path(str(EXAMPLES / name), run_name="__main__")
        assert stream.getvalue()
