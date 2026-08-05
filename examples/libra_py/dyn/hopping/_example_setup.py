"""Directory-local setup shared by the hopping examples."""

from pathlib import Path
import sys


EXAMPLE_DIR = Path(__file__).resolve().parent


def use_repo_sources():
    """Make ``src`` importable when launched from this examples directory."""

    if Path.cwd().resolve() != EXAMPLE_DIR:
        raise SystemExit(
            "Run this example from examples/libra_py/dyn/hopping, for example:\n"
            "  cd examples/libra_py/dyn/hopping\n"
            "  python example_direct_probabilities.py"
        )
    source = EXAMPLE_DIR.parents[3] / "src"
    sys.path.insert(0, str(source))
