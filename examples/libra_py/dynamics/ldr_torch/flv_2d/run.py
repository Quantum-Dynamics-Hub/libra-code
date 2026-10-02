"""Run both LDR and DVR with the local settings. No arguments or environment setup required."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import run_cases
from settings import CASES, OUTPUT_DIR


if __name__ == "__main__":
    run_cases(CASES, OUTPUT_DIR)
