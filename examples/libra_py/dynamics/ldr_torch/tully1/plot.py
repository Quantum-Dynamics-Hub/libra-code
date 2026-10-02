"""Plot the local LDR and DVR results together. No arguments or environment setup required."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _plotting import plot_cases
from settings import CASES, OUTPUT_DIR


if __name__ == "__main__":
    plot_cases(CASES, OUTPUT_DIR)
