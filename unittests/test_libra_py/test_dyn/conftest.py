"""Make dyn source imports work from individual test directories."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT/'src') not in sys.path:
    sys.path.insert(0, str(ROOT/'src'))

sys.path.insert(0, str(Path(__file__).resolve().parent))
