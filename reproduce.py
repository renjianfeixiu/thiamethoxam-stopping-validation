from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
for folder, script in [('baseline','run.py'), ('extension','experiments.py'), ('extension','verify.py'), ('external','validate.py'), ('external','check.py')]:
    subprocess.run([sys.executable, script], cwd=root / folder, check=True)
