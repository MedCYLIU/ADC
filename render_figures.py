"""Render six main and two supplementary figures from deposited-data estimates."""

from pathlib import Path
import subprocess, sys, json

base = Path(__file__).resolve().parent
for rel in json.loads((base / "render_steps.json").read_text()):
    print("Rendering", rel, flush=True)
    subprocess.run([sys.executable, str(base / rel)], cwd=base, check=True)
print("Main figures: work/adc_manu_revision/figures")
print("Supplementary figures: work/adc_v6/figures and work/adc_v6/functional")
