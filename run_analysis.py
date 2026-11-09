"""Run the deposited-data analyses with their original numerical code."""

from pathlib import Path
import subprocess, sys, json

base = Path(__file__).resolve().parent
subprocess.run([sys.executable, str(base / "download_data.py"), "--check"], check=True)
for rel in json.loads((base / "analysis_steps.json").read_text()):
    script = base / rel
    args = []
    if script.name in ["audit_sampling.py", "cellularity.py"]:
        args = [
            str(base / "work/adc/source.xlsx"),
            "--outdir",
            str(base / "work/adc_v2/audit"),
        ]
    print("Running", rel, flush=True)
    subprocess.run([sys.executable, str(script), *args], cwd=script.parent, check=True)
print("Analysis complete. Run render_figures.py to render the submitted figures.")
