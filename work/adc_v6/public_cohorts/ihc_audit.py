import pandas as pd, sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
r = Path(__file__).parent
p = pd.read_excel(r / "JCI176454.sdt3.xlsx")
print(p.to_string(index=False))
d = pd.read_csv(r / "jci_paired_change.csv")
print("RNA IDs notinIHC", set(d["Patient ID"]) - set(p["Patient ID"]))
m = d.merge(p, on="Patient ID", how="left")
print("merged change", m["HER2 IHC Change"].value_counts(dropna=False).to_dict())
m.to_csv(r / "jci_paired_rna_ihc.csv", index=False)
