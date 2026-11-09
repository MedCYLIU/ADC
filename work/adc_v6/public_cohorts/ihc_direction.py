from pathlib import Path
import pandas as pd, numpy as np, json

r = Path(__file__).parent
d = pd.read_csv(r / "jci_paired_rna_ihc.csv")
d["RNA_direction"] = np.where(
    d.delta > 0, "Increased", np.where(d.delta < 0, "Decreased", "Unchanged")
)
t = pd.crosstab(d["HER2 IHC Change"], d.RNA_direction)
print(t.to_string())
print(
    d.groupby("HER2 IHC Change")
    .delta.agg(["count", "median", "min", "max"])
    .to_string()
)
t.to_csv(r / "rna_ihc_direction_crosstab.csv")
a = d[d["HER2 IHC Change"].isin(["Increased", "Decreased"])]
print("opposite", sum(a["HER2 IHC Change"] != a.RNA_direction), len(a))
print(d[d["HER2 IHC Change"] == "Not Evaluated"]["Patient ID"].tolist())
