from pathlib import Path
import pandas as pd, numpy as np, re, json, sys

sys.stdout.reconfigure(encoding="utf-8")
r = Path(__file__).parent
m = pd.read_excel(r / "JCI176454.sdt1.xlsx")
g = pd.read_csv(r / "GSE243375_All_Samples_TMM_log2_CPM.csv.gz", index_col=0)
print(
    "expression shape",
    g.shape,
    "index unique",
    g.index.is_unique,
    "columnsunique",
    g.columns.is_unique,
)


def norm(s):
    return re.sub("[^a-z0-9]", "", re.sub("^X", "", str(s)).lower())


gcols = {norm(c): c for c in g.columns}
print("norm unique", len(gcols) == len(g.columns))
m["expression_id"] = m["Sample ID"].map(lambda v: gcols.get(norm(v)))
print("missingmatched", m[m.expression_id.isna()].to_dict("records"))
print("extra", sorted(set(g.columns) - set(m.expression_id)))
print(
    "patientconsistency", m.groupby("Patient ID").pCR.nunique().value_counts().to_dict()
)
print(
    "HRconsistency",
    m[m.Treatment_Group == "Pre"]
    .groupby("Patient ID")
    .HR.nunique()
    .value_counts()
    .to_dict(),
)
for gene in ["ERBB2", "TACSTD2", "TOP1", "ABCG2", "ABCC1", "SLC46A3", "SNX10"]:
    m[gene] = m.expression_id.map(g.loc[gene])
m.to_csv(r / "jci_sample_joined.csv", index=False)
b = m[m.Treatment_Group == "Pre"]
p = b.groupby("Patient ID").agg(
    pcr=("pCR", lambda x: int(x.iloc[0] == "Yes")),
    hr=("HR", lambda x: int(x.iloc[0] == "Positive")),
    het=("HER2_Heterogeneity", lambda x: int(x.iloc[0] == "Yes")),
    ERBB2=("ERBB2", "mean"),
    ERBB2_min=("ERBB2", "min"),
    ERBB2_max=("ERBB2", "max"),
    ncores=("ERBB2", "size"),
    FISH=("ERBB2_FISH_Amp%", "mean"),
)
p["range"] = p.ERBB2_max - p.ERBB2_min
q = pd.read_excel(r / "JCI176454.sdval.xlsx", sheet_name="Fig.1D").rename(
    columns={"PatientID": "Patient ID", "ERBB2": "reported_ERBB2"}
)
p = p.join(q.set_index("Patient ID"))
print("diff vs reported mean", np.max(np.abs(p.ERBB2 - p.reported_ERBB2)))
p.to_csv(r / "jci_patient_baseline.csv")
print(p.groupby(["pcr", "hr"]).size().to_dict())
print(
    "pairedcores",
    (p.ncores == 2).sum(),
    "range",
    p.loc[p.ncores == 2, "range"].describe().to_dict(),
)
print("Nhet", p.het.sum())
(r / "linkage_audit.json").write_text(
    json.dumps(
        {
            "matrix_shape": list(g.shape),
            "clinical_sample_n": len(m),
            "unmatched_clinical": m.loc[m.expression_id.isna(), "Sample ID"].tolist(),
            "unannotated_expression": sorted(set(g.columns) - set(m.expression_id)),
            "baseline_patients": len(p),
            "baseline_samples": len(b),
            "two_core_patients": int((p.ncores == 2).sum()),
            "post_samples": int((m.Treatment_Group == "Post").sum()),
            "maximum_ERBB2_discrepancy_from_published_values": float(
                np.max(np.abs(p.ERBB2 - p.reported_ERBB2))
            ),
        },
        indent=2,
    )
)
