"""Independent coordinate maximization and direct likelihood checks."""

from pathlib import Path
import pandas as pd, numpy as np, json, math

r = Path(__file__).parent
p = pd.read_csv(r / "jci_patient_baseline.csv")
res = json.loads((r / "jci_statistics.json").read_text())
mu = res["mean_ERBB2"]
sd = res["SD_ERBB2"]
q = p[p.ncores == 2]


def coordinate(X, y, offset=0):
    b = np.zeros(X.shape[1])
    last = -np.inf
    for k in range(4000):
        for j in range(X.shape[1]):
            partial = offset + X @ b - X[:, j] * b[j]
            col = X[:, j]

            def score(v):
                z = partial + v * col
                pr = np.exp(-np.logaddexp(0, -z))
                return col @ (y - pr)

            lo, hi = -1.0, 1.0
            while score(lo) < 0:
                lo *= 2
            while score(hi) > 0:
                hi *= 2
            for _ in range(50):
                mid = (lo + hi) / 2
                if score(mid) > 0:
                    lo = mid
                else:
                    hi = mid
            b[j] = (lo + hi) / 2
        eta = offset + X @ b
        ll = float(np.sum(y * eta - np.logaddexp(0, eta)))
        if abs(ll - last) < 1e-12:
            break
        last = ll
    return b, ll, k


out = []
for i, (data, col) in enumerate(
    [(p, "ERBB2"), (q, "ERBB2"), (q, "ERBB2_min"), (q, "ERBB2_max")]
):
    y = data.pcr.to_numpy()
    z = (data[col].to_numpy() - mu) / sd
    X = np.column_stack([np.ones(len(y)), data.hr, z])
    b, ll, k = coordinate(X, y)
    expected = res["models"][i]
    checks = []
    for bound in ["lo", "hi"]:
        offset = np.log(expected[bound]) * z
        _, lp, _ = coordinate(X[:, :2], y, offset)
        checks.append(2 * (ll - lp))
    diff = abs(float(np.exp(b[-1])) - expected["OR"])
    assert diff < 1e-5
    assert max(abs(np.array(checks) - 3.841458820694124)) < 1e-8
    out.append(
        {
            "model": expected["model"],
            "independent_OR": float(np.exp(b[-1])),
            "absolute_OR_difference": diff,
            "coordinate_iterations": k + 1,
            "profile_endpoint_LR": checks,
        }
    )
# Source sanity assertions separate from model code.
s = pd.read_excel(r / "JCI176454.sdt1.xlsx")
pre = s[s.Treatment_Group == "Pre"]
assert len(pre) == 242 and pre["Patient ID"].nunique() == 129
assert (pre.groupby("Patient ID").pCR.nunique() == 1).all()
assert (pre.groupby("Patient ID").HR.nunique() == 1).all()
assert (pre.groupby("Patient ID").HER2_Heterogeneity.nunique() == 1).all()
(r / "independent_model_audit.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
