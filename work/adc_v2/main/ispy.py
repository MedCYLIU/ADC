"""Secondary analysis of I-SPY2 GSE181574. Run in the directory with the downloaded files.
Reads NCBI processed gene values and series annotations, and verifies against article Supplementary Data 1.
"""

import csv, gzip, hashlib, json, math, re
from pathlib import Path
from statistics import NormalDist
import numpy as np


def savecsv(path, data):
    with open(path, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)


specimens = {}
for matrix in sorted(Path(".").glob("GSE181574-*_series_matrix.txt.gz")):
    annotation = {}
    with gzip.open(matrix, "rt", encoding="utf-8") as inp:
        for row in csv.reader(inp, delimiter="\t"):
            if not row:
                continue
            tag = row[0]
            if tag == "!series_matrix_table_begin":
                break
            if tag == "!Sample_title":
                annotation["id"] = [
                    re.search(r"(\d+)\s*$", s).group(1) for s in row[1:]
                ]
            elif tag == "!Sample_geo_accession":
                annotation["gsm"] = row[1:]
            elif tag == "!Sample_characteristics_ch1":
                fields = [s.split(": ", 1) for s in row[1:]]
                assert len(set(f[0] for f in fields)) == 1
                annotation[fields[0][0]] = [f[1] for f in fields]
    for i, key in enumerate(annotation["id"]):
        assert key not in specimens
        specimens[key] = {k: v[i] for k, v in annotation.items()}
        specimens[key]["platform"] = matrix.name.split("-")[1].split("_")[0]
with gzip.open("gene.txt.gz", "rt", encoding="utf-8") as inp:
    records = csv.reader(inp, delimiter="\t")
    names = next(records)
    assert len(set(names)) == len(names) == 127 and set(names) == set(specimens)
    for row in records:
        if row[0] in ("ERBB2", "TACSTD2"):
            for key, v in zip(names, row[1:], strict=True):
                specimens[key][row[0]] = float(v)
ann = list(
    csv.DictReader(open("41467_2021_26019_MOESM4_ESM.csv", encoding="utf-8-sig"))
)
ind = {r["ResearchID"]: r for r in ann}
assert len(ind) == len(ann)
omap = {
    "T-DM1/P": "T-DM1 + Pertuzumab",
    "THP": "Paclitaxel + Pertuzumab + Trastuzumab",
    "TH Control": "Paclitaxel + Trastuzumab",
}
assert len({r["gsm"] for r in specimens.values()}) == 127
for key, row in specimens.items():
    orig = ind[key]
    assert (
        row["her2"] == "1"
        and row["pcr"] == orig["pCR (0:no pCR; 1:pCR)"]
        and row["hr"] == orig["HR status (0:HR-;1:HR+)"]
    )
    assert orig["Arm"].lower().replace(" ", "") == omap[row["arm"]].lower().replace(
        " ", ""
    )
rows = list(specimens.values())
savecsv("ispy_joined.csv", rows)
expr = np.array([r["ERBB2"] for r in rows])
assert np.isfinite(expr).all()
mu, sd = expr.mean(), expr.std(ddof=1)


def design(arm):
    armdata = [r for r in rows if r["arm"] == arm]
    y = np.array([int(r["pcr"]) for r in armdata])
    z = np.array([(r["ERBB2"] - mu) / sd for r in armdata])
    hr = np.array([int(r["hr"]) for r in armdata])
    return y, z, hr


def logistic(X, y, offset=0):
    b = np.zeros(X.shape[1])
    for _ in range(100):
        eta = offset + X @ b
        prob = np.exp(-np.logaddexp(0, -eta))
        score = X.T @ (y - prob)
        inf = X.T @ (X * (prob * (1 - prob))[:, None])
        if max(abs(score)) < 1e-9:
            break
        step = np.linalg.solve(inf, score)
        b += step
    assert max(abs(score)) < 1e-7
    eta = offset + X @ b
    return b, np.sum(y * eta - np.logaddexp(0, eta))


def regression(arm):
    y, z, h = design(arm)
    base = np.stack([np.ones(len(y)), h], 1)
    X = np.column_stack([base, z])
    beta, ml = logistic(X, y)

    def lr(b):
        _, ll = logistic(base, y, b * z)
        return 2 * (ml - ll)

    roots = []
    for direction in [-1, 1]:
        a = beta[-1]
        b = a + 0.5 * direction
        while lr(b) < 3.841458820694124:
            b = a + 2 * (b - a)
        for _ in range(50):
            mid = (a + b) / 2
            if lr(mid) < 3.841458820694124:
                a = mid
            else:
                b = mid
        roots.append(math.exp((a + b) / 2))
    _, null = logistic(base, y)
    p = math.erfc(math.sqrt(max(0, 2 * (ml - null)) / 2))
    return dict(
        arm=arm,
        n=len(y),
        pcr=int(y.sum()),
        non_pcr=int(len(y) - y.sum()),
        OR=math.exp(beta[-1]),
        lo=roots[0],
        hi=roots[1],
        p_lr=p,
        mean_ERBB2_pcr=np.mean(z[y == 1]) * sd + mu,
        mean_ERBB2_residual=np.mean(z[y == 0]) * sd + mu,
    )


r = [regression(a) for a in ["T-DM1/P", "THP", "TH Control"]]
savecsv("ispy_models.csv", r)
metadata = dict(
    publication_n=len(ann),
    matched_n=len(rows),
    missing_from_geo=sorted(set(ind) - set(specimens)),
    sd_ERBB2=sd,
    mean_ERBB2=mu,
    input_sha256={
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [
            Path("gene.txt.gz"),
            *sorted(Path(".").glob("GSE181574-*_series_matrix.txt.gz")),
            Path("41467_2021_26019_MOESM4_ESM.csv"),
        ]
    },
    models=r,
)
Path("ispy_stats.json").write_text(json.dumps(metadata, indent=2))
print(json.dumps(metadata, indent=2))
