"""Retrieve exact public inputs, or import them from the supplied local archive."""

from pathlib import Path
from zipfile import ZipFile
from urllib.request import Request, urlopen
import argparse, hashlib, json, shutil

BASE = Path(__file__).resolve().parent


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def validate(path, item):
    return (
        path.is_file()
        and path.stat().st_size == item["bytes"]
        and sha(path) == item["sha256"]
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--from-archive",
        type=Path,
        help="Import pinned inputs from ADC_public_data_complete_package.zip instead of downloading",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Check existing inputs without changing files or using the network",
    )
    args = parser.parse_args()
    items = json.loads((BASE / "data_sources.json").read_text())
    archive = ZipFile(args.from_archive) if args.from_archive else None
    try:
        for item in items:
            dest = (BASE / item["path"]).resolve()
            if not dest.is_relative_to(BASE):
                raise ValueError("Input path is outside repository")
            if validate(dest, item):
                print("Verified", item["path"], flush=True)
                continue
            if args.check:
                raise SystemExit("Missing or changed input: " + item["path"])
            dest.parent.mkdir(parents=True, exist_ok=True)
            tmp = dest.with_name(dest.name + ".download")
            if archive:
                with archive.open(item["path"]) as src, tmp.open("wb") as out:
                    shutil.copyfileobj(src, out)
            else:
                print("Downloading", item["path"], flush=True)
                req = Request(
                    item["url"], headers={"User-Agent": "ADC-reproducibility/1.0"}
                )
                with urlopen(req, timeout=60) as src, tmp.open("wb") as out:
                    shutil.copyfileobj(src, out)
            if not validate(tmp, item):
                raise SystemExit(
                    "Checksum or size mismatch: "
                    + item["path"]
                    + ". The expected snapshot was not replaced. Use the pinned local archive or inspect source changes; do not alter expected hashes silently."
                )
            tmp.replace(dest)
            print("Verified", item["path"], flush=True)
    finally:
        if archive:
            archive.close()


if __name__ == "__main__":
    main()
