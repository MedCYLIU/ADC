"""Download and validate the two official GSE209998 source files."""

from concurrent.futures import ThreadPoolExecutor
import gzip
from pathlib import Path
import time
from urllib.request import Request, urlopen

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE209nnn/GSE209998/"
MATRIX = "GSE209998_AUR_129_UQN.txt.gz"
META = "GSE209998_series_matrix.txt.gz"


def validate_gzip(path):
    """Read the whole stream to check for truncation and CRC errors."""
    with gzip.open(path, "rb") as handle:
        while handle.read(4 * 1024 * 1024):
            pass


def download():
    sources = [(GEO + "matrix/" + META, META), (GEO + "suppl/" + MATRIX, MATRIX)]
    for url, filename in sources:
        destination = Path(filename)
        temporary = Path(filename + ".tmp")
        if destination.exists():
            try:
                validate_gzip(destination)
                continue
            except (OSError, EOFError):
                pass

        request = Request(url, headers={"Range": "bytes=0-0"})
        with urlopen(request, timeout=35) as response:
            partial_content = response.status == 206
            if partial_content:
                total_size = int(response.headers["Content-Range"].split("/")[-1])
            else:
                with temporary.open("wb") as output:
                    while block := response.read(65536):
                        output.write(block)

        if partial_content:
            chunk_size = 1024 * 1024

            def request_chunk(index):
                start = index * chunk_size
                end = min((index + 1) * chunk_size, total_size) - 1
                for attempt in range(3):
                    try:
                        request = Request(
                            url, headers={"Range": f"bytes={start}-{end}"}
                        )
                        with urlopen(request, timeout=105) as response:
                            assert response.status == 206
                            expected_range = f"bytes {start}-{end}/{total_size}"
                            assert response.headers["Content-Range"] == expected_range
                            block = response.read()
                            assert len(block) == end - start + 1
                            return block
                    except Exception:
                        if attempt == 2:
                            raise
                        time.sleep(1 + attempt)

            chunk_count = (total_size + chunk_size - 1) // chunk_size
            with ThreadPoolExecutor(max_workers=12) as pool, temporary.open(
                "wb"
            ) as output:
                for block in pool.map(request_chunk, range(chunk_count)):
                    output.write(block)
            assert temporary.stat().st_size == total_size

        validate_gzip(temporary)
        temporary.replace(destination)


if __name__ == "__main__":
    download()
