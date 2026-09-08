"""Download the public source files for the South Korea feasibility analysis."""

from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path


RAW_DIR = Path(__file__).resolve().parent / "raw"
KOSIS_ENDPOINT = "https://kosis.kr/file_mass/file_down.jsp"
KDCA_ENDPOINT = "https://chs.kdca.go.kr/cscdhpfile/abs/fileCmmn/fileDown.do"
BOUNDARY_URL = "http://www.gisdeveloper.co.kr/download/admin_shp/SIG_201905.zip"

KOSIS_FILES = {
    "municipal_mortality": {
        2015: ("4895", "101_DT_1B34E13_Y_2015"),
        2016: ("5887", "101_DT_1B34E13_Y_2016"),
        2017: ("6364", "101_DT_1B34E13_Y_2017"),
        2018: ("6982", "101_DT_1B34E13_Y_2018"),
        2019: ("7195", "101_DT_1B34E13_Y_2019"),
    },
    "midyear_population": {
        2015: ("4828", "101_DT_1B040M5_Y_2015"),
        2016: ("5398", "101_DT_1B040M5_Y_2016"),
        2017: ("5992", "101_DT_1B040M5_Y_2017"),
        2018: ("6769", "101_DT_1B040M5_Y_2018"),
        2019: ("7032", "101_DT_1B040M5_Y_2019"),
    },
    "national_mortality": {
        2015: ("4896", "101_DT_1B34E11_Y_2015"),
        2016: ("5886", "101_DT_1B34E11_Y_2016"),
        2017: ("6371", "101_DT_1B34E11_Y_2017"),
        2018: ("6981", "101_DT_1B34E11_Y_2018"),
        2019: ("7194", "101_DT_1B34E11_Y_2019"),
    },
}

KOSIS_OPTIONS = {
    "municipal_mortality": {"list_id": "", "vw_cd": "MT_ZTITLE", "language": "eng"},
    "midyear_population": {"list_id": "A_7", "vw_cd": "MT_ZTITLE", "language": "kor"},
    "national_mortality": {"list_id": "", "vw_cd": "", "language": "eng"},
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _is_complete_archive(path: Path) -> bool:
    return path.exists() and zipfile.is_zipfile(path)


def _stream_response(response, destination: Path) -> None:
    temporary = destination.with_suffix(destination.suffix + ".part")
    with temporary.open("wb") as output:
        while chunk := response.read(1024 * 1024):
            output.write(chunk)
    temporary.replace(destination)


def _download_get(url: str, destination: Path, force: bool) -> None:
    if not force and _is_complete_archive(destination):
        print(f"Using existing archive: {destination.name}")
        return
    print(f"Downloading: {destination.name}")
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        _stream_response(response, destination)
    if not _is_complete_archive(destination):
        raise RuntimeError(f"Downloaded file is not a complete archive: {destination}")


def _download_post(url: str, fields: dict[str, str], destination: Path, force: bool) -> None:
    if not force and _is_complete_archive(destination):
        print(f"Using existing archive: {destination.name}")
        return
    print(f"Downloading: {destination.name}")
    payload = urllib.parse.urlencode(fields).encode("ascii")
    request = urllib.request.Request(
        url,
        data=payload,
        headers={"User-Agent": "Mozilla/5.0"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        _stream_response(response, destination)
    if not _is_complete_archive(destination):
        raise RuntimeError(f"Downloaded file is not a complete archive: {destination}")


def download_kosis(force: bool) -> None:
    for group, files in KOSIS_FILES.items():
        options = KOSIS_OPTIONS[group]
        for year, (file_number, base_name) in files.items():
            destination = RAW_DIR / f"kosis_{group}_{year}.zip"
            language = options["language"]
            fields = {
                "tbl_id": base_name.rsplit("_Y_", 1)[0].replace("101_", "", 1),
                "org_id": "101",
                "filename": f"{base_name}_ENG" if language == "eng" else base_name,
                "file_no": file_number,
                "file_type": "ONE",
                "vw_cd": options["vw_cd"],
                "list_id": options["list_id"],
                "usrId": "null",
                "usrName": "null",
                "down_cnt": "1",
                "use_no": "0",
                "page": "kosis",
            }
            if language == "eng":
                fields["scr_at"] = "eng"
            _download_post(KOSIS_ENDPOINT, fields, destination, force)


def download_kdca(force: bool) -> None:
    destination = RAW_DIR / "kdca_health_database_v1_7.xlsx"
    _download_post(KDCA_ENDPOINT, {"SEQ": "19d700db9c62"}, destination, force)


def download_boundaries(force: bool) -> None:
    _download_get(BOUNDARY_URL, RAW_DIR / "SIG_201905.zip", force)


def write_checksums() -> None:
    rows = ["filename,sha256,size_bytes"]
    for path in sorted(RAW_DIR.iterdir()):
        if path.is_file() and path.suffix.lower() in {".zip", ".xlsx"}:
            rows.append(f"{path.name},{_sha256(path)},{path.stat().st_size}")
    (RAW_DIR / "checksums.csv").write_text("\n".join(rows) + "\n", encoding="ascii")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="Replace complete existing downloads.")
    parser.add_argument("--skip-kosis", action="store_true")
    parser.add_argument("--skip-kdca", action="store_true")
    parser.add_argument("--skip-boundaries", action="store_true")
    args = parser.parse_args()

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    if not args.skip_kosis:
        download_kosis(args.force)
    if not args.skip_kdca:
        download_kdca(args.force)
    if not args.skip_boundaries:
        download_boundaries(args.force)
    write_checksums()
    return 0


if __name__ == "__main__":
    sys.exit(main())
