"""Check that an sdist carries the paper sources without generated results.

Run from the source checkout after building the archive. This compares the
archive with the checkout; --source-root can select another source tree.
"""

import argparse
from pathlib import Path, PurePosixPath
import tarfile


def check_archive(archive: Path, source_root: Path) -> list[str]:
    required = {
        "LICENSE",
        "README.md",
        "CITATION.cff",
        "CMakeLists.txt",
        "cmake/rebound_bridgeConfig.cmake.in",
        "tests/paper/common/validation_timer.h",
        "tests/paper/CMakeLists.txt",
        "tests/paper/README.md",
        "tests/paper/run_all.py",
        "tests/paper/requirements.txt",
        "tests/paper/mnras.mplstyle",
        "tests/paper/common/galilean_horizons_state.h",
        "tests/paper/checks/check_paper_sdist.py",
        "rebound/src/rebound.h",
    }
    for directory, extensions in (
        ("src", {".c", ".h"}),
        ("rebound/src", {".c", ".h"}),
        ("tests/paper", {".c", ".h", ".py", ".md", ".txt", ".mplstyle"}),
    ):
        for path in (source_root / directory).rglob("*"):
            relative = path.relative_to(source_root)
            if path.is_file() and path.suffix in extensions:
                if not {"results", "__pycache__"}.intersection(relative.parts):
                    required.add(relative.as_posix())

    errors = []
    with tarfile.open(archive, "r:gz") as package:
        files = {}
        for member in package.getmembers():
            if not member.isfile():
                continue
            path = PurePosixPath(member.name)
            relative = PurePosixPath(*path.parts[1:])
            files[relative.as_posix()] = member
            if relative.parts[:1] == ("tests",) and relative.parts[:2] != ("tests", "paper"):
                errors.append(f"Non-paper test source: {relative}")
            if {"results", "__pycache__", ".pytest_cache", "build"}.intersection(relative.parts):
                errors.append(f"Generated content: {relative}")
            if relative.suffix.lower() in {".pyc", ".pyo", ".o", ".obj", ".exe", ".dll", ".pyd", ".so"}:
                errors.append(f"Compiled output: {relative}")
        for name in sorted(required):
            if name not in files:
                errors.append(f"Missing source: {name}")
            elif not (source_root / name).is_file():
                errors.append(f"Missing comparison source: {name}")
            elif package.extractfile(files[name]).read() != (source_root / name).read_bytes():
                errors.append(f"Content differs from source: {name}")
    if not errors:
        print(f"PASS: {len(required)} paper/build source files match; no generated results or binaries.")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[3])
    arguments = parser.parse_args()
    errors = check_archive(arguments.archive, arguments.source_root.resolve())
    for error in errors:
        print(error)
    return bool(errors)


if __name__ == "__main__":
    raise SystemExit(main())
