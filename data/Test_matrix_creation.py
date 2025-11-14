"""Hand-crafted alignments for exercising :func:`P3.ComputeDistMatrix`."""

import sys
from pathlib import Path
from typing import Dict, Tuple


# Ensure that the repository root is on sys.path so that ``import P3`` works
# even when the script is executed from inside ``data`` or any other folder.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


AlignmentMap = Dict[Tuple[int, int], Tuple[str, str]]

ALIGNMENTS: AlignmentMap = {
    (1, 2): (
        "ACGTCGTAACAA",
        "ACGTCGTTACGT",
    ),
    (1, 3): (
        "ACGTACGT--ACGT",
        "ACGTTCGTATGCGT",
    ),
    (1, 4): (
        "ACGTACGTACACGTACGT--ACGTACGTACGTAAACGTTCGTATGCGT",
        "ACGTACGTAAACGTTCGTATGCGTACGTACGTACACGTACGT--ACGT",
    ),
    (2, 3): (
        "ACGTACGT--ACGT",
        "ACGTTCGTATGCGT",
    ),
    (2, 4): (
        "ACGTACGT--ACGT",
        "ACGTTCGTATGCGT",
    ),
    (3, 4): (
        "ACGTACGTACACGTACGT--ACGTACGTACACGTACGTGTAAACGTTCGTATGCGT",
        "ACGTACGTAAACGTTCGTATGCACGTACGTGTACGTACGTACACGTACGT--ACGT",
    ),
}


if __name__ == "__main__":
    from P3 import ComputeDistMatrix

    matrix = ComputeDistMatrix(ALIGNMENTS)
    for row in matrix:
        print(", ".join(f"{value:.6f}" for value in row))
