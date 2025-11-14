from typing import Dict, List, Tuple
from itertools import combinations
from math import log


VALID_CHARS = {"A", "C", "G", "T", "a", "c", "g", "t", "-"}


def ensure(condition):
    """Raise Exeption "malformed input" when condition fails."""

    if not condition:
        raise Exception("malformed input")


def validate_alignment_map(
    alignment_map: Dict[Tuple[int, int], Tuple[str, str]],
) -> Tuple[List[int], Dict[Tuple[int, int], Tuple[str, str]]]:
    """Normalise and validate the input alignment dictionary.

    Ensures keys are ordered integer pairs, values are equal-length strings
    containing only valid characters, every pair occurs exactly once, and
    returns the sorted list of sequence IDs alongside the cleaned mapping.
    """

    # input is a dict that is not empty
    ensure(isinstance(alignment_map, dict) and alignment_map)

    # creates dict for ouput
    normalised: Dict[Tuple[int, int], Tuple[str, str]] = {}
    # creates a set to track unique sequences
    identifiers: set[int] = set()

    for dict_key, dict_value in alignment_map.items():
        # --- Key Validation ---
        ensure(isinstance(dict_key, tuple) and len(dict_key) == 2)
        i, j = dict_key
        ensure(all(isinstance(idx, int) for idx in dict_key))
        # Reject self-alignment i.e. (1,1) or (2,2)
        ensure(i != j)
        # Normalize key to (min_index, max_index) to handle symmetry.
        ordered_key = (i, j) if i < j else (j, i)

        # --- Value Validation ---
        ensure(isinstance(dict_value, tuple) and len(dict_value) == 2)
        seq_a, seq_b = dict_value
        ensure(isinstance(seq_a, str) and isinstance(seq_b, str))
        ensure(len(seq_a) == len(seq_b) and len(seq_a) > 0)

        ensure(set(seq_a).issubset(VALID_CHARS) and set(seq_b).issubset(VALID_CHARS))

        """
        should the loop stop if we have duplicates or should duplicates 
        be ignored so that we get all the data stored inside the variable ??
        """

        # Consistency Check (Avoiding Duplicates)
        ensure(ordered_key not in normalised)
        # Storing the Clean Data
        normalised[ordered_key] = (seq_a, seq_b)
        identifiers.update(ordered_key)

    ids_sorted = sorted(identifiers)
    # Must have at least two unique sequences for pairwise comparison.
    ensure(len(ids_sorted) >= 2)

    # Calculate the full set of pairs expected from all found identifiers.
    expected_pairs = set(combinations(ids_sorted, 2))
    # Final check: ensure the input map contains the exact, complete set of expected pairs.
    ensure(expected_pairs == set(normalised))

    return ids_sorted, normalised


def compute_dist(aligned_a, aligned_b):
    """Return the Jukes-Cantor distance between two aligned sequences."""

    comparable = 0
    mismatches = 0

    for char_a, char_b in zip(aligned_a, aligned_b):
        if "-" in (char_a, char_b):
            continue
        comparable += 1
        if char_a != char_b:
            mismatches += 1

    if comparable == 0:
        return 0.0

    p_distance = mismatches / comparable
    if p_distance >= 0.75:
        return 30.0

    correction = 1.0 - (4.0 / 3.0) * p_distance
    ensure(correction > 0.0)

    return -0.75 * log(correction)


def ComputeDistMatrix(alignment_map):
    """Create the symmetric distance matrix for all sequence pairs.

    Parameters:
    alignment_map:
            Mapping from pairs of sequence IDs to their aligned nucleotide strings.

    Returns:
    list[list[float]]
            Symmetric matrix "D" with "D[i][j]" storing the evolutionary
            distance between sequence "ids[i]" and "ids[j]".
    """

    ids_sorted, normalised = validate_alignment_map(alignment_map)

    size = len(ids_sorted)
    index_map = {seq_id: pos for pos, seq_id in enumerate(ids_sorted)}
    matrix: List[List[float]] = [[0.0] * size for _ in range(size)]

    for (i, j), alignment in normalised.items():
        dist = compute_dist(*alignment)
        idx_i, idx_j = index_map[i], index_map[j]
        matrix[idx_i][idx_j] = dist
        matrix[idx_j][idx_i] = dist

    return matrix
