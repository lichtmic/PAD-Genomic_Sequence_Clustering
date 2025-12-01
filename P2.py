import numpy as np

def input_check(list_of_pairs):
    """
    Validate input format and clean sequences.
    
    Checks that input is a list of tuples where each tuple contains
    two strings (label, sequence). Cleans sequences by removing spaces
    and converting to uppercase.
    
    Args:
        list_of_pairs (list[tuple[str, str]]): List of (label, sequence) pairs
        
    Returns:
        list[tuple[str, str]]: Cleaned list of (label, sequence) pairs
        
    Raises:
        Exception: If input format is invalid, raises "malformed input"
    """
    if not isinstance(list_of_pairs, list):
        raise Exception("malformed input")
    if not all(isinstance(pair, tuple) and len(pair) == 2 for pair in list_of_pairs):
        raise Exception("malformed input")
    if not all(isinstance(pair[0], str) and isinstance(pair[1], str) for pair in list_of_pairs):
        raise Exception("malformed input")
    else:
        cleaned = [(label, seq.replace(" ", "").upper()) for label, seq in list_of_pairs]
        return cleaned
    
def generate_all_pairs(list_of_pairs):
    """
    Generate all unique index pairs (i, j) for given list of pairs.
    
    Args:
        list_of_pairs (list[tuple[str, str]]): List of (label, sequence) pairs
        
    Returns:
        list[tuple[int, int]]: List of unique index pairs (i, j) with i < j
    """
    n = len(list_of_pairs)
    pairs = []
    # Generate all unique index pairs (i, j) with i < j
    for i in range(n):
        for j in range(i+1, n):
            pairs.append((i, j))
    return pairs

def fill_matrix(seq1, seq2, match=5, mismatch=-2, gap=-6):
    """
    Fill dynamic programming matrix for sequence alignment.
    
    Implements a dynamic programming algorithm to find
    the optimal global alignment between two sequences.
    
    Args:
        seq1 (str): First nucleotide sequence
        seq2 (str): Second nucleotide sequence
        match (int, optional): Score for matching nucleotides. Default: 5
        mismatch (int, optional): Score for mismatched nucleotides. Default: -2
        gap (int, optional): Penalty for gaps. Default: -6
        
    Returns:
        tuple[np.ndarray, np.ndarray]: A tuple containing:
            - score: (m+1 × n+1) matrix of alignment scores
            - traceback: (m+1 × n+1) matrix of traceback directions
                         (0=diagonal, 1=up, 2=left)
    """
    m, n = len(seq1), len(seq2)
    score = np.zeros((m + 1, n + 1))
    traceback = np.zeros((m + 1, n + 1), dtype=int)
    
    # Initialize first column
    for i in range(1, m + 1):
        score[i, 0] = i * gap
        traceback[i, 0] = 1  # up
    
    # Initialize first row
    for j in range(1, n + 1):
        score[0, j] = j * gap
        traceback[0, j] = 2  # left
    
    # Fill the matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                diag = score[i-1, j-1] + match
            else:
                diag = score[i-1, j-1] + mismatch
            
            up = score[i-1, j] + gap
            left = score[i, j-1] + gap
            
            max_score = max(diag, up, left)
            score[i, j] = max_score
            
            # Store traceback direction
            if max_score == diag:
                traceback[i, j] = 0  # diagonal
            elif max_score == up:
                traceback[i, j] = 1  # up
            else:
                traceback[i, j] = 2  # left
    
    return score, traceback

def traceback(seq1, seq2, traceback): 
    """
    Reconstruct optimal alignment by tracing back through dynamic programming matrix.
    
    Starts from the bottom-right corner of the traceback matrix and
    follows the path of traceback directions back to the origin, building
    the aligned sequences with gaps as needed.
    
    Args:
        seq1 (str): First nucleotide sequence
        seq2 (str): Second nucleotide sequence
        traceback (np.ndarray): Matrix of traceback directions where
                                0=diagonal, 1=up, 2=left
        
    Returns:
        tuple[str, str]: Aligned versions of seq1 and seq2 with gaps ('-')
                        inserted to show optimal alignment. Both aligned
                        sequences have the same length.

    """
    aligned1 = []
    aligned2 = []
    i, j = len(seq1), len(seq2)
    
    while i > 0 or j > 0:
        if i == 0:  # Top edge (only horizontal moves left possible)
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
        elif j == 0:  # Left edge (only vertical moves up possible)
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        elif traceback[i, j] == 0:  # Diagonal
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif traceback[i, j] == 1:  # Up
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        else:  # traceback[i, j] == 2, Left
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
    
    aligned1.reverse()
    aligned2.reverse()
    
    return ''.join(aligned1), ''.join(aligned2)

def align_sequences_with_DP(seq1, seq2):
    """
    Align two sequences using dynamic programming.
    
    Combines the matrix filling and traceback steps to produce
    the optimal global alignment of two nucleotide sequences.
    
    Args:
        seq1 (str): First nucleotide sequence
        seq2 (str): Second nucleotide sequence
        
    Returns:
        tuple[str, str]: Aligned versions of seq1 and seq2 with gaps
        
    """

    score, traceback = fill_matrix(seq1, seq2)
    return traceback(seq1, seq2, traceback)

def AlignByDP(list_of_pairs):   # Main function to align sequences by dynamic programming
    """
    Align all pairwise combinations of sequences using dynamic programming.
    
    Takes a list of (label, sequence) pairs and performs optimal global
    alignment on each unique pair using a dynamic programming algorithm
    with scoring parameters: match=+5, mismatch=-2, gap=-6.
    
    Args:
        list_of_pairs (list[tuple[str, str]]): List of (label, sequence) pairs
        
    Returns:
        dict[tuple[int, int], tuple[str, str]]: Dictionary where keys are
            index pairs (i, j) and values are tuples of aligned sequences
            with gaps inserted. Both sequences in each pair have equal length.
            
    Raises:
        Exception: If input format is invalid, raises "malformed input"
    """
    sequences = input_check(list_of_pairs)
    results = {}
    for i, j in generate_all_pairs(sequences):
        seq1, seq2 = sequences[i][1], sequences[j][1]
        aligned1, aligned2 = align_sequences_with_DP(seq1, seq2)
        results[(i, j)] = (aligned1, aligned2)
    return results
