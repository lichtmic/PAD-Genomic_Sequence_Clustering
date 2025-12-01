def read_file(filename):
    """
    Read the contents of a sequence file.
    
    Args:
        filename (str): Path to the sequence file
        
    Returns:
        str: Contents of the file as a string
        
    Raises:
        Exception: If file is not found, raises "malformed input"
    """
    try:
        with open(filename, 'r') as file:
            return file.read()
    except FileNotFoundError:
        raise Exception("malformed input")

def extract_lines_starting_with_greater_than(content):
    """
    Extract and parse lines that start with '>' character.
    
    Lines starting with '>' contain sequence data (label and nucleotides).
    Empty lines are skipped. Non-empty lines not starting with '>' are invalid.
    
    Args:
        content (str): File content as a single string
        
    Returns:
        list[str]: List of parsed data lines (without the '>' prefix)
        
    Raises:
        Exception: If non-empty line doesn't start with '>', raises "malformed input"
    """
    lines = content.split('\n')
    result = []
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            line = line[1:].strip()
            if not line:  # Check if empty after removing '>'
                raise Exception("malformed input")
            result.append(line)
        elif line:  # Non-empty line that doesn't start with '>'
            raise Exception("malformed input")
    return result

def extract_organism(data):
    """
    Extract organism labels and raw sequences from parsed lines.
    
    Each line should contain a label followed by a sequence.
    The sequence may contain spaces which will be preserved for later cleaning.
    
    Args:
        data (list[str]): List of parsed data lines (label + sequence)
        
    Returns:
        tuple[list[str], list[str]]: A tuple containing:
            - List of organism labels (capitalized)
            - List of raw sequences (with spaces)
            
    Raises:
        Exception: If line doesn't contain both label and sequence, raises "malformed input"
    """
    org_names = []
    seq_raw = []
    for line in data:
        parts = line.split()
        if len(parts) < 2:  # Must have label and sequence
            raise Exception("malformed input")
        
        name = parts[0]
        seq_raw.append(' '.join(parts[1:]))
        org_names.append(name.lower().capitalize())  # normalize case
    return org_names, seq_raw

def clean_and_check_seqs(raw):
    """
    Clean sequences and validate that they contain only valid nucleotides (both upper and lower case).
    
    Removes all whitespace and converts to uppercase. Validates that all
    characters are valid nucleotides (A, T, C, G).
    
    Args:
        raw (list[str]): List of raw sequences potentially containing spaces
        
    Returns:
        list[str]: List of cleaned, uppercase sequences
        
    Raises:
        Exception: If sequence is empty or contains invalid characters, 
                   raises "malformed input"
    """
    clean_seqs = []
    acceptableChars = set('ATCGatcg')  # Accept both upper and lower case Sequences
    
    for seq in raw:
        clean_seq = seq.replace(' ', '').upper()
        
        if not clean_seq:  # Empty sequence
            raise Exception("malformed input")
        
        # Check if all characters are valid
        if not all(chars in acceptableChars for chars in clean_seq):
            raise Exception("malformed input")
        
        clean_seqs.append(clean_seq)
    
    return clean_seqs

def generate_output_list(organisms, sequences):
    """
    Combine organism labels and sequences into list of tuples.
    
    Args:
        organisms (list[str]): List of organism labels
        sequences (list[str]): List of cleaned sequences
        
    Returns:
        list[tuple[str, str]]: List of (label, sequence) pairs
        
    Raises:
        Exception: If lists have different lengths, raises "malformed input"
    """
    if len(organisms) != len(sequences):
        raise Exception("malformed input")
    
    return [(org, seq) for org, seq in zip(organisms, sequences)]

def ParseSeqFile(filename):
    """
    Parse a genomic sequence file and extract organism-sequence pairs.
    
    Args:
        filename (str): Path to the sequence file
        
    Returns:
        list[tuple[str, str]]: List of (label, sequence) pairs where
                               labels are capitalized and sequences are
                               uppercase without spaces
    """
    content = read_file(filename)
    lines = extract_lines_starting_with_greater_than(content)
    org_names, seq_raw = extract_organism(lines)
    seqs = clean_and_check_seqs(seq_raw)
    list_of_pairs = generate_output_list(org_names, seqs)
    return list_of_pairs
