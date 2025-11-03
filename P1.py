def read_file(filename):
    try:
        with open(filename, 'r') as file:
            return file.read()
    except FileNotFoundError:
        raise Exception("malformed input")

def extract_lines_starting_with_greater_than(content):
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
    clean_seqs = []
    acceptableChars = set('ATCGatcg')  # Accept both cases
    
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
    if len(organisms) != len(sequences):
        raise Exception("malformed input")
    
    return [(org, seq) for org, seq in zip(organisms, sequences)]

def ParseSeqFile(filename):
    content = read_file(filename)
    lines = extract_lines_starting_with_greater_than(content)
    org_names, seq_raw = extract_organism(lines)
    seqs = clean_and_check_seqs(seq_raw)
    list_of_pairs = generate_output_list(org_names, seqs)
    return list_of_pairs
