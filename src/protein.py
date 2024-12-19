import lzma
import multiprocessing
from typing import List, Dict, Tuple, Any, Union
import pickle
from .annotation import compress_annotation, decompress_annotation

AMINO_ACID_GROUPS = {
    'A': 0b00000000,
    'I': 0b00000001,
    'L': 0b00000010,
    'M': 0b00000011,
    'F': 0b00000100,
    'W': 0b00000101,
    'V': 0b00000110,
    'S': 0b00001000,
    'T': 0b00001001,
    'N': 0b00001010,
    'Q': 0b00001011,
    'K': 0b00010000,
    'R': 0b00010001,
    'H': 0b00010010,
    'D': 0b00011000,
    'E': 0b00011001,
    'Y': 0b00100000,
    'G': 0b00100001,
    'C': 0b00100010,
    'P': 0b00100011,
    '*': 0b00101000,
    'B': 0b00110000,  
    'J': 0b00110001,
    'O': 0b00110010,
    'U': 0b00110011,  
    'X': 0b00110100,  
    'Z': 0b00110101,  
    '-': 0b00111000,
}

def read_fasta(file_path: str) -> List[Tuple[str, str]]:
    sequences = []
    current_annotation = ""
    current_sequence = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences.append((current_annotation, ''.join(current_sequence)))
                    current_sequence = []
                current_annotation = line
            else:
                current_sequence.append(line)

        if current_sequence:
            sequences.append((current_annotation, ''.join(current_sequence)))

    return sequences

def protein_encode_sequence(sequence: str) -> bytes:
    sequence = ''.join(char.upper() if char.upper() in AMINO_ACID_GROUPS else 'X' for char in sequence)
    result_bytes = bytearray()
    for char in sequence:
        result_bytes.append(AMINO_ACID_GROUPS[char])
    
    return bytes(result_bytes)

def protein_decode_sequence(compressed: bytes) -> str:
    reverse_mapping = {code: aa for aa, code in AMINO_ACID_GROUPS.items()}
    return ''.join(reverse_mapping.get(byte, 'A') for byte in compressed)

def protein_compress_single(args: Tuple[str, str]) -> Tuple[str, bytes]:
    annotation, sequence = args
    compressed_sequence = protein_encode_sequence(sequence)
    return annotation, compressed_sequence

def protein_decompress_single(compressed_data: Tuple[str, bytes]) -> Tuple[str, str]:
    annotation, compressed_sequence = compressed_data
    decompressed_sequence = protein_decode_sequence(compressed_sequence)
    return annotation, decompressed_sequence

def protein_parallel_compress(sequences: List[Tuple[str, str]], processes: int = None) -> List[Tuple[str, bytes]]:
    with multiprocessing.Pool(processes=processes) as pool:
        return pool.map(protein_compress_single, sequences)

def protein_parallel_decompress(compressed_data: List[Tuple[str, bytes]], processes: int = None) -> List[Tuple[str, str]]:
    with multiprocessing.Pool(processes=processes) as pool:
        return pool.map(protein_decompress_single, compressed_data)

def protein_write_compressed(data: List[Tuple[str, bytes]], output_file: str, preset: int = 3):
    annotations = [item[0] for item in data]
    annotation_patterns, variable_parts = compress_annotation(annotations)

    compressed_data = {
        'format': 'PROTEIN_FASTA',
        'annotation_patterns': annotation_patterns,
        'variable_parts': variable_parts,
        'sequences': [seq for _, seq in data]
    }

    with open(output_file, 'wb') as f_out:
        compressed = lzma.compress(pickle.dumps(compressed_data), preset=preset)
        f_out.write(compressed)

def protein_read_compressed(file_path: str) -> List[Tuple[str, bytes]]:
    with open(file_path, 'rb') as f_in:
        compressed_data = pickle.loads(lzma.decompress(f_in.read()))
    
    if compressed_data.get('format') != 'PROTEIN_FASTA':
        raise ValueError("Unsupported file format")

    annotations = decompress_annotation(
        compressed_data['annotation_patterns'],
        compressed_data['variable_parts']
    )

    return list(zip(annotations, compressed_data['sequences']))
