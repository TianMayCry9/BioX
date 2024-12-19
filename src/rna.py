import lzma
import multiprocessing
from typing import List, Dict, Tuple, Any, Union
import pickle
from annotation import compress_annotation, decompress_annotation
from quality import  encode_quality_smart, decode_quality_smart


def detect_file_format(file_path: str) -> str:
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()
        if first_line.startswith('>'):
            return 'fasta'
        elif first_line.startswith('@'):
            return 'fastq'
        else:
            raise ValueError("Error")


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

def read_fastq(file_path: str) -> List[Tuple[str, str, str, str]]: 
    sequences = []
    with open(file_path, 'r') as file:
        current_record = []
        for line in file:
            line = line.strip()
            if not line:  
                continue
            if line.startswith('@') and len(current_record) == 0:
                current_record = [line]
            elif line.startswith('@') and len(current_record) == 4:
                sequences.append((current_record[0], current_record[1], current_record[2], current_record[3]))  
                current_record = [line]
            elif len(current_record) < 4:
                current_record.append(line)

        if len(current_record) == 4:
            sequences.append((current_record[0], current_record[1], current_record[2], current_record[3]))  

    return sequences

def encode_length(length: int) -> List[int]:
    bits = []
    if length < 255:
        bits.append(0)
        for i in range(7, -1, -1):
            bits.append((length >> i) & 1)
        return bits
    elif length < 65535:
        bits.extend([1, 0])
        for i in range(15, -1, -1):
            bits.append((length >> i) & 1)
        return bits
    else:
        bits.extend([1, 1])
        for i in range(31, -1, -1):
            bits.append((length >> i) & 1)
        return bits

def decode_length(bits: List[int], start: int) -> Tuple[int, int]:
    if bits[start] == 0:
        length = 0
        for i in range(8):
            length = (length << 1) | bits[start + 1 + i]
        return length, start + 9
    elif bits[start] == 1 and bits[start + 1] == 0:
        length = 0
        for i in range(16):
            length = (length << 1) | bits[start + 2 + i]
        return length, start + 18
    else:
        length = 0
        for i in range(32):
            length = (length << 1) | bits[start + 2 + i]
        return length, start + 34


def dna_encode_sequence(sequence: str) -> bytes:
    sequence = sequence.upper()
    segments = []
    current_type = None 
    current_start = 0
    current_segment = []
    
    for i, char in enumerate(sequence):
        is_special = char not in 'AUCG'
        if current_type is not None and is_special != current_type:
            segments.append((current_type, current_start, ''.join(current_segment)))
            current_start = i
            current_segment = []
        
        current_type = is_special
        current_segment.append(char)
    
    if current_segment:
        segments.append((current_type, current_start, ''.join(current_segment)))

    merged_segments = []
    current_merged = None
    
    for segment in segments:
        is_special, start, data = segment
        
        if current_merged is None:
            current_merged = list(segment)
        elif current_merged[0] == is_special and len(current_merged[2]) < 10:
            current_merged[2] += data
        else:
            merged_segments.append(tuple(current_merged))
            current_merged = list(segment)

    if current_merged:
        merged_segments.append(tuple(current_merged))

    bits = []
    
    total_segments = len(merged_segments)
    for i in range(31, -1, -1):
        bits.append((total_segments >> i) & 1)

    for is_special, _, data in merged_segments:
        bits.append(1 if is_special else 0)

        bits.extend(encode_length(len(data)))
        
        if is_special:
            for char in data:
                ascii_val = ord(char)
                for i in range(8):
                    bits.append((ascii_val >> (7-i)) & 1)
        else:
            for char in data:
                if char == 'A':
                    bits.extend([0, 0])
                elif char == 'U':
                    bits.extend([1, 1])
                elif char == 'C':
                    bits.extend([0, 1])
                elif char == 'G':
                    bits.extend([1, 0])

    while len(bits) % 8 != 0:
        bits.append(0)
    
    result_bytes = bytearray()
    for i in range(0, len(bits), 8):
        byte = 0
        for j in range(8):
            if i + j < len(bits):
                byte = (byte << 1) | bits[i + j]
        result_bytes.append(byte)
    
    return bytes(result_bytes)

def dna_decode_sequence(compressed: bytes) -> str:
    decoded = []
    
    bits = []
    for byte in compressed:
        for i in range(7, -1, -1):
            bits.append((byte >> i) & 1)

    total_segments = 0
    for i in range(32):
        total_segments = (total_segments << 1) | bits[i]
    i = 32

    for _ in range(total_segments):
        if i >= len(bits):
            break

        is_special = bits[i]
        i += 1

        length, i = decode_length(bits, i)
        
        if is_special:
            for _ in range(length):
                if i + 8 > len(bits):
                    break
                ascii_val = 0
                for j in range(8):
                    ascii_val = (ascii_val << 1) | bits[i + j]
                decoded.append(chr(ascii_val))
                i += 8
        else:
            for _ in range(length):
                if i + 2 > len(bits):
                    break
                two_bits = bits[i:i+2]
                if two_bits == [0, 0]:
                    decoded.append('A')
                elif two_bits == [1, 1]:
                    decoded.append('U')
                elif two_bits == [0, 1]:
                    decoded.append('C')
                elif two_bits == [1, 0]:
                    decoded.append('G')
                i += 2
    
    return ''.join(decoded)

DNA_CODES = {
    'A': 0b00000000,
    'U': 0b00000001,
    'C': 0b00000010,
    'G': 0b00000011,
    'N': 0b00010011,
    'R': 0b00010100,
    'Y': 0b00010101,
    'M': 0b00010110,
    'K': 0b00010111,
    'S': 0b00011000,
    'W': 0b00011001,
    'H': 0b00011010,
    'B': 0b00011011,
    'V': 0b00011100,
    'D': 0b00011101,
}

CHUNK_SIZE = 8 * 1024 * 1024 * 1024

def dna_encode_sequence_plant(sequence: str) -> bytes:
    sequence = sequence.upper() 
    result_bytes = bytearray()
    for char in sequence:
        if char in DNA_CODES:
            result_bytes.append(DNA_CODES[char])
        else:
            result_bytes.append(ord(char))
    return bytes(result_bytes)

def dna_decode_sequence_plant(compressed: bytes) -> str:
    reverse_mapping = {code: base for base, code in DNA_CODES.items()}
    return ''.join(reverse_mapping.get(byte, chr(byte)) for byte in compressed)

def process_sequence_chunks(sequence: str, encode_func) -> bytes:
    if len(sequence) <= CHUNK_SIZE:
        return encode_func(sequence)
    
    result = bytearray()
    for i in range(0, len(sequence), CHUNK_SIZE):
        chunk = sequence[i:i + CHUNK_SIZE]
        result.extend(encode_func(chunk))
    
    return bytes(result)

def dna_compress_single_fasta_sequence(args: Tuple[str, str, bool]) -> Tuple[str, str]:
    annotation, sequence, is_plant = args
    if is_plant:
        compressed_sequence = process_sequence_chunks(sequence, dna_encode_sequence_plant)
    else:
        compressed_sequence = process_sequence_chunks(sequence, dna_encode_sequence)
    return annotation, compressed_sequence

def dna_compress_single_sequence(args: Tuple[str, str, bool, str, str]) -> Tuple[str, str, bool, str, Tuple[bytes, bytes, int]]:
    annotation, sequence, is_plant, plus_line, quality = args
    if is_plant:
        compressed_sequence = process_sequence_chunks(sequence, dna_encode_sequence_plant)
    else:
        compressed_sequence = process_sequence_chunks(sequence, dna_encode_sequence)
    compressed_quality, metadata, encoding_type = encode_quality_smart(quality)
    return annotation, compressed_sequence, plus_line, (compressed_quality, metadata, encoding_type)

def dna_parallel_compress_fasta(sequences: List[Tuple[str, str]], 
                              processes: int = None,
                              is_plant: bool = False) -> List[Tuple[str, str]]:
    processes = processes or multiprocessing.cpu_count()
    sequences_with_flag = [(ann, seq, is_plant) for ann, seq in sequences]
    with multiprocessing.Pool(processes=processes) as pool:
        return pool.map(dna_compress_single_fasta_sequence, sequences_with_flag)


def dna_parallel_compress(sequences: List[Tuple[str, str, str, str]], 
                        processes: int = None, 
                        is_plant: bool = False) -> List[Tuple[str, str, str, Tuple[bytes, bytes, int]]]:
    processes = processes or multiprocessing.cpu_count()
    sequences_with_flag = [(ann, seq, is_plant, plus, qual) for ann, seq, plus, qual in sequences]
    with multiprocessing.Pool(processes=processes) as pool:
        return pool.map(dna_compress_single_sequence, sequences_with_flag)

def dna_write_compressed_fastq(data: List[Tuple[str, str, str, Tuple[bytes, bytes, int]]], output_file: str, compress_plus: bool = False, preset: int = 3, is_plant: bool = False):
    annotations = [item[0] for item in data]
    annotation_patterns, variable_parts = compress_annotation(annotations)

    plus_patterns = None
    plus_variables = None
    if compress_plus:
        plus_lines = [item[2] for item in data]
        plus_patterns, plus_variables = compress_annotation(plus_lines)

    sequences_data = [(seq, qual_data) for _, seq, _, qual_data in data]
 
    compressed_data = {
        'format': 'FASTQ',
        'is_plant': is_plant,  
        'annotation_patterns': annotation_patterns,
        'variable_parts': variable_parts,
        'plus_compressed': compress_plus,
        'plus_patterns': plus_patterns,
        'plus_variables': plus_variables,
        'sequences': sequences_data
    }

    with open(output_file, 'wb') as f_out:
        compressed = lzma.compress(pickle.dumps(compressed_data), preset=preset)
        f_out.write(compressed)


def dna_write_compressed_fasta(data: List[Tuple[str, str]], output_file: str, preset: int = 3, is_plant: bool = False):
    annotations = [item[0] for item in data]

    annotation_patterns, variable_parts = compress_annotation(annotations)

    compressed_data = {
        'format': 'FASTA',
        'is_plant': is_plant,  
        'annotation_patterns': annotation_patterns,
        'variable_parts': variable_parts,
        'sequences': [seq for _, seq in data]
    }

    with open(output_file, 'wb') as f_out:
        compressed = lzma.compress(pickle.dumps(compressed_data), preset=preset)
        f_out.write(compressed)

def dna_read_compressed_file(file_path: str) -> Tuple[str, Union[List[Tuple[str, str]], List[Tuple[str, str, str, Tuple[bytes, bytes, int]]]]]:
    with open(file_path, 'rb') as f_in:
        compressed_data = pickle.loads(lzma.decompress(f_in.read()))
    
    file_format = compressed_data.get('format')
    is_plant = compressed_data.get('is_plant', False)
    
    if not file_format:
        raise ValueError("Error")

    if file_format == 'FASTA':
        annotations = decompress_annotation(
            compressed_data['annotation_patterns'],
            compressed_data['variable_parts']
        )
        return 'FASTA', list(zip(annotations, compressed_data['sequences'], [is_plant] * len(annotations)))
    
    elif file_format == 'FASTQ':
        annotations = decompress_annotation(
            compressed_data['annotation_patterns'],
            compressed_data['variable_parts']
        )

        if compressed_data.get('plus_compressed', False):
            plus_lines = decompress_annotation(
                compressed_data['plus_patterns'],
                compressed_data['plus_variables']
            )
        else:
            plus_lines = ['+'] * len(annotations)

        sequences = []
        for ann, (seq, qual_data), plus in zip(annotations, compressed_data['sequences'], plus_lines):
            sequences.append((ann, seq, plus, qual_data)) 
        
        return 'FASTQ', sequences

def dna_decompress_single_fasta_sequence(compressed_data: Tuple[str, str, bool]) -> Tuple[str, str]:
    annotation, compressed_sequence, is_plant = compressed_data
    if is_plant:
        decompressed_sequence = dna_decode_sequence_plant(compressed_sequence)
    else:
        decompressed_sequence = dna_decode_sequence(compressed_sequence)
    return annotation, decompressed_sequence


def dna_parallel_decompress_fasta(compressed_data: List[Tuple[str, str, bool]], 
                                processes: int = None,
                                is_plant: bool = False) -> List[Tuple[str, str]]:
    if len(compressed_data[0]) == 3:
        compressed_data_with_flag = compressed_data
    else:
        compressed_data_with_flag = [(ann, seq, is_plant) for ann, seq in compressed_data]
    
    with multiprocessing.Pool(processes=processes) as pool:
        return pool.map(dna_decompress_single_fasta_sequence, compressed_data_with_flag)

def dna_decompress_single_sequence(compressed_data: Tuple[str, str, str, Tuple[bytes, bytes, int], bool]) -> Tuple[str, str, str, str]:
    try:
        annotation, compressed_sequence, plus_line, quality_data, is_plant = compressed_data
        compressed_quality, metadata, encoding_type = quality_data

        if is_plant:
            decompressed_sequence = dna_decode_sequence_plant(compressed_sequence)
        else:
            decompressed_sequence = dna_decode_sequence(compressed_sequence)
        
        decompressed_quality = decode_quality_smart(compressed_quality, metadata, encoding_type)
        return annotation, decompressed_sequence, plus_line, decompressed_quality
    except Exception as e:
        print(f"Error in decompress_single_sequence: {str(e)}")
        raise

def dna_parallel_decompress(compressed_data: List[Tuple[str, bytes, str, Tuple[bytes, bytes, int]]], 
                          processes: int = None, 
                          is_plant: bool = False) -> List[Tuple[str, str, str, str]]:
    try:
        if len(compressed_data[0]) == 5: 
            data_with_plant = compressed_data
        else:  
            data_with_plant = [(ann, seq, plus, qual_data, is_plant) 
                              for ann, seq, plus, qual_data in compressed_data]
        
        with multiprocessing.Pool(processes=processes) as pool:
            results = pool.map(dna_decompress_single_sequence, data_with_plant)
        
        return results
    except Exception as e:
        print(f"Error in parallel_decompress: {str(e)}")
        raise

