import lzma
import multiprocessing
from typing import List, Dict, Tuple, Any, Union
import pickle
from annotation import compress_annotation, decompress_annotation

# 按物理化学特性对氨基酸分组并编码为8位
AMINO_ACID_GROUPS = {
    # 疏水性
    'A': 0b00000000,
    'I': 0b00000001,
    'L': 0b00000010,
    'M': 0b00000011,
    'F': 0b00000100,
    'W': 0b00000101,
    'V': 0b00000110,
    
    # 极性非带电
    'S': 0b00001000,
    'T': 0b00001001,
    'N': 0b00001010,
    'Q': 0b00001011,
    
    # 带正电
    'K': 0b00010000,
    'R': 0b00010001,
    'H': 0b00010010,
    
    # 带负电
    'D': 0b00011000,
    'E': 0b00011001,
    
    # 特殊结构
    'Y': 0b00100000,
    'G': 0b00100001,
    'C': 0b00100010,
    'P': 0b00100011,
    
    # 终止密码子
    '*': 0b00101000,
    
    # 其他字母
    'B': 0b00110000,  
    'J': 0b00110001,
    'O': 0b00110010,
    'U': 0b00110011,  
    'X': 0b00110100,  
    'Z': 0b00110101,  
    
    # 间隙符号
    '-': 0b00111000,
}

def read_fasta(file_path: str) -> List[Tuple[str, str]]:
    """读取FASTA文件"""
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
    """压缩蛋白质序列"""
    # 预处理：转大写，未知字符转为X
    sequence = ''.join(char.upper() if char.upper() in AMINO_ACID_GROUPS else 'X' for char in sequence)
    
    # 直接转换为字节
    result_bytes = bytearray()
    for char in sequence:
        result_bytes.append(AMINO_ACID_GROUPS[char])
    
    return bytes(result_bytes)

def protein_decode_sequence(compressed: bytes) -> str:
    """解压蛋白质序列"""
    # 创建反向映射
    reverse_mapping = {code: aa for aa, code in AMINO_ACID_GROUPS.items()}
    
    # 直接从字节解码
    return ''.join(reverse_mapping.get(byte, 'A') for byte in compressed)

def protein_compress_single(args: Tuple[str, str]) -> Tuple[str, bytes]:
    """压缩单个蛋白质序列"""
    annotation, sequence = args
    compressed_sequence = protein_encode_sequence(sequence)
    return annotation, compressed_sequence

def protein_decompress_single(compressed_data: Tuple[str, bytes]) -> Tuple[str, str]:
    """解压单个蛋白质序列"""
    annotation, compressed_sequence = compressed_data
    decompressed_sequence = protein_decode_sequence(compressed_sequence)
    return annotation, decompressed_sequence

def protein_parallel_compress(sequences: List[Tuple[str, str]], processes: int = None) -> List[Tuple[str, bytes]]:
    """并行压缩处理蛋白质序列"""
    with multiprocessing.Pool(processes=processes) as pool:
        return pool.map(protein_compress_single, sequences)

def protein_parallel_decompress(compressed_data: List[Tuple[str, bytes]], processes: int = None) -> List[Tuple[str, str]]:
    """并行解压处理蛋白质序列"""
    with multiprocessing.Pool(processes=processes) as pool:
        return pool.map(protein_decompress_single, compressed_data)

def protein_write_compressed(data: List[Tuple[str, bytes]], output_file: str, preset: int = 3):
    """写入压缩的蛋白质序列文件"""
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
    """读取压缩的蛋白质序列文件"""
    with open(file_path, 'rb') as f_in:
        compressed_data = pickle.loads(lzma.decompress(f_in.read()))
    
    if compressed_data.get('format') != 'PROTEIN_FASTA':
        raise ValueError("Unsupported file format")

    annotations = decompress_annotation(
        compressed_data['annotation_patterns'],
        compressed_data['variable_parts']
    )

    return list(zip(annotations, compressed_data['sequences']))
