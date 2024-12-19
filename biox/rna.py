import lzma
import multiprocessing
from typing import List, Dict, Tuple, Any, Union
import pickle
from annotation import compress_annotation, decompress_annotation
from quality import  encode_quality_smart, decode_quality_smart

# 文件格式识别与解析
def detect_file_format(file_path: str) -> str:
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()
        if first_line.startswith('>'):
            return 'fasta'
        elif first_line.startswith('@'):
            return 'fastq'
        else:
            raise ValueError("无法识别的文件格式")

# 读取FASTA文件内容
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

def read_fastq(file_path: str) -> List[Tuple[str, str, str, str]]:  # 修改返回类型
    sequences = []
    with open(file_path, 'r') as file:
        current_record = []
        for line in file:
            line = line.strip()
            if not line:  # Skip empty lines
                continue
            if line.startswith('@') and len(current_record) == 0:
                current_record = [line]
            elif line.startswith('@') and len(current_record) == 4:
                sequences.append((current_record[0], current_record[1], current_record[2], current_record[3]))  # 添加+行
                current_record = [line]
            elif len(current_record) < 4:
                current_record.append(line)

        if len(current_record) == 4:
            sequences.append((current_record[0], current_record[1], current_record[2], current_record[3]))  # 添加+行

    return sequences

def encode_length(length: int) -> List[int]:
    """变长编码长度
    < 255: [0][8位长度]
    < 65535: [10][16位长度]
    其他: [11][32位长度]
    """
    bits = []
    if length < 255:
        bits.append(0)
        # 8位长度
        for i in range(7, -1, -1):
            bits.append((length >> i) & 1)
        return bits
    elif length < 65535:
        bits.extend([1, 0])
        # 16位长度
        for i in range(15, -1, -1):
            bits.append((length >> i) & 1)
        return bits
    else:
        bits.extend([1, 1])
        # 32位长度
        for i in range(31, -1, -1):
            bits.append((length >> i) & 1)
        return bits

def decode_length(bits: List[int], start: int) -> Tuple[int, int]:
    """解码变长长度，返回(长度, 新的开始位置)"""
    if bits[start] == 0:
        # 8位长度
        length = 0
        for i in range(8):
            length = (length << 1) | bits[start + 1 + i]
        return length, start + 9
    elif bits[start] == 1 and bits[start + 1] == 0:
        # 16位长度
        length = 0
        for i in range(16):
            length = (length << 1) | bits[start + 2 + i]
        return length, start + 18
    else:
        # 32位长度
        length = 0
        for i in range(32):
            length = (length << 1) | bits[start + 2 + i]
        return length, start + 34


def rna_encode_sequence(sequence: str) -> bytes:
    """DNA序列编码"""
    sequence = sequence.upper()
    
    # 第一步：预扫描，获取初始段表
    segments = []
    current_type = None  # None=开始, True=特殊字符段, False=ATCG段
    current_start = 0
    current_segment = []
    
    for i, char in enumerate(sequence):
        is_special = char not in 'AUCG'
        
        # 如果类型改变，保存当前段
        if current_type is not None and is_special != current_type:
            segments.append((current_type, current_start, ''.join(current_segment)))
            current_start = i
            current_segment = []
        
        current_type = is_special
        current_segment.append(char)
    
    # 添加最后一段
    if current_segment:
        segments.append((current_type, current_start, ''.join(current_segment)))
    
    # 第二步：合并短段（<10bp）
    merged_segments = []
    current_merged = None
    
    for segment in segments:
        is_special, start, data = segment
        
        if current_merged is None:
            current_merged = list(segment)
        elif current_merged[0] == is_special and len(current_merged[2]) < 10:
            # 合并相同类型的短段
            current_merged[2] += data
        else:
            # 保存当前合并段，开始新的段
            merged_segments.append(tuple(current_merged))
            current_merged = list(segment)
    
    # 添加最后一个合并段
    if current_merged:
        merged_segments.append(tuple(current_merged))
    
    # 第三步：编码
    bits = []
    
    # 添加总段数（32位）
    total_segments = len(merged_segments)
    for i in range(31, -1, -1):
        bits.append((total_segments >> i) & 1)
    
    # 编码每个段
    for is_special, _, data in merged_segments:
        # 添加类型标记（1位）
        bits.append(1 if is_special else 0)
        
        # 添加长度（变长编码）
        bits.extend(encode_length(len(data)))
        
        # 添加数据
        if is_special:
            # 特殊字符：8位ASCII
            for char in data:
                ascii_val = ord(char)
                for i in range(8):
                    bits.append((ascii_val >> (7-i)) & 1)
        else:
            # ATCG：2位编码
            for char in data:
                if char == 'A':
                    bits.extend([0, 0])
                elif char == 'U':
                    bits.extend([1, 1])
                elif char == 'C':
                    bits.extend([0, 1])
                elif char == 'G':
                    bits.extend([1, 0])
    
    # 填充到8位的倍数
    while len(bits) % 8 != 0:
        bits.append(0)
    
    # 转换为字节
    result_bytes = bytearray()
    for i in range(0, len(bits), 8):
        byte = 0
        for j in range(8):
            if i + j < len(bits):
                byte = (byte << 1) | bits[i + j]
        result_bytes.append(byte)
    
    return bytes(result_bytes)

def rna_decode_sequence(compressed: bytes) -> str:
    """RNA序列解码"""
    decoded = []
    
    # 转换为位列表
    bits = []
    for byte in compressed:
        for i in range(7, -1, -1):
            bits.append((byte >> i) & 1)
    
    # 读取总段数
    total_segments = 0
    for i in range(32):
        total_segments = (total_segments << 1) | bits[i]
    i = 32
    
    # 解码每个段
    for _ in range(total_segments):
        if i >= len(bits):
            break
        
        # 读取类型标记
        is_special = bits[i]
        i += 1
        
        # 读取长度
        length, i = decode_length(bits, i)
        
        if is_special:
            # 特殊字符段
            for _ in range(length):
                if i + 8 > len(bits):
                    break
                ascii_val = 0
                for j in range(8):
                    ascii_val = (ascii_val << 1) | bits[i + j]
                decoded.append(chr(ascii_val))
                i += 8
        else:
            # ATCG段
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


# 修改为DNA编码
RNA_CODES = {
    # 标准碱基
    'A': 0b00000000,
    'U': 0b00000001,
    'C': 0b00000010,
    'G': 0b00000011,
    
    # 模糊碱基
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

# 设置分块大小为8GB
CHUNK_SIZE = 8 * 1024 * 1024 * 1024

def rna_encode_sequence_plant(sequence: str) -> bytes:
    """RNA序列编码"""
    sequence = sequence.upper()  # 转换为大写
    result_bytes = bytearray()
    for char in sequence:
        if char in RNA_CODES:
            result_bytes.append(RNA_CODES[char])
        else:
            result_bytes.append(ord(char))
    return bytes(result_bytes)

def rna_decode_sequence_plant(compressed: bytes) -> str:
    """DNA序列解码"""
    # 创建反向映射
    reverse_mapping = {code: base for base, code in RNA_CODES.items()}
    return ''.join(reverse_mapping.get(byte, chr(byte)) for byte in compressed)

def process_large_sequence(sequence: str) -> bytes:
    """处理大型序列"""
    if len(sequence) <= CHUNK_SIZE:
        return rna_encode_sequence_plant(sequence)
    
    result = bytearray()
    # 分块处理
    for i in range(0, len(sequence), CHUNK_SIZE):
        chunk = sequence[i:i + CHUNK_SIZE]
        result.extend(rna_encode_sequence_plant(chunk))
    
    return bytes(result)

def rna_compress_single_fasta_sequence(args: Tuple[str, str, bool]) -> Tuple[str, str]:
    """压缩单个FASTA序列"""
    annotation, sequence, is_plant = args
    if is_plant:
        compressed_sequence = process_large_sequence(sequence)  # 使用植物基因组方案
    else:
        compressed_sequence = rna_encode_sequence(sequence)    # 使用默认方案
    return annotation, compressed_sequence

def rna_parallel_compress_fasta(sequences: List[Tuple[str, str]], 
                              processes: int = None,
                              is_plant: bool = False) -> List[Tuple[str, str]]:
    """并行压缩处理FASTA"""
    processes = processes or multiprocessing.cpu_count()
    # 为每个序列添加is_plant参数
    sequences_with_flag = [(ann, seq, is_plant) for ann, seq in sequences]
    with multiprocessing.Pool(processes=processes) as pool:
        return pool.map(rna_compress_single_fasta_sequence, sequences_with_flag)

def rna_write_compressed_fastq(data: List[Tuple[str, str, str, Tuple[bytes, bytes, int]]], output_file: str, compress_plus: bool = False, preset: int = 3):
    """写入压缩的FASTQ文件"""

    # 提取所有注释行
    annotations = [item[0] for item in data]
    # 压缩注释行
    annotation_patterns, variable_parts = compress_annotation(annotations)

    # 如果需要压缩+行
    plus_patterns = None
    plus_variables = None
    if compress_plus:
        plus_lines = [item[2] for item in data]
        plus_patterns, plus_variables = compress_annotation(plus_lines)

    # 准备完整的压缩数据
    sequences_data = [(seq, qual_data) for _, seq, _, qual_data in data]
 
    compressed_data = {
        'format': 'FASTQ',
        'annotation_patterns': annotation_patterns,
        'variable_parts': variable_parts,
        'plus_compressed': compress_plus,
        'plus_patterns': plus_patterns,
        'plus_variables': plus_variables,
        'sequences': sequences_data
    }

    # 使用lzma压缩
    with open(output_file, 'wb') as f_out:
        compressed = lzma.compress(pickle.dumps(compressed_data), preset=preset)
        f_out.write(compressed)


def rna_write_compressed_fasta(data: List[Tuple[str, str]], 
                             output_file: str, 
                             preset: int = 3,
                             is_plant: bool = False):
    """写入压缩的FASTA文件"""
    # 提取所有注释行
    annotations = [item[0] for item in data]
    # 压缩注释行
    annotation_patterns, variable_parts = compress_annotation(annotations)

    # 准备完整的压缩数据
    compressed_data = {
        'format': 'FASTA',
        'is_plant': is_plant,  # 添加压缩方案标记
        'annotation_patterns': annotation_patterns,
        'variable_parts': variable_parts,
        'sequences': [seq for _, seq in data]
    }

    # 使用lzma压缩
    with open(output_file, 'wb') as f_out:
        compressed = lzma.compress(pickle.dumps(compressed_data), preset=preset)
        f_out.write(compressed)

def rna_read_compressed_file(file_path: str) -> Tuple[str, Union[List[Tuple[str, str]], List[Tuple[str, str, str, Tuple[bytes, bytes, int]]]]]:
    """读取压缩文件"""
    with open(file_path, 'rb') as f_in:
        compressed_data = pickle.loads(lzma.decompress(f_in.read()))
    
    file_format = compressed_data.get('format')
    is_plant = compressed_data.get('is_plant', False)  # 获取压缩方案标记
    
    if not file_format:
        raise ValueError("无法识别的文件格式")

    if file_format == 'FASTA':
        # 解压注释行
        annotations = decompress_annotation(
            compressed_data['annotation_patterns'],
            compressed_data['variable_parts']
        )

        # 重建完整结构，包含压缩方案信息
        return 'FASTA', list(zip(annotations, compressed_data['sequences'], [is_plant] * len(annotations)))
    
    elif file_format == 'FASTQ':
        # 解压注释行
        annotations = decompress_annotation(
            compressed_data['annotation_patterns'],
            compressed_data['variable_parts']
        )

        # 解压+行
        if compressed_data.get('plus_compressed', False):
            plus_lines = decompress_annotation(
                compressed_data['plus_patterns'],
                compressed_data['plus_variables']
            )
        else:
            plus_lines = ['+'] * len(annotations)

        # 重建完整数据结构，确保包含质量值数据
        return 'FASTQ', [
            (ann, seq, plus_line, qual_data)
            for ann, (seq, qual_data), plus_line in zip(annotations, compressed_data['sequences'], plus_lines)
        ]
    
    else:
        raise ValueError(f"不支持的文件格式: {file_format}")


def rna_decompress_single_fasta_sequence(compressed_data: Tuple[str, str, bool]) -> Tuple[str, str]:
    """解压单个FASTA序列"""
    annotation, compressed_sequence, is_plant = compressed_data
    if is_plant:
        decompressed_sequence = rna_decode_sequence_plant(compressed_sequence)
    else:
        decompressed_sequence = rna_decode_sequence(compressed_sequence)
    return annotation, decompressed_sequence


def rna_parallel_decompress_fasta(compressed_data: List[Tuple[str, str]], 
                                processes: int = None,
                                is_plant: bool = False) -> List[Tuple[str, str]]:
    """并行解压处理FASTA"""
    # 为每个压缩数据添加is_plant参数
    compressed_data_with_flag = [(ann, seq, is_plant) for ann, seq in compressed_data]
    with multiprocessing.Pool(processes=processes) as pool:
        return pool.map(rna_decompress_single_fasta_sequence, compressed_data_with_flag)

def rna_decompress_single_sequence(compressed_data: Tuple[str, str, str, Tuple[bytes, bytes, int]]) -> Tuple[str, str, str, str]:
    """解压单个序列"""
    try:
        
        annotation, compressed_sequence, plus_line, quality_data = compressed_data
        compressed_quality, metadata, encoding_type = quality_data  # 解包质量值数据
        
        decompressed_sequence = rna_decode_sequence(compressed_sequence)
        decompressed_quality = decode_quality_smart(compressed_quality, metadata, encoding_type)
        
        return annotation, decompressed_sequence, plus_line, decompressed_quality
    except Exception as e:
        print(f"Error in decompress_single_sequence: {str(e)}")
        raise

def rna_parallel_decompress(compressed_data: List[Tuple[str, bytes, str, Tuple[bytes, bytes, int]]], processes: int = None) -> List[Tuple[str, str, str, str]]:
    """并行解压处理"""
    try:
        
        with multiprocessing.Pool(processes=processes) as pool:
            results = pool.map(rna_decompress_single_sequence, compressed_data)
        
        return results
    except Exception as e:
        print(f"Error in parallel_decompress: {str(e)}")
        raise


def rna_compress_single_sequence(args: Tuple[str, str, str, str]) -> Tuple[str, str, str, Tuple[bytes, bytes, int]]:
    """压缩单个FASTQ序列"""
    annotation, sequence, plus_line, quality = args
    # 使用分块处理大序列
    compressed_sequence = process_large_sequence(sequence)
    compressed_quality, metadata, encoding_type = encode_quality_smart(quality)
    return annotation, compressed_sequence, plus_line, (compressed_quality, metadata, encoding_type)

def rna_parallel_compress(sequences: List[Tuple[str, str, str, str]], 
                        processes: int = None) -> List[Tuple[str, str, str, Tuple[bytes, bytes, int]]]:
    """并行压缩处理FASTQ"""
    processes = processes or multiprocessing.cpu_count()
    with multiprocessing.Pool(processes=processes) as pool:
        return pool.map(rna_compress_single_sequence, sequences)