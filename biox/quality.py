from dataclasses import dataclass
from collections import Counter
from typing import Any, Dict, Tuple
import pickle
import heapq

def analyze_quality_pattern(quality: str, sample_size: int = 100) -> float:
    """分析质量值模式
    返回重复率 (0-1之间，1表示全部相同)
    """
    if not quality:
        return 0.0
    
    # 取样本进行分析
    sample = quality[:min(sample_size, len(quality))]
    
    # 计算连续相同字符的长度
    runs = []
    current_char = sample[0]
    current_run = 1
    
    for char in sample[1:]:
        if char == current_char:
            current_run += 1
        else:
            runs.append(current_run)
            current_char = char
            current_run = 1
    runs.append(current_run)
    
    # 计算平均运行长度
    avg_run_length = sum(runs) / len(runs)
    
    # 返回标准化的得分 (1表示非常适合RLE，0表示不适合)
    return min(1.0, avg_run_length / 10.0)

def encode_quality_smart(quality: str) -> Tuple[bytes, bytes, int]:
    """智能质量值编码"""
    if not quality:
        return b'', b'', 0
    
    # 分析模式
    rle_score = analyze_quality_pattern(quality)
    
    if rle_score > 0.7:  # 如果有较多重复模式
        # 尝试RLE编码
        rle_data = rle_encode_quality(quality)
        
        # 计算压缩比
        rle_ratio = len(rle_data) / len(quality)

        if rle_ratio < 0.5:  # 如果RLE压缩效果好
            print("Quality encoding method: Run-Length Encoding")
            return rle_data, b'', 1
    
    # 默认使用Huffman编码
    print("Quality encoding method: Huffman Encoding")
    encoded, metadata = huffman_encode(quality)
    return encoded, metadata, 0

def decode_quality_smart(encoded: bytes, metadata: bytes, encoding_type: int) -> str:
    """智能解码质量值"""
    if encoding_type == 1:  # RLE
        return rle_decode_quality(encoded)
    else:  # Huffman
        return huffman_decode(encoded, metadata)

@dataclass
class HuffmanNode:
    char: str = ''
    freq: int = 0
    left: Any = None
    right: Any = None

    def __lt__(self, other):
        return self.freq < other.freq

def build_huffman_tree(text: str) -> HuffmanNode:
    """构建哈夫曼树"""
    frequency = Counter(text)
    heap = []
    for char, freq in frequency.items():
        heapq.heappush(heap, HuffmanNode(char=char, freq=freq))

    while len(heap) > 1:
        left = heapq.heappop(heap)
        right = heapq.heappop(heap)
        internal = HuffmanNode(freq=left.freq + right.freq, left=left, right=right)
        heapq.heappush(heap, internal)

    return heap[0] if heap else HuffmanNode()

def build_huffman_codes(node: HuffmanNode, code: str = '', codes: Dict[str, str] = None) -> Dict[str, str]:
    """生成哈夫曼编码表"""
    if codes is None:
        codes = {}

    if node.char:  # 叶子节点
        codes[node.char] = code if code else '0'
        return codes

    if node.left:
        build_huffman_codes(node.left, code + '0', codes)
    if node.right:
        build_huffman_codes(node.right, code + '1', codes)

    return codes

def huffman_encode(text: str) -> Tuple[bytes, bytes]:
    """哈夫曼编码"""
    # 构建哈夫曼树和编码表
    tree = build_huffman_tree(text)
    codes = build_huffman_codes(tree)

    # 编码文本
    encoded = ''.join(codes[char] for char in text)

    # 处理编码结果，确保长度是8的倍数
    padding = (8 - len(encoded) % 8) % 8
    encoded += '0' * padding

    # 转换为字节
    encoded_bytes = bytearray()
    for i in range(0, len(encoded), 8):
        encoded_bytes.append(int(encoded[i:i + 8], 2))

    # 序列化树结构用于解码
    tree_data = pickle.dumps((codes, padding))

    return bytes(encoded_bytes), tree_data

def huffman_decode(encoded_bytes: bytes, tree_data: bytes) -> str:
    """哈夫曼解码"""
    # 加载编码表和填充信息
    codes, padding = pickle.loads(tree_data)

    # 反转编码表用于解码
    decode_dict = {code: char for char, code in codes.items()}

    # 转换字节为二进制字符串
    binary = ''.join(format(byte, '08b') for byte in encoded_bytes)
    if padding:
        binary = binary[:-padding]

    # 解码
    decoded = []
    current = ''
    for bit in binary:
        current += bit
        if current in decode_dict:
            decoded.append(decode_dict[current])
            current = ''

    return ''.join(decoded) 

def rle_encode_quality(quality: str) -> bytes:
    """游程编码质量值
    格式：质量值,长度;质量值,长度;...
    """
    if not quality:
        return b''
    
    result = []
    current_char = quality[0]
    count = 1
    
    for char in quality[1:]:
        if char == current_char:
            count += 1
        else:
            result.append(f"{current_char},{count}")
            current_char = char
            count = 1
    
    # 添加最后一个运行
    result.append(f"{current_char},{count}")
    
    # 用分号连接所有运行
    return ';'.join(result).encode()

def rle_decode_quality(encoded: bytes) -> str:
    """解码RLE编码的质量值"""
    if not encoded:
        return ''
    
    result = []
    # 解码为字符串并按分号分割
    runs = encoded.decode().split(';')
    
    for run in runs:
        char, count = run.split(',')
        result.extend([char] * int(count))
    
    return ''.join(result)
