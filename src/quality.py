from dataclasses import dataclass
from collections import Counter
from typing import Any, Dict, Tuple
import pickle
import heapq

def analyze_quality_pattern(quality: str, sample_size: int = 100) -> float:
    if not quality:
        return 0.0

    sample = quality[:min(sample_size, len(quality))]

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

    avg_run_length = sum(runs) / len(runs)

    return min(1.0, avg_run_length / 10.0)

def encode_quality_smart(quality: str) -> Tuple[bytes, bytes, int]:
    if not quality:
        return b'', b'', 0

    rle_score = analyze_quality_pattern(quality)
    
    if rle_score > 0.7: 
        rle_data = rle_encode_quality(quality)

        rle_ratio = len(rle_data) / len(quality)

        if rle_ratio < 0.5:  
            return rle_data, b'', 1

    encoded, metadata = huffman_encode(quality)
    return encoded, metadata, 0

def decode_quality_smart(encoded: bytes, metadata: bytes, encoding_type: int) -> str:
    if encoding_type == 1:  
        return rle_decode_quality(encoded)
    else:  
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
    if codes is None:
        codes = {}

    if node.char:  
        codes[node.char] = code if code else '0'
        return codes

    if node.left:
        build_huffman_codes(node.left, code + '0', codes)
    if node.right:
        build_huffman_codes(node.right, code + '1', codes)

    return codes

def huffman_encode(text: str) -> Tuple[bytes, bytes]:
    tree = build_huffman_tree(text)
    codes = build_huffman_codes(tree)

    encoded = ''.join(codes[char] for char in text)

    padding = (8 - len(encoded) % 8) % 8
    encoded += '0' * padding

    encoded_bytes = bytearray()
    for i in range(0, len(encoded), 8):
        encoded_bytes.append(int(encoded[i:i + 8], 2))

    tree_data = pickle.dumps((codes, padding))

    return bytes(encoded_bytes), tree_data

def huffman_decode(encoded_bytes: bytes, tree_data: bytes) -> str:
    codes, padding = pickle.loads(tree_data)

    decode_dict = {code: char for char, code in codes.items()}

    binary = ''.join(format(byte, '08b') for byte in encoded_bytes)
    if padding:
        binary = binary[:-padding]

    decoded = []
    current = ''
    for bit in binary:
        current += bit
        if current in decode_dict:
            decoded.append(decode_dict[current])
            current = ''

    return ''.join(decoded) 

def rle_encode_quality(quality: str) -> bytes:
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

    result.append(f"{current_char},{count}")
   
    return ';'.join(result).encode()

def rle_decode_quality(encoded: bytes) -> str:
    if not encoded:
        return ''
    
    result = []
    runs = encoded.decode().split(';')
    
    for run in runs:
        char, count = run.split(',')
        result.extend([char] * int(count))
    
    return ''.join(result)
