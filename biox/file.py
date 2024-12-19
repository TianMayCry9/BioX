import lzma
import os
from typing import List, Tuple

def file_compress(input_file: str, output_file: str, preset: int = 3) -> None:
    """Compress a regular file using LZMA"""
    with open(input_file, 'rb') as f_in:
        data = f_in.read()
    
    # 使用lzma压缩
    with open(output_file, 'wb') as f_out:
        compressed = lzma.compress(data, preset=preset)
        f_out.write(compressed)

def file_decompress(input_file: str, output_file: str) -> None:
    """Decompress a file compressed with LZMA"""
    with open(input_file, 'rb') as f_in:
        compressed_data = f_in.read()
    
    # 解压数据
    with open(output_file, 'wb') as f_out:
        decompressed = lzma.decompress(compressed_data)
        f_out.write(decompressed) 