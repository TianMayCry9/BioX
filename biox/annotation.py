from typing import List, Dict, Tuple
import os

def compress_annotation(annotations: List[str]) -> Tuple[Dict[str, str], List[str]]:
    """压缩注释行，识别共同模式并提取可变部分"""
    # 找到共同前缀和后缀
    if not annotations:
        return {}, []

    # 找最长公共前缀
    prefix = os.path.commonprefix(annotations)
    # 去除前缀后的字符串
    stripped_prefix = [ann[len(prefix):] for ann in annotations]

    # 找最长公共后缀
    reversed_anns = [s[::-1] for s in stripped_prefix]
    suffix = os.path.commonprefix(reversed_anns)[::-1]

    # 提取变化部分
    variable_parts = []
    suffix_len = len(suffix)
    for ann in stripped_prefix:
        if suffix_len:
            variable_part = ann[:-suffix_len]
        else:
            variable_part = ann
        variable_parts.append(variable_part)

    patterns = {
        'prefix': prefix,
        'suffix': suffix
    }

    return patterns, variable_parts


def decompress_annotation(patterns: Dict[str, str], variable_parts: List[str]) -> List[str]:
    """重建完整的注释行"""
    prefix = patterns.get('prefix', '')
    suffix = patterns.get('suffix', '')

    return [f"{prefix}{var_part}{suffix}" for var_part in variable_parts]