from typing import List, Dict, Tuple
import os

def compress_annotation(annotations: List[str]) -> Tuple[Dict[str, str], List[str]]:
    if not annotations:
        return {}, []

    prefix = os.path.commonprefix(annotations)
    stripped_prefix = [ann[len(prefix):] for ann in annotations]

    reversed_anns = [s[::-1] for s in stripped_prefix]

    suffix = os.path.commonprefix(reversed_anns)[::-1]

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
    
    prefix = patterns.get('prefix', '')
    suffix = patterns.get('suffix', '')

    return [f"{prefix}{var_part}{suffix}" for var_part in variable_parts]