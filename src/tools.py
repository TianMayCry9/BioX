import os
from itertools import combinations
from .dna import (
    read_fasta as dna_read_fasta,
    dna_parallel_compress_fasta,
    dna_write_compressed_fasta
)
from .rna import (
    read_fasta as rna_read_fasta,
    rna_parallel_compress_fasta,
    rna_write_compressed_fasta
)
from .protein import (
    read_fasta as protein_read_fasta,
    protein_parallel_compress,
    protein_write_compressed
)

def get_sequence_files(input_dir):
    sequence_files = []
    valid_extensions = ('.fasta', '.fa', '.fna')
    
    try:
        for file_name in os.listdir(input_dir):
            if file_name.endswith(valid_extensions):
                if os.path.isfile(os.path.join(input_dir, file_name)):
                    sequence_files.append(file_name)
        return sequence_files
    except Exception as e:
        print(f"Error reading directory: {e}")
        return []

def add_class_prefix(input_dir):
    renamed_files = []
    valid_extensions = ('.fasta', '.fa', '.fna')
    
    try:
        for file_name in os.listdir(input_dir):
            if file_name.endswith(valid_extensions):
                if file_name.startswith('class_'):
                    renamed_files.append(file_name)
                    continue
                    
                file_path = os.path.join(input_dir, file_name)
                if os.path.isfile(file_path):
                    new_file_name = f"class_{file_name}"
                    new_file_path = os.path.join(input_dir, new_file_name)
                    os.rename(file_path, new_file_path)
                    renamed_files.append(new_file_name)

        if renamed_files:
            print("All files prefixed with class_")
        return renamed_files

    except Exception as e:
        print(f"Error during renaming: {e}")
        return []

def get_file_extension(filename):
    for ext in ('.fasta', '.fa', '.fna'):
        if filename.endswith(ext):
            return ext
    return ''

def pairwise_merge_fasta(class_files, input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    file_combinations = combinations(class_files, 2)
    merged_count = 0

    for file1, file2 in file_combinations:
        try:
            ext1 = get_file_extension(file1)
            ext2 = get_file_extension(file2)
            
            species1 = file1[6:-len(ext1)]
            species2 = file2[6:-len(ext2)]

            file1_path = os.path.join(input_dir, file1)
            file2_path = os.path.join(input_dir, file2)

            if not os.path.exists(file1_path):
                print(f"Error: File not found {file1_path}")
                continue
            if not os.path.exists(file2_path):
                print(f"Error: File not found {file2_path}")
                continue

            output_filename = f"{species1}-{species2}{ext1}"
            output_path = os.path.join(output_dir, output_filename)

            if os.path.exists(output_path):
                continue

            with open(file1_path, 'r') as f1, open(file2_path, 'r') as f2:
                content1 = f1.read()
                content2 = f2.read()

            with open(output_path, 'w') as out_file:
                out_file.write(content1 + "\n" + content2)

            merged_count += 1

        except Exception as e:
            print(f"Error merging files ({file1}, {file2}): {e}")
            continue

    print(f"\nSuccessfully merged {merged_count} sequence pairs")

def batch_compress_sequences(input_dir, output_dir, sequence_type='dna'):
    os.makedirs(output_dir, exist_ok=True)
    
    compression_handlers = {
        'dna': {
            'read': dna_read_fasta,
            'compress': dna_parallel_compress_fasta,
            'write': dna_write_compressed_fasta,
        },
        'rna': {
            'read': rna_read_fasta,
            'compress': rna_parallel_compress_fasta,
            'write': rna_write_compressed_fasta,
        },
        'protein': {
            'read': protein_read_fasta,
            'compress': protein_parallel_compress,
            'write': protein_write_compressed,
        }
    }
    
    handler = compression_handlers.get(sequence_type)
    if not handler:
        raise ValueError(f"Unsupported sequence type: {sequence_type}")
    
    compressed_count = 0
    error_count = 0
    
    for file_name in os.listdir(input_dir):
        if file_name.endswith(('.fasta', '.fa', '.fna')):
            input_path = os.path.join(input_dir, file_name)
            output_path = os.path.join(output_dir, f"{file_name}.biox")
            
            if os.path.exists(output_path):
                continue
                
            try:
                sequences = handler['read'](input_path)
                if sequence_type == 'protein':
                    compressed_data = handler['compress'](sequences)
                    handler['write'](compressed_data, output_path, preset=3)
                else:
                    compressed_data = handler['compress'](sequences)
                    handler['write'](compressed_data, output_path, preset=3)
                compressed_count += 1
                
            except Exception as e:
                print(f"Compression failed: {file_name}, error: {e}")
                error_count += 1
                continue
    
    print(f"\nCompression completed:")
    print(f"- Success: {compressed_count} files")
    print(f"- Failed: {error_count} files")

def cleanup_temp_files(input_dir):
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.biox') or '-' in file_name:
            file_path = os.path.join(input_dir, file_name)
            try:
                os.remove(file_path)
            except Exception as e:
                print(f"Error removing temporary file {file_name}: {e}")

def restore_original_names(input_dir):
    valid_extensions = ('.fasta', '.fa', '.fna')
    for file_name in os.listdir(input_dir):
        if file_name.startswith('class_') and file_name.endswith(valid_extensions):
            old_path = os.path.join(input_dir, file_name)
            new_name = file_name[6:]  
            new_path = os.path.join(input_dir, new_name)
            try:
                os.rename(old_path, new_path)
            except Exception as e:
                print(f"Error restoring original name for {file_name}: {e}")

def process_all(input_dir, output_dir, sequence_type='dna'):
    sequence_files = get_sequence_files(input_dir)
    if not sequence_files:
        print("No valid sequence files found")
        return False

    print(f"Found {len(sequence_files)} sequence files")

    print("\nProcessing class_ prefix...")
    class_files = add_class_prefix(input_dir)
    if not class_files:
        print("Failed to process class_ prefix")
        return False

    print("\nCompressing original sequence files...")
    for file_name in class_files:
        input_path = os.path.join(input_dir, file_name)
        output_path = os.path.join(input_dir, f"{file_name}.biox")
        batch_compress_sequences(input_dir, input_dir, sequence_type)

    print("\nMerging sequence files...")
    pairwise_merge_fasta(class_files, input_dir, input_dir)

    print("\nCompressing merged sequence files...")
    batch_compress_sequences(input_dir, input_dir, sequence_type)

    print("\nAll processing completed!")
    return True

def finally_cleanup(input_dir):
        cleanup_temp_files(input_dir)
        restore_original_names(input_dir)
