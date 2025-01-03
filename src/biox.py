import argparse
import os
import time
from .dna import (detect_file_format, read_fasta, read_fastq,
                dna_parallel_compress_fasta, dna_parallel_compress,
                dna_write_compressed_fasta, dna_write_compressed_fastq,
                dna_read_compressed_file, dna_parallel_decompress_fasta,
                dna_parallel_decompress)
from .rna import (rna_parallel_compress_fasta, rna_parallel_compress,
                rna_write_compressed_fasta, rna_write_compressed_fastq,
                rna_read_compressed_file, rna_parallel_decompress_fasta,
                rna_parallel_decompress)

from .protein import (
    protein_parallel_compress, protein_parallel_decompress,
    protein_write_compressed, protein_read_compressed,
    read_fasta as protein_read_fasta
)

from .file import file_compress, file_decompress

def main():
    parser = argparse.ArgumentParser(description='Biological Sequence Compression Tool')
    parser.add_argument('input_file', help='Input file path (FASTQ/FASTA)')
    parser.add_argument('-c', '--compress', action='store_true', help='Compression mode')
    parser.add_argument('-d', '--decompress', action='store_true', help='Decompression mode')
    parser.add_argument('-t', '--type', choices=['dna', 'rna', 'protein', 'file'], 
                       required=True, help='Sequence type (dna/rna/protein) or regular file')
    parser.add_argument('-l', '--level', type=int, choices=range(1, 10), 
                       default=9, help='Compression level (1-9, default: 9)')
    parser.add_argument('-ps', '--plus_line', choices=['ignore', 'compress'], 
                       default='ignore', help='FASTQ plus line handling')
    parser.add_argument('-n', '--num_processes', type=int, default=None, 
                       help=f'Number of parallel processes')
    parser.add_argument('-p', '--plant', action='store_true', help='Use plant genome compression scheme')
    parser.add_argument('-s', '--split', type=int, choices=range(2, 11),
                       help='Split output into N volumes (2-10)')
    parser.add_argument('-o', '--output', help='Output file path (default: input_file.biox)')
    
    args = parser.parse_args()

    if not args.output:
        if args.decompress:
            if args.input_file.endswith('.biox'):
                args.output = args.input_file[:-5] + '.decoded'
            elif args.input_file.endswith('.biox.001'):
                args.output = args.input_file[:-9] + '.decoded'
            else:
                args.output = args.input_file + '.decoded'
        else:
            args.output = args.input_file + '.biox'

    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' does not exist")
        return

    if not args.compress and not args.decompress:
        parser.error("Must specify either -c (compress) or -d (decompress) mode")

    try:
        if args.type == 'file':
            if args.compress:
                print(f"\nCompressing file: {args.input_file}")
                original_size = os.path.getsize(args.input_file)
                print(f"Original file size: {original_size / 1024 / 1024:.2f} MB")
                
                start_time = time.time()
                file_compress(args.input_file, args.output, preset=args.level)
                end_time = time.time()
                compress_time = end_time - start_time
                
                compressed_size = os.path.getsize(args.output)
                compression_ratio = compressed_size / original_size
                
                print(f"\nCompression completed:")
                print(f"- Time taken: {compress_time:.2f} seconds")
                print(f"- Compressed size: {compressed_size / 1024 / 1024:.2f} MB")
                print(f"- Compression ratio: {compression_ratio:.2f}x")
                
            else:  
                print(f"\nDecompressing file: {args.input_file}")
                start_time = time.time()
                
                file_decompress(args.input_file, args.output)
                
                end_time = time.time()
                print(f"\nDecompression completed! Time taken: {end_time - start_time:.2f} seconds")
            
            return

        if args.compress:
            file_format = detect_file_format(args.input_file)
            print(f"Detected file format: {file_format}")

            original_size = os.path.getsize(args.input_file)
            print(f"Original file size: {original_size / 1024 / 1024:.2f} MB")

            if args.type == 'protein':
                if file_format != 'fasta':
                    print("Error: Protein sequences only support FASTA format")
                    return
                
                print(f"\nReading FASTA file: {args.input_file}")
                sequences = protein_read_fasta(args.input_file)
                total_sequences = len(sequences)
                print(f"Total sequences: {total_sequences}")
                
                print("\nStarting compression...")
                start_time = time.time()
                compressed_data = protein_parallel_compress(sequences, processes=args.num_processes)
                protein_write_compressed(compressed_data, args.output, preset=args.level)
                end_time = time.time()
                compress_time = end_time - start_time

            else:
                handlers = {
                    'dna': {
                        'fastq': (dna_parallel_compress, dna_write_compressed_fastq),
                        'fasta': (dna_parallel_compress_fasta, dna_write_compressed_fasta)
                    },
                    'rna': {
                        'fastq': (rna_parallel_compress, rna_write_compressed_fastq),
                        'fasta': (rna_parallel_compress_fasta, rna_write_compressed_fasta)
                    }
                }

                compress_plus = (args.plus_line == 'compress')
                try:
                    parallel_compress, write_compressed = handlers[args.type][file_format]
                except KeyError:
                    print(f"Error: Unsupported file format: {file_format} for {args.type}")
                    return

                if file_format == 'fastq':
                    print(f"Reading FASTQ file: {args.input_file}")
                    sequences = read_fastq(args.input_file)
                else:
                    print(f"Reading FASTA file: {args.input_file}")
                    sequences = read_fasta(args.input_file)

                total_sequences = len(sequences)
                print(f"Total sequences: {total_sequences}")

                print("\nStarting compression...")
                start_time = time.time()
                if args.type == 'dna' or args.type == 'rna':
                    compressed_data = parallel_compress(sequences, 
                                                     processes=args.num_processes,
                                                     is_plant=args.plant)
                else:
                    compressed_data = parallel_compress(sequences, 
                                                     processes=args.num_processes)
                
                if file_format == 'fastq':
                    write_compressed(compressed_data, args.output, 
                                  compress_plus=compress_plus, preset=args.level,
                                  is_plant=args.plant if args.type in ['dna', 'rna'] else False,
                                  num_volumes=args.split)
                else:
                    write_compressed(compressed_data, args.output, 
                                  preset=args.level,
                                  is_plant=args.plant if args.type in ['dna', 'rna'] else False,
                                  num_volumes=args.split)
                
                end_time = time.time()
                compress_time = end_time - start_time

                if args.split:
                    compressed_size = 0
                    for i in range(1, args.split + 1):
                        volume_path = f"{args.output}.{i:03d}"
                        if os.path.exists(volume_path):
                            compressed_size += os.path.getsize(volume_path)
                else:
                    compressed_size = os.path.getsize(args.output)

            compression_ratio = compressed_size / original_size 

            print(f"\nCompression completed:")
            print(f"- Time taken: {compress_time:.2f} seconds")
            print(f"- Compressed size: {compressed_size / 1024 / 1024:.2f} MB")
            print(f"- Compression ratio: {compression_ratio:.2f}x")
            print(f"- Processing speed: {total_sequences / compress_time:.2f} sequences/second")
            
            if args.split:
                print("\nVolume sizes:")
                for i in range(1, args.split + 1):
                    volume_path = f"{args.output}.{i:03d}"
                    if os.path.exists(volume_path):
                        volume_size = os.path.getsize(volume_path)
                        print(f"- Volume {i:03d}: {volume_size / 1024 / 1024:.2f} MB")

        else:
            print(f"\nReading compressed file: {args.input_file}")
            start_time = time.time()

            if args.type == 'protein':
                compressed_data = protein_read_compressed(args.input_file)
                decompressed_data = protein_parallel_decompress(compressed_data, processes=args.num_processes)
                file_format = 'fasta'
            else:
                read_func = dna_read_compressed_file if args.type == 'dna' else rna_read_compressed_file
                file_format, compressed_data = read_func(args.input_file)
                
                if file_format.endswith('FASTQ'):
                    decompress_func = dna_parallel_decompress if args.type == 'dna' else rna_parallel_decompress
                    
                    decompressed_data = decompress_func(compressed_data, 
                                                      processes=args.num_processes,
                                                      is_plant=args.plant)
                    file_format = 'fastq'
                else:
                    decompress_func = dna_parallel_decompress_fasta if args.type == 'dna' else rna_parallel_decompress_fasta
                    is_plant = compressed_data[0][2] if len(compressed_data[0]) == 3 else args.plant
                    decompressed_data = decompress_func(compressed_data, 
                                                      processes=args.num_processes,
                                                      is_plant=is_plant)
                    file_format = 'fasta'
                    
            print(f"\nWriting {file_format.upper()} file: {args.output}")
            with open(args.output, 'w') as f:
                if file_format == 'fastq':
                    for annotation, sequence, plus_line, quality in decompressed_data:
                        f.write(f"{annotation}\n{sequence}\n{plus_line}\n{quality}\n")
                else:
                    for annotation, sequence in decompressed_data:
                        f.write(f"{annotation}\n{sequence}\n")

            end_time = time.time()
            print(f"\nDecompression completed! Time taken: {end_time - start_time:.2f} seconds")

    except Exception as e:
        print(f"\nError: {str(e)}")
        return

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nProgram interrupted by user")
    except Exception as e:
        print(f"\nProgram error: {str(e)}")
    finally:
        print("\nProgram finished")
