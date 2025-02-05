import argparse
import os
import time
from sklearn.metrics import f1_score
from .dna import (
    detect_file_format, read_fasta, read_fastq,
    dna_parallel_compress, dna_parallel_decompress,
    dna_parallel_compress_fasta, dna_parallel_decompress_fasta,
    dna_write_compressed_fastq, dna_write_compressed_fasta,
    dna_read_compressed_file
)

from .rna import (
    rna_parallel_compress, rna_parallel_decompress,
    rna_parallel_compress_fasta, rna_parallel_decompress_fasta,
    rna_write_compressed_fastq, rna_write_compressed_fasta,
    rna_read_compressed_file
)

from .protein import (
    protein_parallel_compress, protein_parallel_decompress,
    protein_write_compressed, protein_read_compressed,
    read_fasta as protein_read_fasta
)

from .file import file_compress, file_decompress

from .ncd import NCDClassifier
from .lzjd import LZJDKNNCorrector, read_fasta_folder
from .wbncd import BCDClassifier

from .tools import process_all, finally_cleanup

def perform_distance_analysis(args):
    output_dir = args.output or f'biox_{args.method}_output'
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        if args.method == 'lzjd':
            sequences = read_fasta_folder(args.input)
            if not sequences:
                raise ValueError("No valid sequence files found")
            
            corrector = LZJDKNNCorrector(
                sequences=sequences,
                taxonomy_level=args.taxonomy_level,
                confidence_threshold=args.confidence
            )
            
            predictions, true_labels = corrector.perform_knn_classification(k=args.neighbors)
            corrected_distances = corrector.correct_distances(base_alpha=args.alpha)

            macro_f1 = f1_score(true_labels, predictions, average='macro')
            
            corrector.save_results(
                predictions, true_labels, corrector.probability_matrix,
                macro_f1, output_dir, tree_method=args.tree
            )
        else:
            if not process_all(args.input, args.input, args.type):
                raise ValueError("Failed to process input sequences")
            
            classifier_class = NCDClassifier if args.method == 'ncd' else BCDClassifier
            classifier = classifier_class(
                compressed_dir=args.input,
                taxonomy_level=args.taxonomy_level
            )
            
            sequence_names, labels = classifier.load_compressed_sizes()
            if not sequence_names:
                raise ValueError("No valid sequence files found")
            
            classifier.labels = labels
            classifier.build_distance_matrix(sequence_names, n_jobs=args.jobs)
            
            predictions, true_labels, probabilities, macro_f1 = \
                classifier.perform_knn_classification(k=args.neighbors)
            
            classifier.confidence_threshold = args.confidence
            classifier.correct_distances(base_alpha=args.alpha)
            
            classifier.save_results(
                predictions, true_labels, probabilities, macro_f1,
                output_dir, tree_method=args.tree
            )
            
            finally_cleanup(args.input)
        
        print(f"Analysis completed. Results saved to {output_dir}")
        
    except Exception as e:
        print(f"Error during {args.method.upper()} analysis: {str(e)}")
        raise

def main():
    parser = argparse.ArgumentParser(description='BioX: A tool for biological sequence compression and analysis')
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('-c', '--compress', action='store_true',
                           help='Compress mode')
    mode_group.add_argument('-d', '--decompress', action='store_true',
                           help='Decompress mode')
    mode_group.add_argument('-a', '--analysis', action='store_true',
                           help='Sequence analysis mode')
    parser.add_argument('-i', '--input',
                       required=True,
                       help='Input file/directory path')
    parser.add_argument('-o', '--output',
                       help='Output file/directory path')
    
    comp_group = parser.add_argument_group('Compression/Decompression options')
    comp_group.add_argument('-t', '--type', 
                           choices=['dna', 'rna', 'protein', 'file'],
                           help='Sequence type (dna/rna/protein) or regular file')
    comp_group.add_argument('-l', '--level', 
                           type=int, 
                           choices=range(1, 10),
                           default=3,
                           help='Compression level (1-9, default: 3)')
    comp_group.add_argument('-ps', '--plus_line',
                           choices=['ignore', 'compress'],
                           default='ignore',
                           help='FASTQ plus line handling')
    comp_group.add_argument('--num_processes',
                           type=int,
                           default=None,
                           help='Number of parallel processes')
    comp_group.add_argument('-p', '--plant',
                           action='store_true',
                           help='Use plant genome compression scheme')
    comp_group.add_argument('-s', '--split', type=int, choices=range(2, 11),
                       help='Split output into N volumes (2-10)')
    
    analysis_group = parser.add_argument_group('Sequence Analysis options')
    analysis_group.add_argument('--method', '-m',
                              choices=['ncd', 'wbncd', 'lzjd'],
                              default='ncd',
                              help='Distance calculation method')
    analysis_group.add_argument('--tax', '--taxonomy-level',
                              dest='taxonomy_level',
                              choices=['kingdom', 'phylum', 'class', 'order', 
                                      'family', 'genus', 'species'],
                              default='class',
                              help='NCBI taxonomy level for classification')
    analysis_group.add_argument('-k', '--neighbors',
                              type=int,
                              default=1,
                              help='Number of neighbors for KNN classification')
    analysis_group.add_argument('--alpha',
                              type=float,
                              default=0.3,
                              help='Distance correction coefficient (0-1)')
    analysis_group.add_argument('--confidence',
                              type=float,
                              default=0.6,
                              help='Confidence threshold for classification')
    analysis_group.add_argument('-j', '--jobs',
                              type=int,
                              default=os.cpu_count(),
                              help='Number of parallel jobs')
    analysis_group.add_argument('--tree',
                              choices=['single', 'average', 'weighted', 'complete'],
                              default='single',
                              help='Method for phylogenetic tree construction')
    
    args = parser.parse_args()

    if args.analysis:
        if not args.type and args.method != 'lzjd':
            args.type = 'dna'
            print("No sequence type specified, using default type: DNA")
        perform_distance_analysis(args)
        return

    if not args.output:
        if args.decompress:
            if args.input.endswith('.biox'):
                args.output = args.input[:-5] + '.decoded'
        else:
            args.output = args.input + '.biox'

    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' does not exist")
        return

    if not args.compress and not args.decompress:
        parser.error("Must specify either -c (compress) or -d (decompress) mode")

    try:
        if args.type == 'file':
            if args.compress:
                print(f"\nCompressing file: {args.input}")
                original_size = os.path.getsize(args.input)
                print(f"Original file size: {original_size / 1024 / 1024:.2f} MB")
                
                start_time = time.time()
                file_compress(args.input, args.output, preset=args.level)
                end_time = time.time()
                compress_time = end_time - start_time
                
                compressed_size = os.path.getsize(args.output)
                compression_ratio = compressed_size / original_size
                
                print(f"\nCompression completed:")
                print(f"- Time taken: {compress_time:.2f} seconds")
                print(f"- Compressed size: {compressed_size / 1024 / 1024:.2f} MB")
                print(f"- Compression ratio: {compression_ratio:.2f}x")
                
            else:  
                print(f"\nDecompressing file: {args.input}")
                start_time = time.time()
                
                file_decompress(args.input, args.output)
                
                end_time = time.time()
                print(f"\nDecompression completed! Time taken: {end_time - start_time:.2f} seconds")
            
            return

        if args.compress:
            file_format = detect_file_format(args.input)
            print(f"Detected file format: {file_format}")

            original_size = os.path.getsize(args.input)
            print(f"Original file size: {original_size / 1024 / 1024:.2f} MB")

            if args.type == 'protein':
                if file_format != 'fasta':
                    print("Error: Protein sequences only support FASTA format")
                    return
                
                print(f"\nReading FASTA file: {args.input}")
                sequences = protein_read_fasta(args.input)
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
                    print(f"Reading FASTQ file: {args.input}")
                    sequences = read_fastq(args.input)
                else:
                    print(f"Reading FASTA file: {args.input}")
                    sequences = read_fasta(args.input)

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
            print(f"\nReading compressed file: {args.input}")
            start_time = time.time()

            if args.type == 'protein':
                compressed_data = protein_read_compressed(args.input)
                decompressed_data = protein_parallel_decompress(compressed_data, processes=args.num_processes)
                file_format = 'fasta'
            else:
                read_func = dna_read_compressed_file if args.type == 'dna' else rna_read_compressed_file
                file_format, compressed_data = read_func(args.input)
                
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
