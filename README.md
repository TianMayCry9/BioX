![image](https://github.com/TianMayCry9/BioX/blob/master/BioX.png)
# BioX - Biological Sequence Compression Tool

BioX is an efficient, lossless compression tool designed specifically for biological sequence data. It utilizes innovative compression algorithms to achieve high compression ratios and fast decompression speeds. It supports compression of DNA, RNA, protein sequences, and FASTQ data, and can handle sequences with or without quality scores. Fully supports FASTA and FASTQ formats, including all special characters, and there are no limitations on sequence length or quantity.

## Features

- Supports multiple biological sequence formats (FASTA, FASTQ)
- Supports compression of DNA, RNA, protein sequences
- Special optimizations for plant genomes
- Multi-process parallel processing
- Automatic sequence type detection
- Supports large file chunk processing
- Cross-platform support (Windows, Linux)

## Installation

### Install via conda (Recommended)

```bash
conda install -c bioconda biox
```

### Build from source

Prerequisites: Python 3.6+, pip, numpy, tqdm

```bash
# Clone the repository
git clone https://github.com/TianxXiao/biox.git

# Navigate into the directory and install
cd biox
pip install .
```

## Usage

### Compress

```bash
# Basic usage
biox -c -t dna input.fasta

# Specify compression level (1-9, default is 5)
biox -c -t dna -l 9 input.fasta

# Use multiprocessing for acceleration
biox -c -t dna -n 4 input.fasta

# Process plant genomes
biox -c -t dna -ps compress input.fasta

# Specify output file
biox -c -t dna -o output.biox input.fasta
```

### Decompress

```bash
# Basic usage
biox -d -t dna input.biox

# Specify output file
biox -d -t dna -o output.fasta input.biox
```

### Parameter Overview
To see the possible options type：
```bash
biox -h
```

This will print the following options:
```bash
usage: main.py [-h] [-c] [-d] -t {dna,rna,protein,file} [-l {1,2,3,4,5,6,7,8,9}] [-ps {ignore,compress}] [-n NUM_PROCESSES] [-p] [-o OUTPUT] input_file

Biological Sequence Compression Tool

positional arguments:
  input_file            Input file path (FASTQ/FASTA)

options:
  -h, --help            show this help message and exit
  -c, --compress        Compression mode
  -d, --decompress      Decompression mode
  -t {dna,rna,protein,file}, --type {dna,rna,protein,file}
                        Sequence type (dna/rna/protein) or regular file
  -l {1,2,3,4,5,6,7,8,9}, --level {1,2,3,4,5,6,7,8,9}
                        Compression level (1-9, default: 5)
  -ps {ignore,compress}, --plus_line {ignore,compress}
                        FASTQ plus line handling
  -n NUM_PROCESSES, --num_processes NUM_PROCESSES
                        Number of parallel processes
  -p, --plant           Use plant genome compression scheme
  -o OUTPUT, --output OUTPUT
                        Output file path (default: input_file.biox)
```

## FAQ

1.**How to choose the right compression level?**
   - Level 1: Fastest compression speed, lower compression ratio
   - Level 5: Balanced compression speed and ratio (default)
   - Level 9: Highest compression ratio, slower speed

2.**How to handle very large files?**
   - It is recommended to use the plant genome compression parameter (-p) and adjust the compression level according to available memory.

3.**How to detect the file type?**
   - BioX automatically detects the file type, but you may need to manually specify the sequence type with the -t parameter.

4.**Will the "+" line in FASTQ be ignored?**
   - BioX automatically ignores the "+" line in FASTQ files, but you can manually retain it using the -ps parameter.

## Citation

If you use BioX in your research, please cite:
```
[Paper information to be added]
```

## License

This project is licensed under the GPL-3.0 License.See the  [LICENSE](LICENSE) for details. 

## Contribution Guidelines

We welcome issues and suggestions! For any issue, let us know at [issues link](https://github.com/TianMayCry9/BioX/issues).
