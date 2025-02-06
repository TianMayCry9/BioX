![image](https://github.com/TianMayCry9/BioX/blob/master/BioX.png)
# BioX - Biological Sequence Compression and Analysis Tool

BioX is an efficient, lossless compression tool designed specifically for biological sequence data. It utilizes innovative compression algorithms to achieve high compression ratios and fast decompression speeds. It can handle all file types, especially DNA, RNA, and protein FASTA and FASTQ sequence files, and can handle sequences with or without quality scores. Fully supports FASTA and FASTQ formats, including all special characters, and there are no limitations on sequence length or quantity.

In addition to compression, BioX offers sequence analysis features, including distance calculation methods (ncd, bcd, lzjd), KNN classification with customizable neighbor settings, taxonomy-based classification at various NCBI levels, distance correction coefficients, and phylogenetic tree construction.

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

### Build from Source

**Prerequisites:**
- Python 3.6 or higher
- pip
- numpy
- tqdm

**Installation Steps:**

1. Clone the repository:
```bash
git clone https://github.com/TianMayCry9/BioX.git
cd biox
```
2. Choose one of the following installation methods:
- Option 1 (Recommended):
```bash
pip install .
```
- Option 2 (Linux/Mac):
```bash
chmod +x install.sh
./install.sh
```
- Option 3 (Windows):
```bash
install.bat
```

### Download the Pre-built .exe File (Windows)
Simply download the .exe application from our [GitHub Releases](https://github.com/TianMayCry9/BioX/releases/tag/v1.1.0), double-click to launch, and start compressing your files with ease!

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

# Volume Compression
biox -c -s 4 input.fasta

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

### Sequence Analysis

```bash
# Basic usage for analysis
biox -a -i /tree_test/

# Specify distance calculation method (ncd, bcd, lzjd)
biox -a -i /tree_test/ --method ncd

# Specify NCBI taxonomy level
biox -a -i /tree_test/ --tax class

# KNN classification with 1 neighbors
biox -a -i /tree_test/ -k 1

# Set distance correction coefficient
biox -a -i /tree_test/ --alpha 0.3

# Set confidence threshold for KNN classification
biox -a -i /tree_test/ --confidence 0.95

# Use parallel jobs for analysis
biox -a -i /tree_test/ -j 4

# Construct phylogenetic tree
biox -a -i /tree_test/ --tree single
```

### Parameter Overview
To see the possible options type：
```bash
biox -h
```

This will print the following options:
```bash
usage: biox [-h] (-c | -d | -a) -i INPUT [-o OUTPUT] [-t {dna,rna,protein,file}] [-l {1,2,3,4,5,6,7,8,9}] [-ps {ignore,compress}] [--num_processes NUM_PROCESSES] [-p]
               [--method {ncd,bcd,lzjd}] [--tax {kingdom,phylum,class,order,family,genus,species}] [-k NEIGHBORS] [--alpha ALPHA] [--confidence CONFIDENCE] [-j JOBS]
               [--tree {single,average,weighted,complete}]

BioX: A tool for biological sequence compression and analysis

options:
  -h, --help            show this help message and exit
  -c, --compress        Compress mode
  -d, --decompress      Decompress mode
  -a, --analysis        Sequence analysis mode
  -i INPUT, --input INPUT
                        Input file/directory path
  -o OUTPUT, --output OUTPUT
                        Output file/directory path

Compression/Decompression options:
  -t {dna,rna,protein,file}, --type {dna,rna,protein,file}
                        Sequence type (dna/rna/protein) or regular file
  -l {1,2,3,4,5,6,7,8,9}, --level {1,2,3,4,5,6,7,8,9}
                        Compression level (1-9, default: 3)
  -ps {ignore,compress}, --plus_line {ignore,compress}
                        FASTQ plus line handling
  --num_processes NUM_PROCESSES
                        Number of parallel processes
  -p, --plant           Use plant genome compression scheme

Sequence Analysis options:
  --method {ncd,bcd,lzjd}, -m {ncd,bcd,lzjd}
                        Distance calculation method
  --tax {kingdom,phylum,class,order,family,genus,species}, --taxonomy-level {kingdom,phylum,class,order,family,genus,species}
                        NCBI taxonomy level for classification
  -k NEIGHBORS, --neighbors NEIGHBORS
                        Number of neighbors for KNN classification
  --alpha ALPHA         Distance correction coefficient (0-1)
  --confidence CONFIDENCE
                        Confidence threshold for classification
  -j JOBS, --jobs JOBS  Number of parallel jobs
  --tree {single,average,weighted,complete}
                        Method for phylogenetic tree construction
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

5.**When should I use the KNN correction feature?**
   - The KNN correction feature should be used when there is a relatively large sample size of the biological data, as the classification is based on the NCBI taxonomy standards.

## Citation

If you use BioX in your research, please cite:
```
[Paper information to be added]
```

## License

This project is licensed under the GPL-3.0 License.See the  [LICENSE](LICENSE) for details. 

## Contribution Guidelines

We welcome issues and suggestions! For any issue, let us know at [issues link](https://github.com/TianMayCry9/BioX/issues).
