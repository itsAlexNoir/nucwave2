# nucwave2

A genomics analysis tool for processing nucleosome positioning data from paired-end sequencing alignments. nucwave2 analyzes Bowtie alignment files and genome FASTA files to generate fragment size histograms and genomic position distributions.

<!-- ## Features

- Processes Bowtie alignment files with paired-end reads
- Filters fragments by size (configurable min/max thresholds)
- Generates histograms for fragment starts, ends, and centers
- Outputs results in WIG format for genome browser visualization
- Fast data processing using Polars -->

## Installation

### Using uv (recommended)

It is recomended to use uv to install the package. If uv is not installed on your machine, you can do so by:

```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Then, you have to clone the repo

```bash
# Clone the repository
git clone https://github.com/itsAlexNoir/nucwave2
cd nucwave2

# Install dependencies
uv sync
```

### Using pip

You can also do it the classic (pip) way:

```bash
# Clone the repository
git clone https://github.com/itsAlexNoir/nucwave2
cd nucwave2

# Create a virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -e .
```

## Usage

Run the main script with:

```bash
# Using uv
uv run python main.py <alignment_file> <genome_file> [OPTIONS]

# Using standard Python
python main.py <alignment_file> <genome_file> [OPTIONS]
```

### Arguments

- `alignment_file`: Path to the Bowtie alignment file
- `genome_file`: Path to the genome FASTA file

### Options

- `--output-dir PATH`: Output directory for results (default: `results/`)
- `--write-intermediate-files`: Write intermediate output files
- `--minsize INT`: Minimum fragment size to filter (default: 40)
- `--maxsize INT`: Maximum fragment size to filter (default: 200)

### Example

```bash
uv run python main.py data/MN_opaque1.bowtie data/C_albicans_WO-1_chromosomes.fasta --output-dir results --minsize 40 --maxsize 200
```

## Output

The tool generates:

- Fragment size histogram in WIG format (`readsize_histogram.wig`)
- Genomic position histograms for each chromosome
- Log files in the `log/` directory

## Requirements

- Python >= 3.13
- biopython >= 1.86
- numpy >= 2.4.1
- polars >= 1.37.1
- typer >= 0.21.1
- rich >= 14.2.0
