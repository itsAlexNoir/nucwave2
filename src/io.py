from pathlib import Path
import polars as pl
from Bio import SeqIO


def load_fasta(filepath: Path):
    """
    Load sequences from a FASTA file and return a mapping of record IDs to sequence strings.

    Parameters
    ----------
    filepath : pathlib.Path or str
        Path to the FASTA file to read.

    Returns
    -------
    dict[str, str]
        Dictionary mapping each FASTA record ID to its sequence as a plain string.
        If the file contains no records, an empty dictionary is returned.
        If multiple records share the same ID, later records will overwrite earlier ones.

    Raises
    ------
    FileNotFoundError
        If the supplied filepath does not exist or cannot be opened.
    IOError
        For other I/O related errors when reading the file.
    Exception
        Any parsing errors raised by Bio.SeqIO may be propagated.

    Notes
    -----
    - Parsing is performed using Bio.SeqIO.parse(..., "fasta") from Biopython.
    - The function opens the file for reading; callers may prefer using a context manager
      if they need explicit control over file handling.
    """
    fasta_seqs = SeqIO.parse(open(filepath), "fasta")
    return {fasta.id: str(fasta.seq) for fasta in fasta_seqs}


def load_bowtie_file(filepath: Path):
    """
    Load a paired Bowtie-style file into a Polars DataFrame by combining even/odd line pairs.

    Each record in the input file is expected to occupy two consecutive lines:
    - the even (0-based) line contains the "plus" read information,
    - the odd line immediately after contains the corresponding "minus" read information.

    Each line is treated as a single tab-separated string and the function extracts:
    - field at index 1 -> chromosome name
    - field at index 2 -> coordinate (cast to Int64)
    - last field -> sequence

    Parameters
    ----------
    filepath : Path
        Path to the input file. The file is read using polars.scan_csv(..., has_header=False)
        into a single-column table and then parsed as described above.

    Returns
    -------
    pl.DataFrame
        A DataFrame with one row per paired record and columns:
        - "chr_name"   : chromosome name from the even (plus) lines
        - "coord_p"    : coordinate from the plus lines (Int64)
        - "sequence_p" : sequence from the plus lines
        - "coord_m"    : coordinate from the minus lines (Int64)
        - "sequence_m" : sequence from the minus lines

    Notes
    -----
    - The implementation relies on 0-based indexing to split even (plus) and odd (minus)
      lines and concatenates the parsed results horizontally.
    - A temporary "chr_name_m" column produced from the minus lines is dropped because
      the chromosome name is expected to be the same for both lines of a pair.

    Raises
    ------
    FileNotFoundError
        If the provided filepath does not exist or cannot be opened.
    ValueError
        If a coordinate field cannot be cast to an integer or if the input contains
        an unexpected number/structure of fields (e.g., unpaired lines).

    Examples
    --------
    >>> from pathlib import Path
    >>> df = load_bowtie_file(Path("reads.bowtie"))
    >>> list(df.columns)
    ['chr_name', 'coord_p', 'sequence_p', 'coord_m', 'sequence_m']
    """
    # Load the file
    df = pl.scan_csv(filepath, has_header=False, new_columns=["single"])
    # Get the even rows
    # row_count = df.height
    row_count = df.select(pl.len()).collect().item()

    df_p = (
        df.filter(pl.arange(0, row_count) % 2 == 0)
        .select(
            pl.col("single").str.split("\t").list.get(1).alias("chr_name"),
            pl.col("single").str.split("\t").list.get(2).alias("coord_p"),
            pl.col("single").str.split("\t").list.get(-1).alias("sequence_p"),
        )
        .with_columns(pl.col("coord_p").cast(pl.Int64))
    )
    # Now it the turn for the odd ones
    df_m = (
        df.filter(pl.arange(0, row_count) % 2 != 0)
        .select(
            pl.col("single").str.split("\t").list.get(1).alias("chr_name_m"),
            pl.col("single").str.split("\t").list.get(2).alias("coord_m"),
            pl.col("single").str.split("\t").list.get(-1).alias("sequence_m"),
        )
        .with_columns(pl.col("coord_m").cast(pl.Int64))
    )

    # Concatenate horizontally both
    return pl.concat([df_p, df_m], how="horizontal").drop(pl.col("chr_name_m"))
