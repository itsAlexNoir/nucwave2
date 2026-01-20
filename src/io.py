from pathlib import Path
import polars as pl
from Bio import SeqIO


def load_fasta(filepath: Path):
    fasta_seqs = SeqIO.parse(open(filepath), "fasta")
    return {fasta.id: str(fasta.seq) for fasta in fasta_seqs}


def load_bowtie_file(filepath: Path):
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
