from typing import Union
from pathlib import Path
from rich.logging import logging
import biopython
import polars as pl
import typer

log = logging.get_name

def load_fasta(filepath: Path):
    return open(filepath)

def load_bowtie_file(filepath: Path):
    # Load the file
    df = pl.scan_csv(filepath, has_header=False, new_columns=["single"])
    # Get the even rows
    df_p = df.filter(pl.arange(0, df.height) % 2 == 0).select(pl.col("single").str.split("\t").list.get(1).alias("chr_name"), pl.col("single").str.split("\t").list.get(2).alias("coord_p"), pl.col("single").str.split("\t").list.get(-1).alias("sequence_p")).with_columns(pl.col("coord_p").cast(pl.Int64))
    # Now it the turn for the odd ones
    df_m = df.filter(pl.arange(0, df.height) % 2 != 0).select(pl.col("single").str.split("\t").list.get(1).alias("chr_name_m"), pl.col("single").str.split("\t").list.get(2).alias("coord_m"), pl.col("single").str.split("\t").list.get(-1).alias("sequence_m")).with_columns(pl.col("coord_m").cast(pl.Int64))

    # Concatenate horizontally both
    return pl.concat([df_p, df_m], how="horizontal").drop(pl.col("chr_name_m"))


def get_intervals(alignments: Union[pl.DataFrame, pl.LazyFrame]):
    intervals = alignments.select(start=pl.col("coord_p"), end=pl.col("coord_m") + pl.col("sequence_m").str.len_chars() - 1, center= (pl.col("coord_p") + pl.col("coord_m") + pl.col("sequence_m").str.len_chars() -1 ) // 2)

    counts = 

def main():
    log.info("Hello from nucwave2!")
    log.info("Loading genome from FASTA file...")

    log.info("")

if __name__ == "__main__":
    main()
