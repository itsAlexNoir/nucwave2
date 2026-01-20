import logging
import time
from typing import Union, Dict
from pathlib import Path
import polars as pl
import typer

from src.utils import set_logging
from src.io import load_fasta, load_bowtie_file

set_logging("log/nucwave2.log")
log = logging.getLogger(__name__)


def get_intervals(
    alignments: Union[pl.DataFrame, pl.LazyFrame],
    genome: Dict[str, str],
    minsize: int = 40,
    maxsize: int = 200,
):
    # Define the extension
    extension = int(minsize * 0.5)

    intervals = alignments.select(
        chr_name=pl.col("chr_name"),
        start=pl.col("coord_p"),
        end=pl.col("coord_m") + pl.col("sequence_m").str.len_chars() - 1,
        center=(
            (
                pl.col("coord_p")
                + (pl.col("coord_m") + pl.col("sequence_m").str.len_chars() - 1)
            )
            * 0.5
        ).cast(pl.Int64),
    )
    intervals = intervals.with_columns(d2d_long=pl.col("end") - pl.col("start") + 1)

    intervals = intervals.filter(pl.col("d2d_long") >= minsize).filter(
        pl.col("d2d_long") <= maxsize
    )
    intervals = intervals.with_columns(
        start_ext=pl.col("start") - extension, end_ext=pl.col("end") + extension
    )
    num_filtered_samples = intervals.select(pl.len()).collect().item()

    histograms = {}
    for chr_name, seq in genome.items():
        subset = intervals.filter(pl.col("chr_name") == chr_name)
        if not subset.collect().is_empty():
            histograms[chr_name] = subset.select(
                starts=pl.col("start").hist(bin_count=len(seq)),
                ends=pl.col("end").hist(bin_count=len(seq)),
                centers=pl.col("center").hist(bin_count=len(seq)),
                starts_ext=pl.col("start_ext").hist(bin_count=len(seq)),
                ends_ext=pl.col("end_ext").hist(bin_count=len(seq)),
            )

    sizereads = intervals.select(pl.col("d2d_long").hist(bin_count=maxsize)).collect()
    print("Here")

    return num_filtered_samples, sizereads, histograms


def cli():
    """
    Entry point that exposes the program's CLI.

    This function delegates execution to typer.run(main), turning the
    callable `main` into a Typer/Click-based command-line interface and
    running it with the current process arguments (sys.argv). It accepts
    no parameters and does not return a value; its primary effect is to
    start and manage the CLI lifecycle, including argument parsing,
    dispatching to the `main` function, and propagating CLI-related
    exceptions or exit codes.

    Notes:
    - The referenced `main` callable must be defined in the same module
        and be compatible with Typer (i.e., it may accept typed function
        parameters that map to CLI options/arguments).
    - Typical usage is to call this function from a module guard:
            if __name__ == "__main__":
                    cli()
    - Any SystemExit or exceptions raised by the CLI or `main` will
        behave as they do under Typer/Click (may terminate the process).

    No return value.
    """
    typer.run(main)


def main(
    alignment_file: Path,
    genome_file: Path,
    output_dir: Path = Path("results"),
    write_intermediate_files: bool = False,
    minsize: int = 40,
    maxsize: int = 200,
):
    log.info("=" * 65)
    log.info(" " * 25 + "Hello from nucwave2!")
    log.info("=" * 65)

    if not genome_file.exists():
        raise FileNotFoundError(f"Genome file {genome_file} does not exist.")
    if not alignment_file.exists():
        raise FileNotFoundError(f"Alignment file {alignment_file} does not exist.")
    output_dir.mkdir(parents=True, exist_ok=True)

    log.info("Loading genome from FASTA file...")
    time0 = time.time()
    genome = load_fasta(genome_file)
    time1 = time.time()
    log.info(f"\tTime taken to load genome: {time1 - time0:.2f} seconds")
    log.info("Genome loaded successfully.")

    log.info("Loading alignments from Bowtie file...")
    time0 = time.time()
    alignments = load_bowtie_file(alignment_file)
    time1 = time.time()
    log.info(f"\tTime taken to load alignments: {time1 - time0:.2f} seconds")
    log.info("Alignments loaded successfully.")

    time0 = time.time()
    num_filtered_samples, sizereads, histograms = get_intervals(
        alignments, minsize=minsize, maxsize=maxsize, genome=genome
    )
    time1 = time.time()
    log.info(f"Number of filtered samples: {num_filtered_samples}")
    log.info(f"\tTime taken to get intervals: {time1 - time0:.2f} seconds")

    log.info("#" * 45)
    log.info("#" * 45)

    if write_intermediate_files:
        time0 = int(time.time())
    log.info("Writing intermediate files: fragment size histogram")
    filename = output_dir / "readsize_histogram.wig"
    sizereads.with_row_index().write_csv(filename, include_header=False)
    log.info("  Time taken: %d seconds" % (time.time() - time0))

    log.info("#" * 45)

    # suffix = 'cut_p'
    # WIG = open(WIG_file + "_%s.wig" % suffix,"w")
    # for chrid in chrseq.keys() :
    #     WIG.write('track type=wiggle_0 name=%s_%s description="%s_%s"\n' % (exper,suffix,exper,suffix))
    #     WIG.write('fixedStep chrom=%s start=1 step=1\n' % chrid)
    # for i in range(chrlen[chrid]) :
    #   WIG.write('%d\n' % iniciosPE[chrid][i] )
    # WIG.close()


if __name__ == "__main__":
    cli()
