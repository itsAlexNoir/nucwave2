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

    # intervals = alignments.select(
    #     chr_name=pl.col("chr_name"),
    #     start=pl.col("coord_p"),
    #     end=pl.col("coord_m") + pl.col("sequence_m").str.len_chars() - 1,
    #     center=(
    #         (
    #             pl.col("coord_p")
    #             + (pl.col("coord_m") + pl.col("sequence_m").str.len_chars() - 1)
    #         )
    #         * 0.5
    #     ).cast(pl.Int64),
    # )
    # intervals = intervals.with_columns(d2d_long=pl.col("end") - pl.col("start") + 1)

    intervals = alignments.select(
        chr_name=pl.col("chr_name"),
        start=pl.col("coord_p"),
        end=pl.col("coord_m") + pl.col("sequence_m").str.len_chars() - 1,
    ).with_columns(d2d_long=pl.col("end") - pl.col("start") + 1)

    intervals = intervals.filter(
        (pl.col("d2d_long") > minsize) & (pl.col("d2d_long") < maxsize)
    )
    intervals = intervals.with_columns(
        start_ext=pl.col("start") - extension,
        end_ext=pl.col("end") + extension,
        center=((pl.col("start") + pl.col("end")) * 0.5).cast(pl.Int32),
    )
    available_chromosomes = [
        row[0] for row in intervals.select("chr_name").unique().collect().rows()
    ]
    histograms = {}
    for chr_name in available_chromosomes:
        log.info(f"Processing chromosome: {chr_name}")
        len_seq = len(genome[chr_name])
        subset = intervals.filter(pl.col("chr_name") == chr_name)
        histograms[chr_name] = subset.select(
            starts=pl.col("start").hist(bin_count=len_seq),
            ends=pl.col("end").hist(bin_count=len_seq),
            centers=pl.col("center").hist(bin_count=len_seq),
            starts_ext=pl.col("start_ext").hist(bin_count=len_seq),
            ends_ext=pl.col("end_ext").hist(bin_count=len_seq),
        )

    sizereads = intervals.select(pl.col("d2d_long").hist(bin_count=maxsize))
    num_filtered_samples = intervals.select(pl.len()).collect().item()

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
    wig_files: bool = False,
    minsize: int = 40,
    maxsize: int = 200,
):
    """
    Main entry point for nucwave2.

    Loads a reference genome and alignments, computes fragment intervals and
    histograms, logs timing information, and optionally writes intermediate
    WIG/CSV files to an output directory.

    Parameters
    ----------
    alignment_file : pathlib.Path
        Path to a Bowtie-format alignment file. Must exist or a FileNotFoundError
        will be raised.
    genome_file : pathlib.Path
        Path to a FASTA-format genome file. Must exist or a FileNotFoundError will
        be raised.
    output_dir : pathlib.Path, optional
        Directory where results and intermediate files are written. Created if it
        does not exist. Default: Path("results").
    wig_files : bool, optional
        When True, write intermediate wiggle-format files and a read-size
        histogram CSV-like file to output_dir. Default: False.
    minsize : int, optional
        Minimum fragment size (inclusive) used to filter alignments when computing
        intervals. Default: 40.
    maxsize : int, optional
        Maximum fragment size (inclusive) used to filter alignments when computing
        intervals. Default: 200.

    Behavior
    --------
    1. Validates existence of genome_file and alignment_file, creates output_dir.
    2. Loads the genome using load_fasta(genome_file) and logs elapsed time.
    3. Loads alignments using load_bowtie_file(alignment_file) and logs elapsed time.
    4. Calls get_intervals(alignments, minsize=minsize, maxsize=maxsize, genome=genome)
       to compute:
         - num_filtered_samples: number of samples passing filters
         - sizereads: an object supporting with_row_index().write_csv(...)
         - histograms: a mapping (chromosome id -> table) where each table exposes
           columns used in .select("starts"), .select("ends"), .select("centers")
           and supports .collect().iter_rows()
       and logs elapsed time and the number of filtered samples.
    5. If wig_files is True, writes the following files to output_dir and logs times:
         - readsize_histogram.wig
           Written by calling sizereads.with_row_index().write_csv(filename, include_header=False)
         - cut_p.wig and cut_m.wig
           For each chromosome, writes a simple WIG fixedStep track using the
           per-chromosome histogram table's "starts" and "ends" selections.
         - depth_trimmed_PE.wig
           Writes counts from the "centers" selection for each chromosome.
         - depth_complete_PE.wig
           Writes cumulative depth computed from (starts - ends) rows for each chromosome.
       Note: the function writes simple track headers of the form:
           track type=wiggle_0 name=<filename> description=<filename.stem>
       and uses fixedStep chrom=<chrid> start=1 step=1 for each chromosome block.

    Returns
    -------
    None

    Raises
    ------
    FileNotFoundError
        If genome_file or alignment_file does not exist.

    Notes
    -----
    - The function relies on helper functions: load_fasta, load_bowtie_file, and
      get_intervals. Their behaviors and return types are assumed as described
      above.
    - The code logs progress and timing information for each major step.
    - The histograms objects are expected to be queryable in the manner used in the
      implementation (select/collect/iter_rows) — e.g., a Polars DataFrame or a
      compatible abstraction.
    - Output WIG files produced here are simple fixedStep tracks suitable for
      downstream visualization or processing by tools that accept plain WIG format.
    """
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

    log.info("Getting intervals from alignments...")
    time0 = time.time()
    num_filtered_samples, sizereads, histograms = get_intervals(
        alignments, minsize=minsize, maxsize=maxsize, genome=genome
    )
    time1 = time.time()
    log.info(f"Number of filtered samples: {num_filtered_samples}")
    log.info(f"\tTime taken to get intervals: {time1 - time0:.2f} seconds")

    log.info("#" * 45)
    log.info("#" * 45)

    if wig_files:

        log.info("Collecting all the data from memory for writing files...")
        # Execute histogram dataframe only once, collecting all data needed for writing files
        time0 = int(time.time())
        histograms = {chr_name: hist.collect() for chr_name, hist in histograms.items()}
        time1 = int(time.time())
        log.info("\tTime taken: %d seconds" % (time1 - time0))

        log.info("#" * 45)

        log.info("Writing intermediate files: fragment size histogram")
        time0 = int(time.time())
        filename = output_dir / "readsize_histogram.wig"
        sizereads.collect().with_row_index().write_csv(filename, include_header=False)
        log.info("\tTime taken: %d seconds" % (time.time() - time0))

        log.info("#" * 45)
        log.info("Writing intermediate files: cut points per strand")
        time0 = int(time.time())
        filename = output_dir / "cut_p.wig"
        with open(filename, "w") as WIG:
            for chrid, hist in histograms.items():
                WIG.write(
                    f"track type=wiggle_0 name={filename} description={filename.stem}\n"
                )
                WIG.write(f"fixedStep chrom={chrid} start=1 step=1\n")
                for row in hist.select("starts").iter_rows():
                    WIG.write(f"{row[0]}\n")

        filename = output_dir / "cut_m.wig"
        with open(filename, "w") as WIG:
            for chrid, hist in histograms.items():
                WIG.write(
                    f"track type=wiggle_0 name={filename} description={filename.stem}\n"
                )
                WIG.write(f"fixedStep chrom={chrid} start=1 step=1\n")
                for row in hist.select("ends").iter_rows():
                    WIG.write(f"{row[0]}\n")
        time1 = int(time.time())
        log.info("\tTime taken: %d seconds" % (time1 - time0))

        log.info("#" * 45)

        log.info("Writing intermediate files: PE center count")
        time0 = int(time.time())
        filename = output_dir / "depth_trimmed_PE.wig"
        with open(filename, "w") as WIG:
            for chrid, hist in histograms.items():
                WIG.write(
                    f"track type=wiggle_0 name={filename} description={filename.stem}\n"
                )
                WIG.write(f"fixedStep chrom={chrid} start=1 step=1\n")
                for row in hist.select("centers").iter_rows():
                    WIG.write(f"{row[0]}\n")

        time1 = int(time.time())
        log.info("\tTime taken: %d seconds" % (time1 - time0))
        log.info("#" * 45)

        log.info("Writing intermediate files: depth coverage for complete PE reads")
        time0 = int(time.time())
        filename = output_dir / "depth_complete_PE.wig"
        with open(filename, "w") as WIG:
            for chrid, hist in histograms.items():
                WIG.write(
                    f'track type=wiggle_0 name={filename} description="{filename.stem}"\n'
                )
                WIG.write(f"fixedStep chrom={chrid} start=1 step=1\n")
                deep = 0
                for row in hist.select(pl.col("starts") - pl.col("ends")).iter_rows():
                    deep += row[0]
                    WIG.write(f"{deep}\n")
        time1 = int(time.time())
        log.info("\tTime taken: %d seconds" % (time1 - time0))

        log.info("#" * 45)

        log.info("#" * 45)

        log.info("Calculation completed successfully.")


if __name__ == "__main__":
    cli()
