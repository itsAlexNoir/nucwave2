import logging
from rich import print as rprint
from rich.logging import RichHandler
from pathlib import Path
import numpy as np
import pandas as pd
import typer

logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler()]
)
log = logging.getLogger(__name__)


def render_welcome_logo() -> None:
    """Render a nucleosome-style banner using Rich.

    Returns:
        None
    """
    logo = (
        "\n\n"
        "[blue]     /\\         /\\           /\\             /\\        /\\     /\\ [/blue]\n"
        "[blue]    /  \\       /  \\         /  \\    /\\     /  \\      /  \\   /  \\ [/blue]\n"
        "[blue]___/    \\_____/    \\_______/    \\__/  \\___/    \\____/    \\_/    \\___[/blue]\n"
        "[blue]-----------------------------------------------------------------------[/blue]\n"
        "[sky_blue3]~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~[/sky_blue3]\n"
        "[gold3](o)   (o)   (o)   (o)   (o)   (o)   (o)   (o)   (o)   (o)[/gold3]\n"
        "[sky_blue3]~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~[/sky_blue3]\n"
        "[blue]=======================================================================[/blue]\n"
        "[white on blue]   W E L C O M E   T O   N D R   D E T E C T I O N   [/white on blue]\n"
        "[blue]=======================================================================[/blue]\n"
    )
    rprint(logo)


def parse_wig(path: Path) -> list[tuple[str, int, int, float]]:
    """Parse a WIG file into genomic intervals.

    Supports both fixedStep and variableStep formats and returns 0-based inclusive
    start/end coordinates.

    Args:
        path: Path to the WIG file.

    Returns:
        A list of tuples in the form (chrom, start, end, value).
    """
    out = []
    chrom = None
    mode = None
    span = 1
    fixed_next_pos = None

    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith(("track", "browser")):
                continue
            if line.startswith("variableStep"):
                # parse attrs like: variableStep chrom=chr1 span=150
                mode = "variable"
                attrs = dict(tok.split("=") for tok in line.split()[1:])
                chrom = attrs.get("chrom")
                span = int(attrs.get("span", 1))
                continue
            if line.startswith("fixedStep"):
                mode = "fixed"
                attrs = dict(tok.split("=") for tok in line.split()[1:])
                chrom = attrs.get("chrom")
                start = int(attrs.get("start", "1"))  # 1-based
                step = int(attrs.get("step", "1"))
                span = int(attrs.get("span", "1"))
                fixed_next_pos = start  # 1-based
                continue

            # data line
            if mode == "variable":
                # format: position value
                pos_str, val_str = line.split()
                pos = int(pos_str) - 1           # convert to 0-based
                val = float(val_str)
                out.append((chrom, pos, pos + span - 1, val))  # end inclusive
            elif mode == "fixed":
                val = float(line)
                pos0 = fixed_next_pos - 1
                out.append((chrom, pos0, pos0 + span - 1, val))
                fixed_next_pos += step
            else:
                # plain two-column wig (chrom,pos,value) fallback
                parts = line.split()
                if len(parts) == 3:
                    chrom = parts[0]
                    pos = int(parts[1]) - 1
                    val = float(parts[2])
                    out.append((chrom, pos, pos + span - 1, val))
    return out

def get_NDRs(a: np.ndarray, thresh: float = 0.4, min_len: int = 70) -> tuple[list[np.ndarray], list[tuple[int, int]]]:
    """Detect NDR segments and return groups with index ranges.

    Args:
        a: Signal array.
        thresh: Threshold above which values are considered NDR.
        min_len: Minimum inclusive length for a segment to be kept.

    Returns:
        A tuple of (groups, ranges) where groups are arrays of values and ranges
        are inclusive (start, end) indices.
    """
    a = np.asarray(a)
    mask = a > thresh
    if not mask.any():
        return [], []
    diff = np.diff(mask.astype(int))
    starts = np.where(diff == 1)[0] + 1
    ends   = np.where(diff == -1)[0] + 1
    if mask[0]:
        starts = np.r_[0, starts]
    if mask[-1]:
        ends = np.r_[ends, a.size]
    ends_inclusive = ends - 1
    ranges = list(zip(starts.tolist(), ends_inclusive.tolist()))  # end inclusive
    groups = [a[s:e+1] for s,e in ranges]
    # filter groups by inclusive length > min_len
    filtered = [(g, (s,e)) for g, (s,e) in zip(groups, ranges) if (e - s + 1) > min_len]
    if not filtered:
        return [], []
    groups_f, ranges_f = zip(*filtered)
    return list(groups_f), list(ranges_f)


def main(src_wigfile: Path, output_path: Path):
        """Run NDR detection for a WIG file and write a BED file.

        Args:
                src_wigfile: Path to the input WIG file.
                output_path: Directory where the BED file will be written.

        Raises:
                FileNotFoundError: If the source WIG file does not exist.
        """
        render_welcome_logo()

        log.info(f"Source WIG file: {src_wigfile}")
        if not src_wigfile.exists():
                log.error(f"Source WIG file does not exist: {src_wigfile}")
                raise FileNotFoundError(f"Source WIG file not found: {src_wigfile}")

        output_bedfile = output_path / (src_wigfile.stem + "_NDRs.bed")
        log.info(f"Output BED file: {output_bedfile}")
        if not output_path.exists():
                log.warning(f"Output directory does not exist and will be created: {output_path}")
                output_path.mkdir(parents=True, exist_ok=True)

        log.info("Parsing WIG file...")
        wig_data = parse_wig(src_wigfile)

        log.info("Transforming WIG data to DataFrame...")
        wig_df = pd.DataFrame(wig_data, columns=["chrom", "start", "end", "value"])

        log.info("List all the supercontigs in the WIG file:")
        supercontigs = wig_df["chrom"].unique().tolist()

        with open(output_bedfile, "w") as out_fh:
                log.info("Processing NDR per supercontig:")
                for supercontig in supercontigs:
                        log.info(f"Processing supercontig: {supercontig}")
                        sub_df = wig_df[wig_df["chrom"] == supercontig]
                        _, indexes = get_NDRs(sub_df["value"], thresh=0.4, min_len=70)

                        for ind in indexes:
                                out_fh.write(f"{supercontig}\t{ind[0]}\t{ind[1]}\n")


if __name__ == "__main__":
    typer.run(main)