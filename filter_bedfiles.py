import logging
from rich import print as rprint
from rich.logging import RichHandler
from pathlib import Path
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
    """Render a colorful BED filtering banner using Rich."""
    logo = (
        "\n\n"
        "[magenta]  ____  ____  ____  ____  _     ___  _     ____[/magenta]\n"
        "[magenta] /  __\\\\/  _ \\\\/  _ \\\\/  _ \\\\/ \\\\  /  _ \\\\/ \\\\  /|/  _ \\\\[/magenta]\n"
        "[magenta] |  \\\\  | / \\\\|| / \\\\|| / \\\\|| |  | / \\\\||  \\\\/ || / \\\\|[/magenta]\n"
        "[magenta] |  /_ | \\_/|| \\_/|| \\_/|| |\\\\ \\ | \\_/|| |\\\\  /|| \\_/|[/magenta]\n"
        "[magenta] \\____/\\____/\\____/\\____/\\_/ \\\\_\\____/\\_/ \\\\_/\\____/[/magenta]\n"
        "[cyan]--------------------------------------------------------------[/cyan]\n"
        "[sky_blue3]Filtering BED intervals by length[/sky_blue3]    "
        "[green]::[/green]  [white]keeping long regions[/white]\n"
        "[cyan]--------------------------------------------------------------[/cyan]\n"
    )
    rprint(logo)


def filter_short_seq(file_path: Path, min_len: int = 10):
    """Filter BED rows by interval length and write a new BED file.

    Reads a 3-column BED file (chr, start, end), keeps rows where
    (end - start) is greater than `min_len`, and writes the filtered
    records to a new file with the suffix `.filtered.bed`.

    Args:
        file_path: Path to the input BED file.
        min_len: Minimum length required to keep a row.
    """
    df = pd.read_csv(file_path, delimiter="\t", header=None, names=["chr", "start", "end"])
    df["valid"] = (df["end"] - df["start"]) > min_len
  
    outpath = file_path.with_suffix(".filtered.bed")
    df[df["valid"]].drop(columns=["valid"]).to_csv(outpath, header=None, index=None, sep="\t")


def main(filepath: Path, min_len: int = 10):
    """Run the default filtering task for the practice dataset."""
    render_welcome_logo()
    filter_short_seq(filepath, min_len)

if __name__ == "__main__":
    typer.run(main)