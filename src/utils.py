from typing import Optional
import logging
from logging.handlers import RotatingFileHandler
from rich.logging import RichHandler


def set_logging(
    log_file: Optional[str] = None,
    max_bytes: int = 10 * 1024 * 1024,
    backup_count: int = 5,
):
    """
    Set up logging configuration for the application.
    This function configures the logging format, level, and handlers.

    Args:
        log_file: Optional path to a log file. If provided, logs are written to this file as well
                  as to the console.
        max_bytes: Maximum size in bytes before rotating the log file.
        backup_count: Number of rotated backup files to keep.
    """
    formatter = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    handlers = [RichHandler()]

    if log_file:
        file_handler = RotatingFileHandler(
            log_file, maxBytes=max_bytes, backupCount=backup_count
        )
        file_handler.setFormatter(logging.Formatter(formatter, datefmt="[%X]"))
        handlers.append(file_handler)

    logging.basicConfig(
        level=logging.INFO, format=formatter, datefmt="[%X]", handlers=handlers
    )
