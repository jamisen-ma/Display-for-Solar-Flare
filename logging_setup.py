"""
Centralized logging configuration for the Solar Flare Display project.

Usage:
    from logging_setup import get_logger
    logger = get_logger(__name__)
    logger.info("message")

Configuration can be controlled via environment variables:
    LOG_LEVEL: DEBUG|INFO|WARNING|ERROR|CRITICAL (default: INFO)
    LOG_FORMAT: python logging format string (optional)
"""

from __future__ import annotations

import logging
import os
from typing import Optional


_DEFAULT_FORMAT = (
    "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
)


def _configure_root_logger(level: str, fmt: Optional[str]) -> None:
    if logging.getLogger().handlers:
        # Already configured
        return
    log_level = getattr(logging, level.upper(), logging.INFO)
    log_format = fmt or _DEFAULT_FORMAT
    logging.basicConfig(level=log_level, format=log_format)


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """Return a logger configured with project-wide defaults.

    Parameters
    ----------
    name: Optional[str]
        The logger name. If None, returns the root logger.
    """
    level = os.getenv("LOG_LEVEL", "INFO")
    fmt = os.getenv("LOG_FORMAT", None)
    _configure_root_logger(level, fmt)
    return logging.getLogger(name)


