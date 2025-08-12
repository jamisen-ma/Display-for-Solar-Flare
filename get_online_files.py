"""Helpers for scraping and assembling online file URLs."""

from __future__ import annotations

from typing import Iterable, List

from bs4 import BeautifulSoup
import requests

from logging_setup import get_logger


logger = get_logger(__name__)


def _find_files(url: str) -> list[str]:
    """Return hrefs from a simple directory listing page."""
    try:
        response = requests.get(url, timeout=20)
        response.raise_for_status()
    except Exception as exc:  # broad but intentional around network
        logger.warning("Failed to GET %s: %s", url, exc)
        return []

    soup = BeautifulSoup(response.text, features="lxml")
    return [a.get("href", "") for a in soup.find_all("a")]


def get_image_urls(dates_list: Iterable[str]) -> List[str]:
    """Build list of RHESSI qlook image URLs for provided YYYY-MM-DD dates."""
    file_path_list: List[str] = []
    for date_text in dates_list:
        date_path = str(date_text).replace("-", "/")
        base = (
            f"https://hesperia.gsfc.nasa.gov/hessidata/metadata/qlook_image/{date_path}/"
        )
        for link in _find_files(base):
            if "fsimg" in str(link):
                file_path_list.append(f"{base}{link}")
    logger.info("Collected %d image URLs", len(file_path_list))
    return file_path_list