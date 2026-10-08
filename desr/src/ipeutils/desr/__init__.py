#!/usr/bin/env python
# coding: utf-8

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("ipeutils-desr")
except PackageNotFoundError:
    __version__ = "unknown"

from .site_response_rest import deconv_site_resp, load_freq_geomean

# from .site_response_rest import cDFT extract_station_name, extract_orientation

__all__ = [
    "deconv_site_resp",
    "load_freq_geomean",
    # "cDFT",
    # "extract_station_name",
    # "extract_orientation",
    ]
