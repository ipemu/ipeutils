from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("ipeutils-cssio")
except PackageNotFoundError:
    __version__ = "unknown"

from .stream2wfdisc import save_stream_to_wfdisc

__all__ = ['save_stream_to_wfdisc']
