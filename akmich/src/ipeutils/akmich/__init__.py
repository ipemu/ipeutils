from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("ipeutils-akmich")
except PackageNotFoundError:
    # Výchozí hodnota, pokud balíček ještě nebyl nainstalován (např. při lokálním vývoji)
    __version__ = "unknown"


from .akmich import FAS, ASD, CSD, trFAS, trASD, ko_smoothing, koc_smoothing

__all__ = ['FAS', 'ASD', 'CSD', 'trFAS', 'trASD',
           'ko_smoothing']
