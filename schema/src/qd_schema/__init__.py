"""
Shared file formats of the QD workflows: written by one package, read by
another, so they live here and nowhere else.

Now:   the ensemble manifest (Wigner samples of one structure) and the CP2K
       result summary (DFT track for the webapp).
Later: the library record (record.json) and properties summary, moved from
       QD_builder's builder.library_record once the repositories are published.
"""
from .formats import (
    CP2K_SCHEMA,
    ENSEMBLE_SCHEMA,
    SchemaError,
    cp2k_summary,
    ensemble_manifest,
    validate,
)

__version__ = "0.1.0"
__all__ = ["ENSEMBLE_SCHEMA", "CP2K_SCHEMA", "SchemaError", "ensemble_manifest", "cp2k_summary", "validate"]
