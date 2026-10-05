"""Ensemble manifest and CP2K result summary: builders and a small validator."""
from __future__ import annotations

from datetime import datetime, timezone

ENSEMBLE_SCHEMA = "qd.ensemble.v1"
CP2K_SCHEMA = "qd.cp2k.v1"


class SchemaError(ValueError):
    pass


# Required keys and their types, per format.
_REQUIRED = {
    ENSEMBLE_SCHEMA: {
        "schema": str, "record_id": str, "fingerprint": str, "components": dict,
        "model": dict, "hessian": dict, "temperatures_K": list, "samples_per_temperature": int,
        "seed": int, "frequency_cutoff_cm1": (int, float), "files": dict, "created": str,
    },
    CP2K_SCHEMA: {
        "schema": str, "record_id": str, "level": dict, "geo_opt": dict, "window": dict,
        "frontier": dict, "files": dict, "created": str,
    },
}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def validate(doc: dict) -> dict:
    """Raise SchemaError when `doc` lacks a required key or has the wrong type; return it otherwise."""
    schema = doc.get("schema")
    if schema not in _REQUIRED:
        raise SchemaError(f"unknown schema {schema!r}")
    for key, typ in _REQUIRED[schema].items():
        if key not in doc:
            raise SchemaError(f"{schema}: missing {key!r}")
        if not isinstance(doc[key], typ):
            raise SchemaError(f"{schema}: {key!r} should be {typ}, is {type(doc[key]).__name__}")
    if schema == ENSEMBLE_SCHEMA:
        if not doc["temperatures_K"] or any(t <= 0 for t in doc["temperatures_K"]):
            raise SchemaError(f"{schema}: temperatures_K must be a non-empty list of positive temperatures")
        if doc["samples_per_temperature"] < 1:
            raise SchemaError(f"{schema}: samples_per_temperature must be at least 1")
    return doc


def ensemble_manifest(
    *, record_id: str, fingerprint: str, components: dict, model: dict, hessian: dict,
    temperatures_K, samples_per_temperature: int, seed: int, frequency_cutoff_cm1: float = 20.0,
    files: dict | None = None,
) -> dict:
    """
    Manifest of one structure's Wigner ensembles.

    model:   {"name": "MACE-MH-1", "head": "omat_pbe", "sha256": ...}
    hessian: {"method": "analytic" | "finite_differences", "sha256": ...}
    files:   {"optimised": "optimised.xyz", "modes": "modes.npz",
              "samples": {"300": "wigner_300K.extxyz", ...},
              "embeddings": {"300": "embeddings_300K.npz", ...},
              "dft": {"300": "dft_300K.extxyz", ...}}
    """
    return validate({
        "schema": ENSEMBLE_SCHEMA, "record_id": record_id, "fingerprint": fingerprint,
        "components": components, "model": model, "hessian": hessian,
        "temperatures_K": [float(t) for t in temperatures_K],
        "samples_per_temperature": int(samples_per_temperature), "seed": int(seed),
        "frequency_cutoff_cm1": float(frequency_cutoff_cm1), "files": files or {}, "created": _now(),
    })


def cp2k_summary(
    *, record_id: str, level: dict, geo_opt: dict, window: dict, frontier: dict, files: dict | None = None,
) -> dict:
    """
    Summary of the CP2K PBE track of one structure (goes to the webapp).

    level:    {"code": "CP2K 2026.2", "functional": "PBE", "method": "GAPW",
               "basis": "DZVP-MOLOPT-GGA-GTH", "potential": "GTH-GGA"}
    geo_opt:  {"converged": bool, "energy_ev": float, "steps": int}
    window:   output of qd_workflows.cp2k.mo_window
    frontier: {"homo_ev": .., "lumo_ev": .., "gap_ev": .., "homo_soc_ev": .., "lumo_soc_ev": .., "gap_soc_ev": ..}
    """
    return validate({
        "schema": CP2K_SCHEMA, "record_id": record_id, "level": level, "geo_opt": geo_opt,
        "window": window, "frontier": frontier, "files": files or {}, "created": _now(),
    })
