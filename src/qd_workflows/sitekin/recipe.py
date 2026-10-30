"""The protocol recipe (YAML): material, builder settings, sizes with their surface treatment, model and step
settings. Paths may carry environment variables."""
from __future__ import annotations

import os
from importlib import resources
from pathlib import Path

import yaml

TEMPLATE = "cdse_cdcl2.yaml"


def template_path(name: str = TEMPLATE) -> Path:
    return Path(str(resources.files("qd_workflows") / "templates" / "sitekin" / name))


def load(path: str | Path | None = None) -> dict:
    """The recipe at path (default: the CdSe/CdCl2 template), with $VARS in material.cif and model.path expanded."""
    p = Path(path) if path else template_path()
    r = yaml.safe_load(p.read_text())
    r["material"]["cif"] = os.path.expandvars(r["material"]["cif"])
    r["model"]["path"] = os.path.expandvars(r["model"]["path"])
    return r


RUN_RECIPE = "recipe.yaml"


def for_run(run: str | Path, path: str | Path | None = None) -> tuple[dict, Path]:
    """The recipe of a run: the given path, else <run>/recipe.yaml (written by the build step), else the template.
    Every step of a run, on the laptop or on the cluster, reads the same file this way."""
    own = Path(run) / RUN_RECIPE
    p = Path(path) if path else (own if own.exists() else template_path())
    return load(p), p


def tag(size: float) -> str:
    """Size label used in every file name: 2.5 -> '2p5'."""
    return f"{float(size):.1f}".replace(".", "p")


def size_entry(r: dict, size_or_tag) -> dict:
    t = size_or_tag if isinstance(size_or_tag, str) else tag(size_or_tag)
    for s in r["sizes"]:
        if tag(s["size"]) == t:
            return s
    raise KeyError(f"size {size_or_tag} is not in the recipe")


def kt(r: dict) -> float:
    return 8.617333262e-5 * float(r["T_K"])
