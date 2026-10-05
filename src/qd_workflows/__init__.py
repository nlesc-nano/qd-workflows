"""Workflows tying QD_builder, Orchestr.AI, CP2K and QDEX together."""
from pathlib import Path

__version__ = "0.1.0"

PACKAGE_DIR = Path(__file__).resolve().parent
TEMPLATES = PACKAGE_DIR / "templates"
REPO_DIR = PACKAGE_DIR.parents[1]


def render(template: str, values: dict) -> str:
    """Fill `{{NAME}}` placeholders; every placeholder must get a value."""
    out = template
    for key, value in values.items():
        out = out.replace("{{" + key + "}}", str(value))
    if "{{" in out:
        start = out.index("{{")
        raise KeyError(f"unfilled placeholder {out[start:out.index('}}', start) + 2]}")
    return out
