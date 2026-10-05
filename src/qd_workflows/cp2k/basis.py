"""Basis sets and potentials from CP2K data files (BASIS_MOLOPT_UZH, POTENTIAL_UZH)."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Kind:
    element: str
    basis: str              # e.g. DZVP-MOLOPT-GGA-GTH-q12
    potential: str          # e.g. GTH-GGA-q12
    valence: int
    functions: int          # spherical basis functions per atom


def _basis_entries(text: str, element: str, family: str):
    """(name, line index) of every basis of `element` whose name starts with `family`-q."""
    lines = text.splitlines()
    pattern = re.compile(rf"^{re.escape(family)}-q(\d+)$")
    for i, line in enumerate(lines):
        tokens = line.split()
        if tokens and tokens[0] == element:
            for tok in tokens[1:]:
                if pattern.match(tok):
                    yield tok, i


def count_functions(lines: list[str], start: int) -> int:
    """Spherical functions of the basis whose header is lines[start] (CP2K basis format)."""
    i = start + 1
    nsets = int(lines[i].split()[0])
    i += 1
    total = 0
    for _ in range(nsets):
        head = [int(x) for x in lines[i].split()]
        _n, lmin, lmax, nexp = head[:4]
        nshell = head[4:4 + lmax - lmin + 1]
        total += sum(ns * (2 * l + 1) for l, ns in zip(range(lmin, lmax + 1), nshell))
        i += 1 + nexp
    return total


def kind_for(element: str, basis_file, potential_file, basis_family: str, potential_family: str) -> Kind:
    """
    The kind for `element`: the one basis of `basis_family` in the basis file fixes the
    valence charge, and the potential `potential_family`-q<that charge> must exist.
    """
    text = Path(basis_file).read_text()
    entries = list(_basis_entries(text, element, basis_family))
    if len(entries) != 1:
        found = ", ".join(n for n, _ in entries) or "none"
        raise ValueError(f"{element}: expected one {basis_family}-qN basis in {basis_file}, found {found}")
    name, line = entries[0]
    valence = int(name.rsplit("-q", 1)[1])
    potential = f"{potential_family}-q{valence}"
    pot_text = Path(potential_file).read_text()
    if not re.search(rf"^{re.escape(element)}\s.*\b{re.escape(potential)}\b", pot_text, re.M):
        raise ValueError(f"{element}: no {potential} in {potential_file}")
    return Kind(element, name, potential, valence, count_functions(text.splitlines(), line))


def kinds_for(elements, basis_file, potential_file, basis_family: str, potential_family: str) -> dict[str, Kind]:
    return {e: kind_for(e, basis_file, potential_file, basis_family, potential_family) for e in sorted(set(elements))}


def kind_blocks(kinds: dict[str, Kind]) -> str:
    return "\n".join(
        f"    &KIND {k.element}\n      BASIS_SET  {k.basis}\n      POTENTIAL  {k.potential}\n    &END"
        for k in kinds.values()
    )
