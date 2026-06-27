from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class CTBlock:
    struc_tag: str    # value of s_fep_struc_tag CT property
    raw_block: str    # full text of the f_m_ct { ... } block


def parse_mae(path: str | Path) -> list[CTBlock]:
    """
    Parse a .mae file and return CT blocks that have an s_fep_struc_tag property.
    CTs without this property are silently skipped.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"MAE file not found: {path}")
    text = path.read_text(encoding="latin-1")
    return _extract_tagged_cts(text)


def _extract_tagged_cts(text: str) -> list[CTBlock]:
    blocks = []
    raw_blocks = _split_ct_blocks(text)
    for raw in raw_blocks:
        tag = _extract_struc_tag(raw)
        if tag is not None:
            blocks.append(CTBlock(struc_tag=tag, raw_block=raw))
    return blocks


def _split_ct_blocks(text: str) -> list[str]:
    """Return list of raw f_m_ct block texts (including the opening line)."""
    blocks = []
    lines = text.splitlines(keepends=True)
    current: list[str] = []
    in_block = False
    depth = 0

    for line in lines:
        stripped = line.strip()
        if stripped == "f_m_ct {":
            if in_block and current:
                blocks.append("".join(current))
            current = [line]
            in_block = True
            depth = 1
        elif in_block:
            current.append(line)
            depth += stripped.count("{") - stripped.count("}")
            if depth <= 0:
                blocks.append("".join(current))
                current = []
                in_block = False
                depth = 0

    if current:
        blocks.append("".join(current))
    return blocks


def _extract_struc_tag(raw_block: str) -> str | None:
    """
    Extract s_fep_struc_tag value from a CT block.

    MAE CT-level properties follow this layout:
        f_m_ct {
          prop_name_1         <- header: one property name per line
          prop_name_2
         :::
          value_1            <- values: same order as names
          value_2
         :::
    We find s_fep_struc_tag in the header, then read the value at the same index.
    """
    lines = raw_block.splitlines()

    prop_names: list[str] = []
    sep_count = 0
    values: list[str] = []

    state = "header"
    for line in lines[1:]:   # skip "f_m_ct {"
        stripped = line.strip()
        if stripped == ":::":
            sep_count += 1
            if sep_count == 1:
                state = "values"
            elif sep_count == 2:
                break
            continue

        if state == "header":
            if stripped and not stripped.startswith("#") and not stripped.startswith("m_"):
                prop_names.append(stripped)
        elif state == "values":
            if stripped:
                values.append(stripped.strip('"'))

    if "s_fep_struc_tag" not in prop_names:
        return None
    idx = prop_names.index("s_fep_struc_tag")
    if idx >= len(values):
        return None
    return values[idx]
