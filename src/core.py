# src/core.py

import base64
import io
import re
from typing import List, Dict, Optional, Tuple

from rdkit.Chem import rdChemReactions
from rdkit.Chem import Draw  # for PNG fallback (ReactionToImage)


# -----------------------------
# Helpers: RDF field extraction
# -----------------------------
def _extract_dtype(entry: str, dtype: str) -> Optional[str]:
    """
    Extracts the first line following a $DTYPE <name> ... $DATUM <value>
    """
    pattern = rf"\$DTYPE {re.escape(dtype)}[\s\S]*?\$DATUM ([^\n\r]*)"
    m = re.search(pattern, entry)
    return m.group(1).strip() if m else None


def _extract_rxn_block(entry: str) -> Optional[str]:
    """
    Extract the RXN block for a single reaction entry.

    Capture from the first '$RXN' up to (but not including) the next line
    that begins with '$' (e.g., $DTYPE, $RIREG, or next $RFMT), or end of entry.
    """
    m = re.search(r"(\$RXN[\s\S]*?)(?=\r?\n\$|$)", entry, flags=re.MULTILINE)
    return m.group(1).strip() if m else None


# ------------------------------------------------------
# Drawer compatibility layer for different RDKit builds
# ------------------------------------------------------
def _get_svg_drawer(width: int, height: int):
    """
    Try multiple ways to obtain MolDraw2DSVG across RDKit versions/builds.
    Returns a MolDraw2DSVG instance or None if unavailable.
    """
    # Try the canonical path
    try:
        from rdkit.Chem.Draw import rdMolDraw2D  # type: ignore
        return rdMolDraw2D.MolDraw2DSVG(width, height)
    except Exception:
        pass

    # Some builds require importing from the submodule path
    try:
        from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG  # type: ignore
        return MolDraw2DSVG(width, height)
    except Exception:
        pass

    # Some conda/pip wheels expose MolDraw2DSVG at top-level rdkit.Chem.Draw
    try:
        from rdkit.Chem import Draw as ChemDraw  # re-alias to avoid shadowing
        if hasattr(ChemDraw, "MolDraw2DSVG"):
            return ChemDraw.MolDraw2DSVG(width, height)  # type: ignore
    except Exception:
        pass

    # If everything fails, return None
    return None


def _png_to_svg_wrapper(pil_img, width: int, height: int) -> str:
    """
    Wrap a PNG (PIL Image) inside an SVG container so callers can still
    consume an SVG string even if native SVG drawing isn't available.
    """
    buf = io.BytesIO()
    pil_img.save(buf, format="PNG")
    png_b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    # Build a minimal SVG embedding the PNG as an <image> href
    # Width/height in pixels; browsers handle this fine.
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
        f'<image href="data:image/png;base64,{png_b64}" '
        f'width="{width}" height="{height}" /></svg>'
    )
    return svg


# ----------------------------------------
# RXN → SVG (preferred) and SMILES fallback
# ----------------------------------------
def reaction_rxnblock_to_svg(
    rxn_block: str,
    width: int = 600,
    height: int = 280,
) -> Optional[str]:
    """
    Create an SVG string for a reaction given a full RXN block.
    Tries native SVG drawer first; falls back to PNG-wrapped-in-SVG.
    """
    try:
        rxn = rdChemReactions.ReactionFromRxnBlock(rxn_block)
        if rxn is None:
            return None

        drawer = _get_svg_drawer(width, height)
        if drawer is not None:
            # Native SVG path
            drawer.DrawReaction(rxn)
            drawer.FinishDrawing()
            return drawer.GetDrawingText()

        # Fallback: render PNG and wrap in SVG
        pil_img = Draw.ReactionToImage(rxn, subImgSize=(width, height))
        return _png_to_svg_wrapper(pil_img, width, height)
    except Exception:
        return None


def reaction_smiles_to_svg(
    rxn_smiles: str,
    width: int = 600,
    height: int = 280,
) -> Optional[str]:
    """
    Create an SVG string for a reaction given reaction SMILES.
    Tries native SVG drawer first; falls back to PNG-wrapped-in-SVG.
    """
    try:
        rxn = rdChemReactions.ReactionFromSmarts(rxn_smiles, useSmiles=True)
        if rxn is None:
            return None

        drawer = _get_svg_drawer(width, height)
        if drawer is not None:
            drawer.DrawReaction(rxn)
            drawer.FinishDrawing()
            return drawer.GetDrawingText()

        pil_img = Draw.ReactionToImage(rxn, subImgSize=(width, height))
        return _png_to_svg_wrapper(pil_img, width, height)
    except Exception:
        return None


# ----------------------
# Main parse entry point
# ----------------------
def parse_rdf_reactions(rdf_text: str) -> List[Dict[str, Optional[str]]]:
    """
    Parse an MDL RDF file and extract:
      - reaction SMILES      (metadata only; unchanged for the UI)
      - reaction conditions  (from $DTYPE Reaction Conditions)
      - reaction class/name  (from $DTYPE Name)
      - reaction SVG         (generated primarily from RXN blocks)
    """
    entries = rdf_text.split("$RFMT")
    results: List[Dict[str, Optional[str]]] = []

    for entry in entries:
        # We’ll try RXN first; still extract SMILES for display/fallback.
        if "$RXN" not in entry and "$DTYPE SMILES" not in entry:
            continue

        smiles = _extract_dtype(entry, "SMILES")
        conditions = _extract_dtype(entry, "Reaction Conditions")
        reaction_class = _extract_dtype(entry, "Name")

        svg: Optional[str] = None

        # Prefer RXN block rendering for generality
        rxn_block = _extract_rxn_block(entry)
        if rxn_block:
            svg = reaction_rxnblock_to_svg(rxn_block)

        # Safe fallback: if RXN is missing or unparsable, use SMILES
        if svg is None and smiles:
            svg = reaction_smiles_to_svg(smiles)

        results.append(
            {
                "smiles": smiles,
                "conditions": conditions,
                "reaction_class": reaction_class,
                "svg": svg,
            }
        )

    return results
