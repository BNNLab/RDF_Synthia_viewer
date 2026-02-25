# src/core.py

import base64
import io
import re
from typing import List, Dict, Optional

from rdkit.Chem import rdChemReactions
from rdkit.Chem import Draw  # needed for fallback PNG drawing


# -----------------------------
# Helper: Extract $DTYPE values
# -----------------------------
def _extract_dtype(entry: str, dtype: str) -> Optional[str]:
    pattern = rf"\$DTYPE {re.escape(dtype)}[\s\S]*?\$DATUM ([^\n\r]*)"
    m = re.search(pattern, entry)
    return m.group(1).strip() if m else None


# ---------------------------------------------------------
# Extract precise RXN block without pulling RDF extra lines
# ---------------------------------------------------------
def _extract_rxn_block(entry: str) -> Optional[str]:
    match = re.search(
        r"(\$RXN[\s\S]*?)(?=\r?\n\$|$)",
        entry,
        flags=re.MULTILINE
    )
    return match.group(1).strip() if match else None


# --------------------------------------------------------
# Robust drawer detection (DOES NOT IMPORT rdMolDraw2D)
# --------------------------------------------------------
def _get_svg_drawer(width: int, height: int):
    """
    Try all known places RDKit may expose MolDraw2DSVG.
    We avoid importing from rdkit.Chem.Draw to stay Streamlit-safe.
    """
    # Option A: rdMolDraw2D available directly under rdkit.Chem
    try:
        from rdkit.Chem import rdMolDraw2D
        return rdMolDraw2D.MolDraw2DSVG(width, height)
    except Exception:
        pass

    # Option B: rdMolDraw2D under rdkit.Chem.Draw (some older conda builds)
    try:
        import rdkit.Chem.Draw.rdMolDraw2D as alt
        return alt.MolDraw2DSVG(width, height)
    except Exception:
        pass

    # Option C: MolDraw2DSVG exposed directly on Draw module
    try:
        if hasattr(Draw, "MolDraw2DSVG"):
            return Draw.MolDraw2DSVG(width, height)
    except Exception:
        pass

    return None  # drawer not available → trigger PNG fallback


# --------------------------------------------------------
# PNG → SVG wrapper (so the app still accepts SVG strings)
# --------------------------------------------------------
def _png_to_svg(pil_img, width: int, height: int) -> str:
    buf = io.BytesIO()
    pil_img.save(buf, format="PNG")
    b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
        f'<image href="data:image/png;base64,{b64}" '
        f'width="{width}" height="{height}" />'
        f'</svg>'
    )


# --------------------------------------------------------
# RXN block → SVG
# --------------------------------------------------------
def reaction_rxnblock_to_svg(
    rxn_block: str,
    width: int = 600,
    height: int = 280,
) -> Optional[str]:

    try:
        rxn = rdChemReactions.ReactionFromRxnBlock(rxn_block)
        if rxn is None:
            return None

        drawer = _get_svg_drawer(width, height)
        if drawer:
            drawer.DrawReaction(rxn)
            drawer.FinishDrawing()
            return drawer.GetDrawingText()

        # PNG fallback for environments without SVG drawer
        pil = Draw.ReactionToImage(rxn, subImgSize=(width, height))
        return _png_to_svg(pil, width, height)

    except Exception:
        return None


# --------------------------------------------------------
# SMILES → SVG (fallback only)
# --------------------------------------------------------
def reaction_smiles_to_svg(
    smiles: str,
    width: int = 600,
    height: int = 280,
) -> Optional[str]:

    try:
        rxn = rdChemReactions.ReactionFromSmarts(smiles, useSmiles=True)
        if rxn is None:
            return None

        drawer = _get_svg_drawer(width, height)
        if drawer:
            drawer.DrawReaction(rxn)
            drawer.FinishDrawing()
            return drawer.GetDrawingText()

        pil = Draw.ReactionToImage(rxn, subImgSize=(width, height))
        return _png_to_svg(pil, width, height)

    except Exception:
        return None


# --------------------------------------------------------
# Main RDF parser (unchanged logic)
# --------------------------------------------------------
def parse_rdf_reactions(rdf_text: str) -> List[Dict[str, Optional[str]]]:

    entries = rdf_text.split("$RFMT")
    results = []

    for entry in entries:

        if "$RXN" not in entry and "$DTYPE SMILES" not in entry:
            continue

        smiles = _extract_dtype(entry, "SMILES")
        conditions = _extract_dtype(entry, "Reaction Conditions")
        reaction_class = _extract_dtype(entry, "Name")

        svg = None

        rxn_block = _extract_rxn_block(entry)
        if rxn_block:
            svg = reaction_rxnblock_to_svg(rxn_block)

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
