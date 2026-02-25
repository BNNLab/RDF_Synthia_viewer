# src/core.py

import re
from typing import List, Dict, Optional

from rdkit.Chem import rdChemReactions

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
# SVG drawer detection (NO import of rdkit.Chem.Draw)
# --------------------------------------------------------
def _get_svg_drawer(width: int, height: int):
    """
    Try to obtain MolDraw2DSVG from rdkit.Chem.rdMolDraw2D only.
    Do NOT import rdkit.Chem.Draw to avoid failures on some builds.
    """
    try:
        from rdkit.Chem import rdMolDraw2D  # safe path if present
        return rdMolDraw2D.MolDraw2DSVG(width, height)
    except Exception:
        return None  # drawer not available in this RDKit build


# --------------------------------------------------------
# RXN block → SVG (preferred path)
# --------------------------------------------------------
def reaction_rxnblock_to_svg(
    rxn_block: str,
    width: int = 600,
    height: int = 280,
) -> Optional[str]:
    """
    Create an SVG string for a reaction given a full RXN block.
    Uses rdMolDraw2D if available; otherwise returns None.
    """
    try:
        rxn = rdChemReactions.ReactionFromRxnBlock(rxn_block)
        if rxn is None:
            return None

        drawer = _get_svg_drawer(width, height)
        if drawer is None:
            return None  # cannot draw in this environment

        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return None


# --------------------------------------------------------
# SMILES → SVG (fallback only, using same drawer)
# --------------------------------------------------------
def reaction_smiles_to_svg(
    smiles: str,
    width: int = 600,
    height: int = 280,
) -> Optional[str]:
    """
    Create an SVG string for a reaction given reaction SMILES.
    Uses rdMolDraw2D if available; otherwise returns None.
    """
    try:
        rxn = rdChemReactions.ReactionFromSmarts(smiles, useSmiles=True)
        if rxn is None:
            return None

        drawer = _get_svg_drawer(width, height)
        if drawer is None:
            return None  # cannot draw in this environment

        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return None


# --------------------------------------------------------
# Main RDF parser (unchanged interface)
# --------------------------------------------------------
def parse_rdf_reactions(rdf_text: str) -> List[Dict[str, Optional[str]]]:

    entries = rdf_text.split("$RFMT")
    results: List[Dict[str, Optional[str]]] = []

    for entry in entries:

        if "$RXN" not in entry and "$DTYPE SMILES" not in entry:
            continue

        smiles = _extract_dtype(entry, "SMILES")
        conditions = _extract_dtype(entry, "Reaction Conditions")
        reaction_class = _extract_dtype(entry, "Name")

        svg: Optional[str] = None

        # Prefer RXN block rendering
        rxn_block = _extract_rxn_block(entry)
        if rxn_block:
            svg = reaction_rxnblock_to_svg(rxn_block)

        # Fallback to SMILES if RXN fails
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
