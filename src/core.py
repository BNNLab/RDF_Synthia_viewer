# src/core.py

import re
from typing import List, Dict, Optional

from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D


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

    We capture from the first '$RXN' up to (but not including) the next line
    that starts with '$' (e.g., $DTYPE, $RIREG, or the next $RFMT), or the end
    of the entry. This avoids pulling in RDF metadata and keeps the RXN block clean.
    """
    m = re.search(r"(\$RXN[\s\S]*?)(?=\r?\n\$|$)", entry, flags=re.MULTILINE)
    return m.group(1).strip() if m else None


def reaction_rxnblock_to_svg(
    rxn_block: str,
    width: int = 600,
    height: int = 280,
) -> Optional[str]:
    """
    Create an SVG string for a reaction given a full RXN block.
    Returns None if parsing or drawing fails.
    """
    try:
        rxn = rdChemReactions.ReactionFromRxnBlock(rxn_block)
        if rxn is None:
            return None
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return None


def reaction_smiles_to_svg(  # kept as a fallback for robustness
    rxn_smiles: str,
    width: int = 600,
    height: int = 280,
) -> Optional[str]:
    """
    Create an SVG string for a reaction given reaction SMILES.
    Returns None if parsing or drawing fails.
    """
    try:
        rxn = rdChemReactions.ReactionFromSmarts(rxn_smiles, useSmiles=True)
        if rxn is None:
            return None
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return None


def parse_rdf_reactions(rdf_text: str) -> List[Dict[str, Optional[str]]]:
    """
    Parse an MDL RDF file and extract:
      - reaction SMILES      (metadata only; unchanged for the UI)
      - reaction conditions  (from $DTYPE Reaction Conditions)
      - reaction class/name  (from $DTYPE Name)
      - reaction SVG         (NOW generated primarily from RXN blocks)
    """
    entries = rdf_text.split("$RFMT")
    results: List[Dict[str, Optional[str]]] = []

    for entry in entries:
        # We’ll try to render from RXN. We still extract SMILES for display.
        if "$RXN" not in entry and "$DTYPE SMILES" not in entry:
            continue

        smiles = _extract_dtype(entry, "SMILES")
        conditions = _extract_dtype(entry, "Reaction Conditions")
        reaction_class = _extract_dtype(entry, "Name")

        # Prefer RXN block rendering for generality
        svg: Optional[str] = None
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
