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


def reaction_smiles_to_svg(
    rxn_smiles: str,
    width: int = 600,
    height: int = 280,
) -> Optional[str]:
    """
    Create an SVG string for a reaction given reaction SMILES.
    Returns None if parsing or drawing fails.
    """
    try:
        # Parse reaction *SMILES* (not SMARTS) correctly:
        rxn = rdChemReactions.ReactionFromSmarts(rxn_smiles, useSmiles=True)
        if rxn is None:
            return None

        # Draw to SVG
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        # Either of these works; DrawReaction is a bit more direct:
        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg
    except Exception:
        return None


def parse_rdf_reactions(rdf_text: str) -> List[Dict[str, Optional[str]]]:
    """
    Parse an MDL RDF file and extract:
      - reaction SMILES
      - reaction conditions
      - reaction class (Name)
      - reaction SVG (as a text string)
    """
    entries = rdf_text.split("$RFMT")
    results: List[Dict[str, Optional[str]]] = []

    for entry in entries:
        if "$DTYPE SMILES" not in entry:
            continue

        smiles = _extract_dtype(entry, "SMILES")
        conditions = _extract_dtype(entry, "Reaction Conditions")
        reaction_class = _extract_dtype(entry, "Name")

        svg = reaction_smiles_to_svg(smiles) if smiles else None

        results.append(
            {
                "smiles": smiles,
                "conditions": conditions,
                "reaction_class": reaction_class,
                "svg": svg,  # <-- pass SVG, not PIL
            }
        )

    return results
