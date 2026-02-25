# app.py

import io
import re
import zipfile
import streamlit as st
from src.core import parse_rdf_reactions

st.set_page_config(page_title="RDF file viewer", layout="centered")
st.title("Chemical Reaction Viewer (RDF → SVG images)")

def sanitize_filename(name: str, default: str = "reaction") -> str:
    """
    Make a safe filename by stripping non-filename characters.
    """
    name = name or default
    # Replace spaces with underscores and remove unsafe characters
    name = re.sub(r"\s+", "_", name.strip())
    name = re.sub(r"[^a-zA-Z0-9._-]", "", name)
    # Avoid empty names
    return name if name else default

def build_svg_zip(reactions):
    """
    Build an in-memory ZIP archive of all available SVGs.
    Returns a BytesIO ready for st.download_button or None if there are no SVGs.
    """
    # Collect only entries with an SVG string
    items = [(idx + 1, r) for idx, r in enumerate(reactions) if r.get("svg")]
    if not items:
        return None

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        used_names = set()
        for i, r in items:
            # Construct a readable file name using class or smiles
            base = r.get("reaction_class") or r.get("smiles") or f"reaction_{i}"
            base = sanitize_filename(base)
            # Ensure uniqueness in case of duplicates
            filename = f"{base}.svg"
            k = 2
            while filename in used_names:
                filename = f"{base}_{k}.svg"
                k += 1
            used_names.add(filename)

            svg_bytes = r["svg"].encode("utf-8")
            zf.writestr(filename, svg_bytes)

    buffer.seek(0)
    return buffer

uploaded_file = st.file_uploader("Upload an RDF file", type=["rdf"])

if uploaded_file is not None:
    # Keep special symbols in conditions like Ti(OiPr)₄
    rdf_text = uploaded_file.read().decode("utf-8", errors="ignore")

    reactions = parse_rdf_reactions(rdf_text)

    st.write(f"### Found {len(reactions)} reactions")

    # ---------- NEW: Download all SVGs as a ZIP ----------
    zip_buffer = build_svg_zip(reactions)
    if zip_buffer:
        st.download_button(
            label="⬇️ Download ALL SVGs as ZIP",
            data=zip_buffer,
            file_name="reactions_svg.zip",
            mime="application/zip",
            type="primary",
        )
    else:
        st.info("No SVGs available to zip (no reactions with rendered SVG).")

    # ---------- Per-reaction rendering as before ----------
    for i, r in enumerate(reactions, start=1):
        st.markdown(f"## Reaction {i}")

        if r.get("svg"):
            st.markdown(
                f"""
                <div style="border:1px solid #ddd; padding:8px; border-radius:6px;">
                    {r["svg"]}
                </div>
                """,
                unsafe_allow_html=True,
            )

            st.download_button(
                label="⬇️ Download SVG",
                data=r["svg"].encode("utf-8"),
                file_name=f"reaction_{i}.svg",
                mime="image/svg+xml",
                key=f"dl_svg_{i}",
            )
        else:
            st.warning("No image generated for this reaction (SVG unavailable).")

        if r.get("smiles"):
            st.markdown(f"**SMILES:** `{r['smiles']}`")

        if r.get("conditions"):
            st.markdown(f"**Reaction Conditions:** {r['conditions']}")

        if r.get("reaction_class"):
            st.markdown(f"**Reaction Class:** {r['reaction_class']}")

        st.markdown("---")

# =========================
# Sidebar controls
# =========================
st.sidebar.image("images/BNNLab_v3.png")
st.sidebar.header("Acknowledgements")
st.sidebar.write("This web tool was built to read and visualised RDF files from Synthia retrosynthetic tools. It may not work with rdf files from other sources.")
st.sidebar.header("Disclaimer")
st.sidebar.write("This software was developed by BNNLab, with all rights reserved. It is offered 'as is', without warranty of any kind, express or implied. The user assumes all risk for any malfunctions, errors, or damages resulting from the use of this software. The creator is not responsible for any direct or indirect loss arising from its use.")