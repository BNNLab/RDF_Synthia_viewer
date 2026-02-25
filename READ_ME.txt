To run this locally, first one needs to create a conda environment:

conda create -n rdf_viewer -c conda-forge python=3.11 rdkit streamlit pillow pandas

Then move to the folder rdf_viewer_v2\ and activate the environment before running streamlit:

conda activate rdf_viewer
streamlit run app.py

If it doesn't automatically open the webpage, you can access it at http://localhost:8501/

There are 2 python files: app.py for the interface, and src\core.py for the image generation underneath. To add the reaction conditions and reaction class, most of the work will be in core.py