USRCAT similarity workflow
==========================


Aiming to build a tool for comparing 3D-shapes of small organic bioactive molecules, we exploited a concept of invariant 3D-fingerprints called [Ultrafast Shape Recognition with Credo Atom Types (USRCAT)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3505738/pdf/1758-2946-4-27.pdf).
Using this approach, we have developed a compherensive workflow, which would allow to carry out 3D shape similarity searches in a fairly large chemical space (8M compounds).

Dependencies
====
Required for Library
===
- `rdkit`
- `ipywidgets`
- `tqdm`
- `notebook`
- `wget`

Setup and Installation
======

To use this workflow under conda environment properly, please, execute the following commands:
 
- `conda create -n usrcat python=3.8`
- `conda install -c conda-forge rdkit`
- `conda install -c anaconda ipywidgets`
- `conda install -c conda-forge tqdm`
- `conda install -c anaconda notebook`
- `conda install -c anaconda wget`

Or use usrcatFP.yml file with comand:

- `conda env create -f usrcatFP.yml`

Databases
======

Also, we prepared a few databases with USRCAT fingerprints to screen Chemspace targeted libraries:
 
- [GPCRs Targeted Library](https://chemspacecom-my.sharepoint.com/:u:/g/personal/m_protopopov_chem-space_com/EeAYk44LOtJJqQwsVPqKDx8BeaCQBizh155_tUfPj5l8aw?e=Gd9eur)
- [Allosteric Protein Kinases Targeted Library](https://chemspacecom-my.sharepoint.com/:u:/g/personal/m_protopopov_chem-space_com/EZniFi17gP9OtUHokSINFA8BNgynKLXRWdocIt6eRK3tgA?e=9aPcVh)
- [General Protein Kinases Targeted Library](https://chemspacecom-my.sharepoint.com/:u:/g/personal/m_protopopov_chem-space_com/EUyd0ux19QtMkzQta991MUABm7eHqkqREOc0ePFmaZZ9OQ?e=IvK6ms)
- [Ion Channels Targeted Library](https://chemspacecom-my.sharepoint.com/:u:/g/personal/m_protopopov_chem-space_com/EeSKrA8uYwZOgo35QNwLgQ8BCDlFHEY7f-TU7pb8OkTE5g?e=VZIQp9)

And 
- [ChEMBL active compounds and their analogs library](https://chemspacecom-my.sharepoint.com/:u:/g/personal/m_protopopov_chem-space_com/EaHyY4pRqeRIs9ykQN3imhEBvl8_0OvQA0PHO0KH1Gue5w?e=MfLJXM)

All these libraries could be a good starting point for related drug discovery projects.

More info about compound sets one can find at [Chemspace website](https://chem-space.com/compounds#screening-compounds)



