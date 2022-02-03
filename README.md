USRCAT similarity workflow
==========================


Using this fingerprint, we discovered a workflow based on [Ultrafast Shape Recognition with Credo Atom Types (USRCAT)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3505738/pdf/1758-2946-4-27.pdf), to allow 3D shape similarity search in a huge part of chemical space (8M compounds).

Dependencies
Required for Library

- `rdkit`
- `ipywidgets`
- `tqdm`
- `notebook`
- `wget`

Setup and Installation
======

To use our workflow at conda environment, use the following command 
 
- `conda create -n usrcat python=3.8`
- `conda install -c conda-forge rdkit`
- `conda install -c anaconda ipywidgets`
- `conda install -c conda-forge tqdm`
- `conda install -c anaconda notebook`
- `conda install -c anaconda wget`

Databases
======

Also, we prepare databases of USRCAT fingerprint to screen Chemspace targeted library
 
- [GPCRs Targeted Library](https://chemspacecom-my.sharepoint.com/:u:/g/personal/m_protopopov_chem-space_com/EeAYk44LOtJJqQwsVPqKDx8BeaCQBizh155_tUfPj5l8aw?e=Gd9eur)
- [Allosteric Protein Kinases Targeted Library](https://chemspacecom-my.sharepoint.com/:u:/g/personal/m_protopopov_chem-space_com/EZniFi17gP9OtUHokSINFA8BNgynKLXRWdocIt6eRK3tgA?e=9aPcVh)
- [General Protein Kinases Targeted Library](https://chemspacecom-my.sharepoint.com/:u:/g/personal/m_protopopov_chem-space_com/EUyd0ux19QtMkzQta991MUABm7eHqkqREOc0ePFmaZZ9OQ?e=IvK6ms)
- [Ion Channels Targeted Library](https://chemspacecom-my.sharepoint.com/:u:/g/personal/m_protopopov_chem-space_com/EeSKrA8uYwZOgo35QNwLgQ8BCDlFHEY7f-TU7pb8OkTE5g?e=VZIQp9)

And 
- [ChEMBL active compounds and their analogs library](https://chemspacecom-my.sharepoint.com/:u:/g/personal/m_protopopov_chem-space_com/EaHyY4pRqeRIs9ykQN3imhEBvl8_0OvQA0PHO0KH1Gue5w?e=MfLJXM)

All these libraries could be a good starting point for related drug discovery projects.

More info about compounds sets you could found on [Chemspace website](https://chem-space.com/compounds#screening-compounds)



