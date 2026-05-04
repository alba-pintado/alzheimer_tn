# README

This folder includes all the Jupyter Notebooks and R scripts necessary to reprocuce the analysis described in the manuscript "Lack of evidence for increased transcriptional noise in aged tissues"
(https://www.biorxiv.org/content/10.1101/2022.05.18.492432v1).

##### Structure ######

# ------ Jupyter notebooks-------

The "notebooks" directory contains the Jupyter Notebooks, where datasets are processed and the transcriptional noise associated with aging is measured using Scallop and Decibel. 
In order to run them, you need to first install both Scallop and Decibel. Installation and usage instructions can be found in their respective official documentation sites. 

---Documentation--- 
Scallop documentation: https://scallop.readthedocs.io/en/latest/index.html
Decibel documentation: https://decibel.readthedocs.io/en/latest/index.html

---Repositories---
Scallop repository: https://gitlab.com/olgaibanez/scallop/
Decibel repository: https://gitlab.com/olgaibanez/decibel/

All notebooks are provided as a reproducible notebook (ipynb) and an html for visualization purposes. 

# --------- R scripts -----------
The "R_scripts" directory contains the R scripts used for the cell type composition analysis using generalized linear models (GLM). 


Directories: 
* mouse_lung: Angelidis, Kimmel, TMS FACS lung and TMS droplet lung
* human_lung: Travaglini (HLCA) and Raredon (human samples)
* seven_datasets: the datasets used in Figure 2 (different mouse/human tissues). NOTE: only 5 datasets are included in this directory (Salzer, Solé-Boldó, Martínez-Jimenez, Ximerakis and Enge) because the remaining two (Angelidis and Kimmel) are already included in the mouse_lung directory. 
* figures: notebook that generates the final figures shown in the article


