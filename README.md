# biosecurity_and_overdispersion
This repository contains anonymized data and R scripts to reproduce the analyses from the manuscript titled "Efficient border biosecurity inspection leverages superspreading to reduce biological invasion risk".

In the **data** folder, you will find the data_all_pathways.csv file which contains border biosecurity inspection data associated with the manuscript. Each row represent the result of a border biosecurity inspection on an incoming consignment. The columns are:
* data_type: The data source of the data which can be one of five distinct categories (Germplasm and seed viroids imports to Australia, fresh produce imports to New Zealand, cut flower imports to New Zealand, and pelleted seed imports to New Zealand).
* group_id: Anonymized pathway name. Each pathway represents a distinct combination of commodity and country of origin.
* k: Number of infested units found in the inspection.
* n: Inspection sample size.
* N: Number of units in the consignement (i.e., consignment size). 

The **R** folder contains 14 R scripts that reproduce the analyses and figures shown in the manuscript and supplementary material. The **functions** folder contains custom R functions that are sourced and used in specific R scrips.
The results of the R scripts are stored in the **figs** and **outputs** folders. These results include figures, tables, and model fits that are used in the paper.
Cite the code and data: 
[![DOI](https://zenodo.org/badge/668994138.svg)](https://zenodo.org/badge/latestdoi/668994138)
