#Code for in vivo Perturb-Seq

This repository contains code for in vivo peturb-seq. 

The main code directory is CodeDirectory_Current, and particularly the RunPipeline.sh code that runs most analysis from the paper automatically.

By subdirectory:
**10X_Perturb_Code:** code for extracting perturb-seq perturbation information from 10X single cell data
**DialoutCode:** code for extracting perturb-seq perturbation information from dialout reads and making them in a nice format
**Dialout_Results:** downstream code to clean up the dialout results
**DECode:** Code to perform DE analysis (with limma) on the perturb-seq data. Not included in the RunPipeline.sh script
**SeuratObjects:** Scripts that go from the raw count matrices to a final Seurat object with the perturb-seq data
**CompareToHuman:** The code that looks for how the WGCNA modules we found are conserved in various human datasets. Not run in that standard pipeline
**UpdateSTM:** Runs code to parse strucutural topic model topics run on the Perturb-seq data for publication
**CellComposition:** Code for comparing cell type composition between perturbations
