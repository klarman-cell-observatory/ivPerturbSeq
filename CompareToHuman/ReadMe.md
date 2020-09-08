# Comparision to human data

Included in this directory are some of the scripts we used for analyasis of WGCNA topics in human data

**ProcessSeurat_CI.R:** This code takes a seurat object with cell types and a table of WGCNA info and checks the associated cell types to see if the module is more highly correlated/conserved than by chance.

**TestEnrichment.DE.R:** This code loads in DE results (both from human autism and from perturb seq) and, for each DE gene, tests to see if the logFC of the perturbations are larger than expected by chance (comparing to genes with similiar expression profiles).
