# DEP_mass_spec_analysis

The DEP package (differential enrichment analysis of proteomics data, see https://www.bioconductor.org/packages/release/bioc/html/DEP.html) provides tools to analyze proteomics data from MaxQuant output.

Included here is a standard pipeline for using DEP to analyze LFQ/IBAQ signal from a proteinGroups.txt file.

Also included is a function (protGroup_DEP_df_merger.py) which allows for a merge between the proteinGroups.txt file and a DEP results file. This allows the user to compare raw peptide data with the statistical analysis of abundance change between samples.


