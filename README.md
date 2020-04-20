# GagePathview
Mapping genes from transcriptomic data onto Kegg pathways using R Gage and Pathview packages.

# Description
The scripts are used to analyze a set of genes from a differential gene expression analysis in order to make pathway analysis using Gage and Pathview packages from R/Bioconductor.

This tool is based on the following tutorial from Stephen Turner: https://www.r-bloggers.com/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/ starting from Pathway analysis part and using as input a list of gene IDs from DESeq2 analysis either as Gene Name or as ENSEMBL ID.

* Gage package description: https://bioconductor.org/packages/release/bioc/html/gage.html
* Pathview package description: https://bioconductor.org/packages/release/bioc/html/pathview.html

# Usage
The main script is launched commandline (from linux terminal usually) with Rscript like this:
'Rscript gp/gage_pathview_analysis_cmd.R' and following instructions on given options.

