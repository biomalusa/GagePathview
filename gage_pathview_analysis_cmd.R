#! /usr/bin/env Rscript

suppressMessages(require(docopt))

#---------------Functions---------------------

source("gage_pathview_analysis_bb.R")

#-------Input parsing----------
'Usage:
gage_pathview_analysis_cmd.R [options]

Options:
-i Name of file with input data/genes
-p Value of p.geomean cutoff in Gage analysis [default: 0.25]
-o Name of output folder where results are saved
' -> doc

opts <- docopt(doc)

input_file <- opts$i
# id <- opts$id
# kegg_pathway <- opt$k
pgeomean <- as.numeric(opts$p)
output_folder <- opts$o

# gage_pathview(input_file = input_file,
#               id = id,
#               kegg_pathway = kegg_pathway,
#               pgeomean = pgeomean,
#               output_folder = output_folder)

gage_pathview(input_file = input_file,
              pgeomean = pgeomean,
              output_folder = output_folder)
