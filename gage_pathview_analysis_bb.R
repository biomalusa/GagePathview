prepare_data <- function(input_file){
  suppressMessages(require(dplyr))
  suppressMessages(require(stringr))
  suppressMessages(require(AnnotationDbi))
  suppressMessages(require(org.Hs.eg.db))
  
  dge <- read.table(file = input_file,
                    sep = "\t",
                    header = T,
                    check.names = F)
                    # row.names = 1)
  
  #--------id analysis-------- 

  if ("translated_names" %in% colnames(dge)) {
    df <- dge %>% dplyr::select(c(translated_names, log2FoldChange, padj))
    df$geneid = mapIds(org.Hs.eg.db,
                       keys=as.character(df$translated_names), 
                       column="ENTREZID",
                       keytype="SYMBOL",
                       multiVals="first")
    df2 <- df %>% dplyr::filter(!str_detect(translated_names, "ENSG0")) %>% dplyr::filter(!str_detect(geneid, "NA"))
    } else {
    df <- dge %>% dplyr::select(c(gene_name, log2FoldChange, padj))
    df <- df %>% dplyr::filter(padj < 0.05)
    df$geneid = mapIds(org.Hs.eg.db,
                       keys=as.character(df$gene_name), 
                       column="ENTREZID",
                       keytype="SYMBOL",
                       multiVals="first")
    df2 <- df %>% dplyr::filter(!str_detect(gene_name, "ENSG0")) %>% dplyr::filter(!str_detect(geneid, "NA"))
    }

  input_data <- as.data.frame(df2)
  return(input_data)
}

# gage_pathview <- function(input, kegg_pathway, pgeomean, output_folder){
gage_pathview <- function(input_file, pgeomean = 0.25, output_folder = output_folder){
  
  message("Gage/Pathview Analysis Start")
  source("KEGG_gene_annotation.R")
  
  #--------Libraries--------
  message("Loading Libraries")
  suppressMessages(require(pathview))
  suppressMessages(require(gage))

  #--------Output directory--------
  message("Creating Output Directory")
  
  mainDir <- setwd(getwd())
  outDir <- file.path(mainDir)
  subDir <- output_folder
  resDir <- file.path(outDir, subDir)
  
  if (file.exists(resDir)){
    message ("Results folder for analysis found")
  } else {
    dir.create(resDir)
    message ("Results folder for analysis created")
  } 
  
  #--------load data--------
  message("Loading data and annotating Ids")
  
  input_data <- prepare_data(input_file)
  # return(input_data)

  #--------load Kegg pathways--------
  message("Loading Kegg Pathways")
  
  ## pathways from GAGE kegg.gsets
  pathways <- kegg.gsets(species = "hsa",
                         id.type = "kegg",
                         check.new = T)
  # sigmetidx <- pathways$sigmet.idx
  # listkegg <- pathways$kg.sets[sigmetidx] 
  listkegg <- pathways$kg.sets
  
  ## prepare vector with geneid + FC
  foldchanges = input_data$log2FoldChange
  names(foldchanges) = input_data$geneid
  # head(foldchanges)
  
  #--------run Gage-------
  message("Gage Analysis")
  
  keggres = gage(exprs = foldchanges,
                 gsets = listkegg,
                 same.dir = T)
  
  lapply(keggres, head)
  
  keggresp <- as.data.frame(keggres$greater)
  keggresp$pathway <- rownames(keggres$greater)
  keggresp2 <- keggresp %>% dplyr::filter(p.geomean < pgeomean)
  
  keggrespathways = data.frame(id=keggresp2$pathway, keggresp2) %>% 
    tbl_df() %>% 
    # filter(row_number()<=5) %>% 
    .$id %>% 
    as.character()
  
  keggresids = substr(keggrespathways, start=1, stop=8)
  # avoid hsa01100 kegg pathway (too big)
  keggresids <- keggresids[!keggresids %in% "hsa01100"]
  
  #--------Pathview analysis--------
  message("Pathview Analysis")
  
  # detach dplyr (in order to avoid issues) as reported here: https://www.biostars.org/p/242157/
  detach("package:dplyr", unload=TRUE)
  
  # set output
  setwd(output_folder)
  
  # Define plotting function for applying later
  plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
  
  # plot multiple pathways (plots saved to disk and returns a throwaway list object)
  tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))
  
  #---------Saving results----------------
  message("Saving Results on Disk and as Rds object")
  
  write.table(keggresp2, "GageResults.tab", sep = "\t", row.names = F)
  
  report_data <- list()
  
  # save results (for Rmd report)
  report_data[["DataFrame"]] <- as.data.frame(input_data)
  report_data[["FoldChanges"]] <- foldchanges
  report_data[["Pgeomean"]] <- pgeomean
  report_data[["PathviewResults"]] <- tmp
  report_data[["GageResults"]] <- keggresp2
  
  # saveRds (tmp)
  saveRDS(report_data, "GagePathviewAnalysis.Rds")
  message("Gage/Pathview analysis completed")
  
}
