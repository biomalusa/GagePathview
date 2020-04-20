#------gene annotat: biomaRt--------
ens_to_name <- function(ens_names) {
  require("biomaRt")
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                    #host="asia.ensembl.org",
                    host="uswest.ensembl.org",
                    # host="www.ensembl.org",
                    # host="grch37.ensembl.org",
                    path="/biomart/martservice",
                    dataset="hsapiens_gene_ensembl")
  mapping <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'hgnc_symbol', 'entrezgene', 'description'),
                   filters = 'ensembl_gene_id',
                   values = ens_names, 
                   mart = ensembl)
  
  gene_names <- list(ids = ens_names, names = mapping[match(ens_names, mapping$ensembl_gene_id),]$external_gene_name)
  gene_final_names <- ifelse(is.na(as.character(gene_names$names)), gene_names$ids, gene_names$names)
  return(gene_final_names)
}

ens_to_geneid <- function(ens_names) {
  require("biomaRt")
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                    #host="asia.ensembl.org",
                    # host="uswest.ensembl.org",
                    host="www.ensembl.org",
                    # host="grch37.ensembl.org",
                    path="/biomart/martservice",
                    dataset="hsapiens_gene_ensembl")
  mapping <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'hgnc_symbol', 'entrezgene', 'description'),
                   filters = 'ensembl_gene_id',
                   values = ens_names, 
                   mart = ensembl)
  
  gene_names <- list(ids = ens_names, names = mapping[match(ens_names, mapping$ensembl_gene_id),]$entrezgene)
  gene_final_names <- ifelse(is.na(as.character(gene_names$names)), gene_names$ids, gene_names$names)
  return(gene_final_names)
}

ens_to_desc <- function(ens_names) {
  require("biomaRt")
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                    #host="asia.ensembl.org",
                    # host="uswest.ensembl.org",
                    host="www.ensembl.org",
                    # host="grch37.ensembl.org",
                    path="/biomart/martservice",
                    dataset="hsapiens_gene_ensembl")
  mapping <- getBM(attributes = c('ensembl_gene_id', 
                                  'external_gene_name', 
                                  'hgnc_symbol', 
                                  'entrezgene', 
                                  # 'description', 
                                  'wikigene_description'),
                   filters = 'ensembl_gene_id',
                   values = ens_names, 
                   mart = ensembl)
  
  gene_names <- list(ids = ens_names, names = mapping[match(ens_names, mapping$ensembl_gene_id),]$wikigene_description)
  gene_final_names <- ifelse(is.na(as.character(gene_names$names)), gene_names$ids, gene_names$names)
  return(gene_final_names)
}

name_to_geneid <- function(ens_names) {
  require("biomaRt")
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                    host="www.ensembl.org",
                    path="/biomart/martservice",
                    dataset="hsapiens_gene_ensembl")
  mapping <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'hgnc_symbol', 'entrezgene', 'description'),
                   filters = 'external_gene_name',
                   values = ens_names, 
                   mart = ensembl)
  
  gene_names <- list(ids = ens_names, names = mapping[match(ens_names, mapping$external_gene_name),]$entrezgene)
  gene_final_names <- ifelse(is.na(as.character(gene_names$names)), gene_names$ids, gene_names$names)
  return(gene_final_names)
}