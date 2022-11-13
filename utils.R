process_obj <- function(obj, seed = 1234, nfeat = 800, assay = "RNA", ndim=10, find_variable_features = T){
  obj <- NormalizeData(obj, verbose = TRUE)
  if(find_variable_features){
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeat , verbose = FALSE)
  }
  obj <- ScaleData(obj, verbose = TRUE)
  obj <- RunPCA(obj, features = (obj@assays[[assay]])@var.features,
                ndims.print = 1:5, nfeatures.print = 5)
  
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:ndim, seed.use=seed, n.neighbors = 30, min.dist=0.3)
  return(obj)
}

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, 
                   attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

  #humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  #print(head(humanx))
  #return(humanx)
  return(genesV2)
}

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

hs.2.mm <- function(hs,keep.na =F){
  if(!keep.na){
    mm <- subset(ProjecTILs::Hs2Mm.convert.table,Gene.HS %in% hs)$Gene.MM
  }else{
    tt <- ProjecTILs::Hs2Mm.convert.table
    ii <- match(hs,tt$Gene.HS)
    mm <- (tt[ii,])$Gene.MM
  }
  return(mm)
}