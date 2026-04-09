###
### 011 Load_preprocess_SN_ATAC.R
###
# Purpose: Load, preprocess, and QC individual snRNA/ATAC-seq data samples (split by batch and lanes).
# Dependencies:
source("/ChP_RNA/010 Load_preprocess_util.R", verbose = FALSE)
curr_dir <- "/ChP_RNA" #output results to specific folder

Metadata <- read_csv(paste0(curr_dir, "/metadata/ex vivo ChP metadata.csv")) #Participant metadata
set.seed(1234)

#Load
##Preprocess data
###Batch1-10 and newCh 12-21 (example usage)
data_dir <- "/data_dir/Cellranger/" #Specifies the data's directory !!!!Modify!!!!
loadfun('ChP_1_1', data_dir, exp_cmo=1:8)
loadfun('ChP_1_2', data_dir, exp_cmo=1:8)
loadfun('ChP_1_3', data_dir, exp_cmo=1:8)
loadfun('ChP_1_4', data_dir, exp_cmo=1:8)
ChP_1 <- merge(ChP_1_1, y=c(ChP_1_2, ChP_1_3, ChP_1_4), add.cell.ids=c("ChP1_1","ChP1_2","ChP1_3","ChP1_4"))

demultiplexfun(ChP_1, assay="CMO", p.quant_CMO1=.95, p.quant_CMO2=.99999)

demultiplexfunm <- function(srt){
  nameobj <- deparse(substitute(srt)) #make char
  Idents(srt) <- "CMO_maxID"
  
  #Rename by participant ID; note that projIDs are masked
  srt <- RenameIdents(object = srt, 'CMO-1'='109', 'CMO-2'='109', 'CMO-3'='153', 'CMO-4'='153', 'CMO-5'='104', 'CMO-6'='104', 'CMO-7'='133', 'CMO-8'='133')
  srt@meta.data$sample <- srt@active.ident
  assign(paste(nameobj), srt, envir = .GlobalEnv)  #to global environment
}
demultiplexfunm(ChP_1)

metafun(ChP_1, "109")
metafun(ChP_1, "153")
metafun(ChP_1, "133")
metafun(ChP_1, "104")

rm(list = grep("ChP", ls(), value=TRUE))
rm(list = grep("s[[:digit:]]", ls(), value=TRUE))
gc()


###Batch 11 (example usage note metafunB and option to load ATAC data)
data_dir <- "/bucket_data/FC103_FC103_HT_Multiome/"
ATAC_data_dir <- "/bucket_data/FC103_FC104_arc_multiome/"
loadfun('C1_GEX', data_dir, ATAC_data_dir = ATAC_data_dir, exp_cmo=1:12)
loadfun('C2_GEX', data_dir, ATAC_data_dir = ATAC_data_dir, exp_cmo=1:12)
loadfun('C3_GEX', data_dir, ATAC_data_dir = ATAC_data_dir, exp_cmo=1:12)
loadfun('C4_GEX', data_dir, ATAC_data_dir = ATAC_data_dir, exp_cmo=1:12)
ChP_11 <- merge(sC1_GEX, y=c(sC2_GEX, sC3_GEX, sC4_GEX), add.cell.ids=c("ChP11_1","ChP11_2","ChP11_3","ChP11_4"))

#demultiplex
demultiplexfun(ChP_11, assay="CMO", p.quant_CMO1=.95, p.quant_CMO2=.99999)

demultiplexfunm <- function(srt){
  nameobj <- deparse(substitute(srt)) #make char
  #CMO metadata !!!!MODIFY!!!!
  Idents(srt) <- "CMO_maxID"
  srt <- RenameIdents(object = srt, 'CMO-1'='167', 'CMO-2'='167', 'CMO-3'='131', 'CMO-4'='131', 'CMO-5'='145', 'CMO-6'='145', 'CMO-7'='163', 'CMO-8'='163', 'CMO-9'='145', 'CMO-10'='145', 'CMO-11'='141', 'CMO-12'='141')
  srt@meta.data$sample <- srt@active.ident
  assign(paste(nameobj), srt, envir = .GlobalEnv)  #to global environment
}
demultiplexfunm(ChP_11)

#Split, add metadata, and save
metafunB(ChP_11, "167")
metafunB(ChP_11, "131")
metafunB(ChP_11, "145")
metafunB(ChP_11, "163")
metafunB(ChP_11, "145")
metafunB(ChP_11, "141")

rm(list = grep("ChP", ls(), value=TRUE))
rm(list = grep("sC[[:digit:]]", ls(), value=TRUE))
gc()