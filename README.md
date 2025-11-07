# A-Multi-omic-Atlas-of-Human-Choroid-Plexus-in-Alzheimer-s-Disease

#Complete snRNA-seq (RNA), snATAC-seq (ATAC), and spatial transcriptomics (ST) analysis pipelines.

util.R suffix include required functions, RNA, ATAC, and ST suffixes indicate dataset.

010 Load_preprocess_util.R - Load, preprocess, and QC individual snRNA/ATAC-seq data samples (split by batch and lanes).

011 Load_preprocess_SN_RNA_ATAC.R


020 Cell_type_annotation_util.R - Cell type annotation and cleaning of each dataset. ...(4 and 5) 

021 Cell_type_annotation_RNA.R

022 Cell_type_annotation_ATAC.R

023 Cell_type_annotation_ST.R

030 limma_pattern_fgsea_util.R - limma-Voom analyses, fgsea GO term analyses, and pattern matching analyses.

031 Cell_type_proportion_RNA.R

032 Cell_type_annotation_differential_RNA.R

033 AD_differential_RNA.R


034 Cell_type_proportion_ATAC.R

035 Cell_type_annotation_differential_ATAC.R

036 cogdx_AD_differential_ATAC.R


037 Cell_type_proportion_ST.R

038 Cell_type_annotation_differential_ST.R

039 cogdx_AD_differential_ST.R


040 CellChat_util_RNA.R

041 CellChat_Major_RNA.R


050 ST_specific_plots_util.R - Spatial molecular expression visualizations

051 Major_cell_type_plots.R

052 Sub_cell_type_plots.R

053 AD_differential_plots.R


060 Spatial_proximity_util.R - Distance and proximity statistics

061 Spatial_proximity.R


070 Spatial_neighborhood_util.R - Neighbourhood-based spatial clustering

071 Spatial_neighborhood.R

072 Spatial_neighborhood_proportions.R

073 Spatial_neighborhood_differential.R

080 Spatial_cellchat_util.R - CellChat based on Neighbourhood spatial clustering
081 Spatial_cellchat.R
