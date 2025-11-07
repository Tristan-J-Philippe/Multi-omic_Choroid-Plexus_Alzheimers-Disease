# A-Multi-omic-Atlas-of-Human-Choroid-Plexus-in-Alzheimer-s-Disease

#Complete snRNA-seq (RNA), snATAC-seq (ATAC), and spatial transcriptomics (ST) analysis pipelines.<br />
<br />
util.R suffix include required functions, RNA, ATAC, and ST suffixes indicate dataset.<br />
010 Load_preprocess_util.R - Load, preprocess, and QC individual snRNA/ATAC-seq data samples (split by batch and lanes).<br />
011 Load_preprocess_SN_RNA_ATAC.R<br />
<br />
020 Cell_type_annotation_util.R - Cell type annotation and cleaning of each dataset. ...(4 and 5) <br />
021 Cell_type_annotation_RNA.R<br />
022 Cell_type_annotation_ATAC.R<br />
023 Cell_type_annotation_ST.R<br />
030 limma_pattern_fgsea_util.R - limma-Voom analyses, fgsea GO term analyses, and pattern matching analyses.<br />
031 Cell_type_proportion_RNA.R<br />
032 Cell_type_annotation_differential_RNA.R<br />
033 AD_differential_RNA.R<br />
<br />
034 Cell_type_proportion_ATAC.R<br />
035 Cell_type_annotation_differential_ATAC.R<br />
036 cogdx_AD_differential_ATAC.R<br />
<br />
037 Cell_type_proportion_ST.R<br />
038 Cell_type_annotation_differential_ST.R<br />
039 cogdx_AD_differential_ST.R<br />
<br />
040 CellChat_util_RNA.R<br />
041 CellChat_Major_RNA.R<br />
<br />
050 ST_specific_plots_util.R - Spatial molecular expression visualizations<br />
051 Major_cell_type_plots.R<br />
052 Sub_cell_type_plots.R<br />
053 AD_differential_plots.R<br />
<br />
060 Spatial_proximity_util.R - Distance and proximity statistics<br />
061 Spatial_proximity.R<br />
<br />
070 Spatial_neighborhood_util.R - Neighbourhood-based spatial clustering<br />
071 Spatial_neighborhood.R<br />
072 Spatial_neighborhood_proportions.R<br />
073 Spatial_neighborhood_differential.R<br />
<br />
080 Spatial_cellchat_util.R - CellChat based on Neighbourhood spatial clustering<br />
081 Spatial_cellchat.R<br />
