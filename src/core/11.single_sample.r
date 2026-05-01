# ==============================================================================
# Single Sample scRNA-seq Analysis Pipeline
# Author: zhang jian
# Description: Main workflow for processing a single scRNA-seq sample.
#              Loads data -> QC -> Normalization -> DimRed -> Clustering -> Annotation.
# ==============================================================================

singlesample_scRNA_seq_analysis <- function(scRNAtype, name, reduceType = 'FALSE', ctx, normalization_method = "LogNormalize"){
  
  # Unpack Context: Load all variables from ctx into the current environment
  list2env(ctx, envir = environment())
  
  info(logger, paste0('>>> Starting Single Sample Analysis: ', name))
  
  # 1. Load Data & QC
  # Using the optimized DNBC4_10X wrapper (handles both formats)
  # Note: DNBC4_10X returns a LIST. For single sample, we take the first element.
  # Or we use PatchDealform10X/DNB directly if the input logic allows.
  # Let's assume 'DNBC4_10X' from 02.cells_qc.r is available and handles list input.
  
  # Construct a temporary list for compatibility with list-based loader
  input_list <- list()
  input_list[[name]] <- Cellranger_dir
  
  info(logger, 'Loading and performing QC...')

  group_data <- cellRangerlist_dataframe |> dplyr::select(c(name,group))

  # This returns a list of Seurat objects
  Seurat_list <- DNBC4_10X(scRNAtype = scRNAtype,scRNAAutofilted = scRNAAutofilted,
                           input_list,group_data,Cellranger_dir,origin_tax_ID,
                           cellRangerlist_dataframe = cellRangerlist_dataframe)
  
  Seurat_matrix <- Seurat_list[[1]]
  
  # 2. Normalization & Dimensionality Reduction
  info(logger, paste0('Normalization & Feature Selection (Method: ', normalization_method, ')'))
  
  # Using optimized NormalFeature wrapper
  # Returns a list, extract first
  Seurat_list_norm <- NormalFeature(list(Seurat_matrix), method = normalization_method, vars.to.regress = "percent.mt")
  Seurat_matrix <- Seurat_list_norm[[1]]
  
  # Run PCA (Wrapper handles ScaleData if needed)
  Seurat_matrix <- RunPCA_Wrapper(Seurat_matrix)
  
  # 3. Determine Dimensions (ElbowPlot)
  pc_cutoff <- DrawElbowPlot(Seurat_matrix, name, figure_dir)
  
  # 4. Clustering & Embedding
  info(logger, 'Clustering & Dimensionality Reduction...')
  Seurat_matrix <- FindNeighbors(Seurat_matrix, dims = 1:pc_cutoff)
  
  # Multi-res clustering
  res_list <- c(0.4, 0.6, 0.8, 1.0, 1.2)
  for(r in res_list) Seurat_matrix <- FindClusters(Seurat_matrix, resolution = r, verbose = FALSE)
  Idents(Seurat_matrix) <- paste0("RNA_snn_res.", 1.2) # Default
  if(normalization_method == "SCT") Idents(Seurat_matrix) <- paste0("SCT_snn_res.", 1.2)
  
  # UMAP
  Seurat_matrix <- RunUMAP(Seurat_matrix, dims = 1:pc_cutoff)
  
  # tSNE (Optional)
  if(reduceType == 'TRUE'){
    Seurat_matrix <- RunTSNE(Seurat_matrix, dims = 1:pc_cutoff)
  }
  
  # Save basic processed object
  saveRDS(Seurat_matrix, file.path(output_dir, 'scrna_seq.rds'))
  
  # 5. Marker Genes
  info(logger, 'Finding Marker Genes...')
  p <- FindClusterMarkersDotplot(Seurat_matrix, group_by = 'seurat_clusters',
                                 save_dir = cluster_marker_gene_dir,
                                 draw_plot = TRUE,
                                 topgene = 10)
  
  # 6. Doublet Detection
  info(logger, 'Running DoubletFinder...')
  Seurat_list_dbl <- Checkdoublet(list(Seurat_matrix), num.cores = NULL,
                                  sct_doubletFinder = FALSE,
                                  save_dir = doublet_dir)
  Seurat_matrix <- Seurat_list_dbl[[1]]
  saveRDS(Seurat_matrix, file.path(doublet_dir, "doublet-Seurat.rds"))
  
  # 7. Ambient RNA Removal
  info(logger, 'Estimating Ambient RNA Contamination...')
  Seurat_matrix <- AmbientRNAContamination(Seurat_matrix, RNAContamination_dir)
  DrawAmbientRNAContamination_UMAP(Seurat_matrix, RNAContamination_dir, 'NULL')
  
  # 8. Visualization
  draw_umap_tsne_plot(Seurat_matrix, reduceType = reduceType, intergetmethods = 'NULL', name = name, figure_dir = figure_dir)
  
  # 9. Automated Annotation (Species Specific)
  # Only run if Human/Mouse, or Custom DB provided
  if(origin_tax_ID %in% c(9606, 10090)){
    info(logger, 'Running Automated Annotation (SingleR)...')
    
    # Cell Level
    Seurat_matrix <- RunSingleR_Unified(Seurat_matrix, origin_tax_ID, annotation_SinglR_dir, level = "cell")
    # Cluster Level
    Seurat_matrix <- RunSingleR_Unified(Seurat_matrix, origin_tax_ID, annotation_SinglR_dir, level = "cluster")
    
    # Visualize SingleR
    # (Assuming visualization functions are updated to handle generic SingleR cols)
  } else {
    info(logger, paste0("Skipping SingleR Annotation for TaxID: ", origin_tax_ID, " (Plant/Other)."))
  }
  
  # CellID / ScType (Optional)
  # RunScType(Seurat_matrix, ...) 
  # info(logger,crayon::blue(crayon::bold('>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<')))
  # info(logger, crayon::bold(crayon::inverse(">>> STEP 14 : Seurat format to Annodata")))
  # Seurat2h5ad(Seurat_obj = Seurat_matrix,
  #            project_name = ctx$project_name,
  #            save_dir = ctx$output_dir,assay_to_use = "RNA",
  #            logger = logger)
  #
  # info(logger, '>>> Single Sample Analysis Complete.')
}
# ==============================================================================
# END 
# ==============================================================================
