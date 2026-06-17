# ==============================================================================
# Report Manifest Generator
# Author: zhang jian
# Description: Generates standardized summary.json for each analysis step and a
#              global manifest.json for downstream Quarto report consumption.
#              Keeps the pipeline decoupled from Quarto rendering logic.
# ==============================================================================

# ==============================================================================
# Helper: .safe_json_write
# Description: Safely writes an R list to JSON file with pretty formatting.
# ==============================================================================
.safe_json_write <- function(data, file_path) {
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(data, path = file_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
}

# ==============================================================================
# Function: write_step_summary
# Description: Writes a per-step summary.json file.
# @ step       : Step identifier (e.g., "qc", "doublet", "ambient_rna").
# @ data       : Named list containing step-specific data (samples, metrics, artifacts).
# @ output_dir : Directory where the step summary should be written.
# @ status     : "success", "warning", or "error".
# @ timestamp  : Optional ISO timestamp. Defaults to current time.
# ==============================================================================
write_step_summary <- function(step, data, output_dir, status = "success", timestamp = NULL) {
  if (is.null(timestamp)) {
    timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  }
  
  summary_data <- list(
    step = step,
    status = status,
    timestamp = timestamp,
    samples = data$samples %||% list(),
    metrics = data$metrics %||% list(),
    artifacts = data$artifacts %||% list()
  )
  
  summary_path <- file.path(output_dir, "summary.json")
  .safe_json_write(summary_data, summary_path)
  
  if (exists("logger")) {
    info(logger, paste0("  Wrote step summary: ", summary_path))
  }
  
  invisible(summary_path)
}

# ==============================================================================
# Function: write_manifest
# Description: Writes the global manifest.json based on Run_Context and discovered
#              step summaries.
# @ ctx        : Pipeline Run_Context list.
# @ output_dir : Root output directory of the pipeline.
# @ steps      : Data frame or list of steps with columns step, summary_path, status.
# @ final_objects : Optional list of final output objects.
# ==============================================================================
write_manifest <- function(ctx, output_dir, steps = NULL, final_objects = NULL, timestamp = NULL) {
  if (is.null(timestamp)) {
    timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  }
  
  # Infer project metadata from ctx
  project_name <- ctx$project_name %||% "unknown"
  species_tax_id <- ctx$origin_tax_ID %||% NA_integer_
  samples <- ctx$infor$infor %||% character(0)
  groups <- unique(ctx$infor$group %||% character(0))
  integration_method <- ctx$intergetmethods %||% "NULL"
  pipeline_version <- if (exists("ScriptsVersion")) ScriptsVersion else "unknown"
  
  # Auto-discover step summaries if not provided
  if (is.null(steps)) {
    steps <- .discover_steps(output_dir)
  }
  
  manifest <- list(
    pipeline_version = pipeline_version,
    project_name = project_name,
    run_timestamp = timestamp,
    species_tax_id = species_tax_id,
    samples = as.character(samples),
    groups = as.character(groups),
    integration_method = as.character(integration_method),
    steps = steps,
    final_objects = final_objects %||% list()
  )
  
  manifest_path <- file.path(output_dir, "manifest.json")
  .safe_json_write(manifest, manifest_path)
  
  if (exists("logger")) {
    info(logger, paste0("  Wrote global manifest: ", manifest_path))
  }
  
  invisible(manifest_path)
}

# ==============================================================================
# Helper: .discover_steps
# Description: Scans the output directory for existing summary.json files and
#              returns a list of steps for the manifest.
# ==============================================================================
.discover_steps <- function(output_dir) {
  step_dirs <- list.dirs(output_dir, recursive = TRUE)
  summary_paths <- step_dirs[file.exists(file.path(step_dirs, "summary.json"))]
  
  if (length(summary_paths) == 0) {
    return(list())
  }
  
  steps <- lapply(summary_paths, function(dir_path) {
    rel_path <- sub(paste0("^", normalizePath(output_dir), "/?"), "", normalizePath(dir_path))
    rel_path <- ifelse(rel_path == "", ".", rel_path)
    summary_json <- tryCatch(
      jsonlite::read_json(file.path(dir_path, "summary.json")),
      error = function(e) list(step = basename(dir_path), status = "unknown")
    )
    list(
      step = summary_json$step %||% basename(dir_path),
      summary = file.path(rel_path, "summary.json"),
      status = summary_json$status %||% "unknown"
    )
  })
  
  steps
}

# ==============================================================================
# Function: build_qc_summary
# Description: Builds QC summary data from existing QC outputs.
# @ qc_dir     : QC output directory.
# @ sample_info : Data frame with columns name, group, library_type.
# ==============================================================================
build_qc_summary <- function(qc_dir, sample_info) {
  cellranger_dir <- file.path(qc_dir, "Cellranger-result")
  if (!dir.exists(cellranger_dir)) {
    return(list(samples = list(), metrics = list(), artifacts = list()))
  }
  
  setting_files <- list.files(cellranger_dir, pattern = "_filted_setting\\.csv$", full.names = TRUE)
  
  samples <- lapply(setting_files, function(f) {
    sample_name <- sub("_filted_setting\\.csv$", "", basename(f))
    settings <- read.csv(f, stringsAsFactors = FALSE)
    
    list(
      name = sample_name,
      origin_cells = NA_integer_,
      filtered_cells = NA_integer_,
      mt_cutoff = settings$MT_cutoff_upper[1],
      nfeature_lower = settings$nFeature_cutoff_lower[1],
      nfeature_upper = settings$nFeature_cutoff_upper[1]
    )
  })
  
  artifacts <- list(
    list(type = "dir", path = "QC/Cellranger-result")
  )
  
  list(samples = samples, metrics = list(), artifacts = artifacts)
}

# ==============================================================================
# Function: build_final_objects
# Description: Builds final_objects list for manifest from common output files.
# @ output_dir : Root output directory.
# ==============================================================================
build_final_objects <- function(output_dir) {
  objects <- list()
  
  merge_obj <- file.path(output_dir, "output", "scrna_seq_merge.qs")
  if (file.exists(merge_obj)) {
    objects <- c(objects, list(list(
      type = "qs",
      path = "output/scrna_seq_merge.qs",
      description = "Final merged Seurat object"
    )))
  }
  
  cell_info <- file.path(output_dir, "output", "Cell-cluster-infor.csv")
  if (file.exists(cell_info)) {
    objects <- c(objects, list(list(
      type = "csv",
      path = "output/Cell-cluster-infor.csv",
      description = "Cell metadata and cluster assignments"
    )))
  }
  
  objects
}

# ==============================================================================
# Helper: %||%
# Description: Null-coalescing operator for concise default values.
# ==============================================================================
`%||%` <- function(x, y) if (is.null(x)) y else x
# ==============================================================================
# END
# ==============================================================================
