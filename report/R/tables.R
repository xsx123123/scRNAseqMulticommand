# ==============================================================================
# tables.R
# Description: Helper functions for rendering tables in the Quarto report.
# ==============================================================================

# Render a QC summary table from one or more QC summary objects.
render_qc_table <- function(summary_list) {
  if (is.null(summary_list)) return(data.frame())
  if (!is.list(summary_list)) summary_list <- list(summary_list)
  
  rows <- lapply(summary_list, function(s) {
    if (is.null(s$samples) || length(s$samples) == 0) return(NULL)
    samples <- s$samples
    if (is.data.frame(samples)) {
      data.frame(
        Sample = samples$name,
        Origin.Cells = samples$origin_cells,
        Filtered.Cells = samples$filtered_cells,
        Retention.Rate = round(samples$filtered_cells / samples$origin_cells * 100, 2),
        MT.Cutoff = samples$mt_cutoff,
        nFeature.Lower = samples$nfeature_lower,
        nFeature.Upper = samples$nfeature_upper,
        stringsAsFactors = FALSE
      )
    } else {
      do.call(rbind, lapply(samples, function(x) data.frame(
        Sample = x$name,
        Origin.Cells = x$origin_cells %||% NA,
        Filtered.Cells = x$filtered_cells %||% NA,
        Retention.Rate = if (!is.null(x$origin_cells) && !is.null(x$filtered_cells) && x$origin_cells > 0) {
          round(x$filtered_cells / x$origin_cells * 100, 2)
        } else NA,
        MT.Cutoff = x$mt_cutoff %||% NA,
        nFeature.Lower = x$nfeature_lower %||% NA,
        nFeature.Upper = x$nfeature_upper %||% NA,
        stringsAsFactors = FALSE
      )))
    }
  })
  
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# Render a steps status table from the manifest.
render_steps_table <- function(manifest) {
  if (is.null(manifest$steps)) return(data.frame())
  steps <- manifest$steps
  if (is.data.frame(steps)) {
    steps[, c("step", "status", "summary")]
  } else {
    data.frame(
      step = sapply(steps, function(x) x$step),
      status = sapply(steps, function(x) x$status),
      summary = sapply(steps, function(x) x$summary),
      stringsAsFactors = FALSE
    )
  }
}

# Null-coalescing helper
`%||%` <- function(x, y) if (is.null(x)) y else x
