# ==============================================================================
# plot_helpers.R
# Description: Common plotting helpers for the Quarto report.
# ==============================================================================

# Display an image file if it exists, otherwise print a placeholder message.
display_image <- function(result_dir, rel_path, caption = NULL) {
  img_path <- file.path(result_dir, rel_path)
  if (!file.exists(img_path)) {
    message("Image not found: ", rel_path)
    return(invisible(NULL))
  }
  knitr::include_graphics(img_path)
}

# Find the first artifact of a given type from a summary object.
find_artifact <- function(summary_obj, type) {
  if (is.null(summary_obj) || is.null(summary_obj$artifacts)) return(NULL)
  artifacts <- summary_obj$artifacts
  if (is.data.frame(artifacts)) {
    match <- artifacts[artifacts$type == type, ]
    if (nrow(match) > 0) return(match$path[1])
  } else if (is.list(artifacts)) {
    for (art in artifacts) {
      if (!is.null(art$type) && art$type == type) return(art$path)
    }
  }
  NULL
}
