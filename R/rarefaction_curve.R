#' @title Generate Rarefaction Curves
#' @description This function generates rarefaction curves for visualizing species richness
#'              at different sequencing depths, based on the OTU table in a phyloseq object.
#' @param data A phyloseq object containing OTU data and associated metadata.
#' @param step An integer specifying the step size for the rarefaction curve. Default is 100.
#' @param ylim A numeric vector of length two specifying the y-axis limits. Default is `NULL`,
#'             which automatically adjusts the y-axis.
#' @param xlim A numeric vector of length two specifying the x-axis limits. Default is `NULL`,
#'             which automatically adjusts the x-axis.
#' @param label A logical value indicating whether sample labels should be displayed on the
#'              plot. Default is `FALSE`.
#' @param color A character string specifying the column in the sample metadata to use for
#'              coloring the samples. Default is `NULL` (no color applied).
#' @param rank An optional rank level for taxa aggregation. Default is `rank`.
#' @param legend A logical value indicating whether to display a legend on the plot.
#'               Default is `TRUE`.
#' @param color.vector An optional character vector specifying custom colors for groups in
#'                     the `color` column. Default is `NULL` (default colors are used).
#' @param legend.position A character string specifying the position of the legend. Default
#'                        is `"topleft"`.
#' @return A rarefaction curve plot.
#' @examples
#' amp_rarecurve(data = physeq, step = 200, color = "sample_type", legend.position = "topright")
#' @export
# Functie voor rarefaction curve (gebaseerd op amp_rarecurve)
amp_rarecurve <- function(data, step = 100, ylim = NULL, xlim = NULL, label = FALSE, color = NULL, rank = rank, legend = TRUE, color.vector = NULL, legend.position = "topleft") {
  abund = otu_table(data)@.Data %>% as.data.frame()
  if (!is.null(color)) {
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    group_vector <- sample_data(data)[, color]@.Data %>% as.data.frame()
    names(group_vector) <- "color_variable"
    group_vector <- as.character(group_vector$color_variable)
    groups <- unique(group_vector)
    n = length(groups)
    cols = gg_color_hue(n)
    if (!is.null(color.vector)) {
      cols <- color.vector
    }
    col_vector <- rep("black", length(group_vector))
    for (i in seq_along(group_vector)) {
      col_vector[i] <- cols[match(group_vector[i], groups)]
    }
  } else {
    col_vector = "black"
  }
  if (is.null(ylim) & is.null(xlim)) {
    vegan::rarecurve(t(abund), rank = rank, step = step, label = label, col = col_vector)
  }
  if (!is.null(ylim) & !is.null(xlim)) {
    vegan::rarecurve(t(abund), rank = rank, step = step, ylim = ylim, xlim = xlim, label = label, col = col_vector)
  }
  if (!is.null(ylim) & is.null(xlim)) {
    vegan::rarecurve(t(abund), rank = rank, step = step, ylim = ylim, label = label, col = col_vector)
  }
  if (is.null(ylim) & !is.null(xlim)) {
    vegan::rarecurve(t(abund), rank = rank, step = step, xlim = xlim, label = label, col = col_vector)
  }
  if (!is.null(color) & legend) {
    legend(legend.position, legend = groups, fill = cols, bty = "n")
  }
}

#' @title Generate and Save Rarefaction Curve
#' @description This function creates a rarefaction curve for a given phyloseq object and
#'              saves the plot as a PDF.
#' @param physeq A phyloseq object containing OTU data and associated metadata.
#' @param rank The rank level at which taxa should be aggregated. Default is `rank`.
#' @param processing_stage A character string describing the processing stage (e.g., "raw" or "filtered").
#'                         This string will be included in the saved PDF file name.
#' @param color A character string specifying the column in the sample metadata to use for
#'              coloring the samples. Default is `NULL`, which automatically sets the color
#'              to `"sample_type"`.
#' @param base_path A character string specifying the base directory for saving output files.
#'                  Default is `"./"`.
#' @return The rarefaction curve plot object.
#' @details This function first checks whether the `sample_type` column exists in the
#'          sample metadata. It then generates a rarefaction curve using the `amp_rarecurve`
#'          function and saves the plot as a PDF file in the specified directory.
#' @examples
#' rarefaction_curve(physeq = physeq, rank = "Genus", processing_stage = "filtered", color = "Treatment")
#' @export

# Functie voor het genereren van een rarefaction curve
rarefaction_curve <- function(physeq = resolved_tree_physeq, rank = "ASV", processing_stage, color = NULL, base_path = "./") {
  # Controleer of sample_type aanwezig is in sample_data
  if (!("sample_type" %in% colnames(sample_data(physeq)))) {
    stop("error: 'sample_type' column is missing from the sample data.")
  }

  # Stel default kleur in als sample_type
  if (is.null(color)) {
    color = "sample_type"
  } else {
    message(paste("message: using", color, "for color in rarefaction plot."))
  }

  # Stel project en folderstructuur in
  project_name = projects
  project_folder = paste0(base_path, project_name)
  figure_folder = paste0(project_folder, "/figures/")
  output_folder_csv_files = paste0(project_folder, "/output_data/csv_files/")
  output_folder_rds_files = paste0(project_folder, "/output_data/rds_files/")

  # Maak de rarefaction plot en sla deze op als PDF
  pdf(glue("{figure_folder}{project_name}_rarefaction_curve_{processing_stage}.pdf"),
      useDingbats = FALSE, width = 7, height = 5)

  rarefaction_plot = amp_rarecurve(physeq, rank = rank, color = color,
                                   legend.position = "bottomright", xlim = c(0, max(sample_sums(physeq))))

  print(rarefaction_plot)
  dev.off()

  return(rarefaction_plot)
}
