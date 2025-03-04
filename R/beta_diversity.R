#' Generate Beta Diversity Ordination Plots
#'
#' This function computes beta diversity using ordination methods (e.g., PCoA) on a phyloseq object and generates corresponding ordination plots. It supports both ASV-level data and data aggregated at various taxonomic levels, and it can handle both relative and, if available, absolute abundance data. Four distance metrics are used for relative abundance plots (Jaccard, Bray-Curtis, Unweighted UniFrac, and Weighted UniFrac), and Manhattan distance is used for absolute abundance plots when a normalization method is provided.
#'
#' @param physeq A phyloseq object containing microbial community data.
#' @param taxrank A character vector specifying the taxonomic levels to process. If the first element (case-insensitive) is \code{"asv"}, ASV-level beta diversity is computed; otherwise, beta diversity is computed for each taxonomic level provided (default: \code{c("Phylum", "Class", "Order", "Family", "Tax_label")}).
#' @param norm_method A character string indicating the normalization method to be used for generating absolute abundance data. Options include:
#'   \itemize{
#'     \item \code{"fcm"}: Flow cytometry normalization.
#'     \item \code{"qpcr"}: qPCR normalization.
#'     \item \code{NULL}: No absolute abundance processing (default).
#'   }
#' @param ordination_method A character string specifying the ordination method to use (default is \code{"PCoA"}).
#' @param color_factor An optional character string specifying the sample metadata column used to color the points in the ordination plot.
#' @param color_continuous A logical value indicating whether the color scale should be continuous (\code{TRUE}) or discrete (\code{FALSE}). Default is \code{TRUE}.
#' @param shape_factor An optional character string specifying the sample metadata column used to assign point shapes in the ordination plot.
#' @param size_factor An optional character string specifying the sample metadata column used to assign point sizes in the ordination plot.
#' @param alpha_factor An optional character string specifying the sample metadata column used to assign point transparency in the ordination plot.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Sets up project directories and folder paths for saving output files.
#'   \item Defines an internal helper function, \code{base_beta_plot}, which performs the ordination (using the specified \code{ordination_method} and distance metric), removes default point layers, and then adds a customized geom_point layer with aesthetics defined by the provided factors.
#'   \item For ASV-level data (when \code{taxrank[1]} is \code{"asv"}):
#'     \itemize{
#'       \item Processes relative abundance data (transformed to percentages) and, if a normalization method is specified, also processes absolute abundance data.
#'       \item Generates ordination plots using multiple distance metrics:
#'         \itemize{
#'           \item \strong{Jaccard}: Based on binary presence/absence.
#'           \item \strong{Bray-Curtis}: Incorporates both presence and abundance.
#'           \item \strong{Unweighted UniFrac}: Considers lineage presence only.
#'           \item \strong{Weighted UniFrac}: Considers both lineage presence and abundance.
#'         }
#'       \item If both DNA and RNA data are present, separate plots are generated for each.
#'     }
#'   \item For other taxonomic levels:
#'     \itemize{
#'       \item Similar processing is performed for each taxonomic rank, with output saved in dedicated subfolders.
#'     }
#'   \item Aesthetic scales (colors, shapes, sizes, alpha) are defined based on the unique levels in the sample data.
#' }
#'
#' @return A ggplot object representing the combined beta diversity ordination plot for the relative (and, if applicable, absolute) data.
#'
#' @examples
#' \dontrun{
#'   # Example: Generate beta diversity plots at the ASV level using PCoA ordination and flow cytometry normalization
#'   beta_plot <- beta_diversity(
#'     physeq = my_physeq,
#'     taxrank = "asv",
#'     norm_method = "fcm",
#'     ordination_method = "PCoA",
#'     color_factor = "Treatment",
#'     shape_factor = "Replica",
#'     size_factor = "Timepoint"
#'   )
#'
#'   # Example: Generate beta diversity plots for Phylum and Class levels without absolute data processing
#'   beta_plot <- beta_diversity(
#'     physeq = my_physeq,
#'     taxrank = c("Phylum", "Class"),
#'     norm_method = NULL,
#'     ordination_method = "PCoA",
#'     color_factor = "Soil_Type"
#'   )
#' }
#'
#' @export
beta_diversity <- function(physeq = physeq,
                           taxrank = c("Phylum", "Class", "Order", "Family", "Tax_label"),
                           norm_method = NULL,
                           ordination_method = "PCoA",
                           color_factor = NULL,
                           color_continuous = TRUE,
                           shape_factor = NULL,
                           size_factor = NULL,
                           alpha_factor = NULL) {

  # Set up project and folder paths
  project_name <- projects
  project_folder <- paste0(base_path, project_name)
  figure_folder <- paste0(project_folder, "/figures/")
  destination_folder <- paste0(project_folder, "/input_data/")
  output_folder_csv_files <- paste0(project_folder, "/output_data/csv_files/")
  output_folder_rds_files <- paste0(project_folder, "/output_data/rds_files/")

  # Function for creating the base beta-diversity plot
  base_beta_plot <- function(psdata, ordination_method, distance_method, title,
                             color_factor, shape_factor, size_factor, alpha_factor) {
    # Perform the ordination
    ordination_res <- ordinate(psdata, method = ordination_method, distance = distance_method)
    base_plot <- plot_ordination(psdata, ordination = ordination_res, axes = c(1, 2))

    # Remove any existing point layers to add our own
    if (length(base_plot$layers) > 0) {
      base_plot$layers <- base_plot$layers[!sapply(base_plot$layers, function(x) inherits(x$geom, "GeomPoint"))]
    }

    # Build the aesthetic mapping list
    aes_params <- list()
    if (!is.null(color_factor)) aes_params$color <- sym(color_factor)
    if (!is.null(shape_factor)) aes_params$shape <- sym(shape_factor)
    if (!is.null(size_factor))  aes_params$size  <- sym(size_factor)
    if (!is.null(alpha_factor)) aes_params$alpha <- sym(alpha_factor)

    # Add points with the specified aesthetics
    base_plot <- base_plot + geom_point(mapping = do.call(aes, aes_params))

    # Apply color scales based on whether the color factor is continuous or discrete
    if (!is.null(color_factor)) {
      if (color_continuous == TRUE) {
        base_plot <- base_plot + scale_color_continuous(low = "lightblue", high = "darkgreen")
      } else if (color_continuous == FALSE) {
        base_plot <- base_plot + scale_color_manual(values = colorset)
      }
    }
    if (!is.null(shape_factor)) {
      base_plot <- base_plot + scale_shape_manual(values = shapeset)
    }
    if (!is.null(size_factor)) {
      base_plot <- base_plot + scale_size_manual(values = sizeset)
    }
    if (!is.null(alpha_factor)) {
      base_plot <- base_plot + scale_alpha_continuous(range = c(0.3, 1))  # Set transparency for alpha
    }

    # Add title and labels; adjust theme and axis settings
    base_plot <- base_plot +
      ggtitle(title) +
      labs(color = color_factor, shape = shape_factor, size = size_factor, alpha = alpha_factor) +
      theme(panel.background = element_rect(fill = "transparent"),
            panel.grid = element_line(colour = "grey90"),
            strip.text = element_text(face = "bold"),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
            axis.line.y = element_line(color = "black", linewidth = 0.3),
            axis.line.x = element_line(color = "black", linewidth = 0.3)) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))  # 5% expansion for y-axis

    return(base_plot)
  }

  # Create the main beta-diversity folder
  beta_div_folder <- paste0(figure_folder, "Beta_diversity/")
  if (!dir.exists(beta_div_folder)){dir.create(beta_div_folder)}

  if (tolower(taxrank[1]) == "asv") {
    log_message("Processing ASV-level beta diversity", log_file)

    if (is.null(norm_method)) {
      psdata_relative <- physeq
    } else if (norm_method == "fcm") {
      psdata_relative <- physeq[["psdata_asv_copy_number_corrected"]]
      psdata_absolute <- physeq[["psdata_asv_fcm_norm_rarefied"]]
    } else if (norm_method == "qpcr") {
      psdata_relative <- physeq[["psdata_asv_copy_number_corrected"]]
      psdata_absolute <- physeq[["psdata_asv_qpcr_norm_rarefied"]]
    }

    # Create the ASV folder under beta-diversity
    asv_folder <- paste0(beta_div_folder, "ASV/")
    if (!dir.exists(asv_folder)){dir.create(asv_folder)}

    # Transform counts to relative abundance (percentage)
    psdata_relative <- transform_sample_counts(psdata_relative, function(x) x / sum(x) * 100)

    # Convert specified variables to factors for absolute data (if available)
    if (!is.null(norm_method)) {
      if (!is.null(color_factor))
        sample_data(psdata_absolute)[[color_factor]] <- as.factor(sample_data(psdata_absolute)[[color_factor]])
      if (!is.null(shape_factor))
        sample_data(psdata_absolute)[[shape_factor]] <- as.factor(sample_data(psdata_absolute)[[shape_factor]])
      if (!is.null(size_factor))
        sample_data(psdata_absolute)[[size_factor]] <- as.factor(sample_data(psdata_absolute)[[size_factor]])
      if (!is.null(alpha_factor))
        sample_data(psdata_absolute)[[alpha_factor]] <- as.factor(sample_data(psdata_absolute)[[alpha_factor]])
    }

    # Define color, shape, and size sets based on psdata_relative
    if (!is.null(color_factor)) {
      sample_data(psdata_relative)[[color_factor]] <- as.factor(sample_data(psdata_relative)[[color_factor]])
      unique_colors <- levels(sample_data(psdata_relative)[[color_factor]])
      n_colors <- length(unique_colors)
      if (color_continuous == FALSE) {
        colorset <<- scales::hue_pal()(n_colors)
      }
    }
    if (!is.null(shape_factor)) {
      sample_data(psdata_relative)[[shape_factor]] <- as.factor(sample_data(psdata_relative)[[shape_factor]])
      n_shapes <- length(unique(sample_data(psdata_relative)[[shape_factor]]))
      shapeset <<- seq_len(n_shapes)
    }
    if (!is.null(size_factor)) {
      sample_data(psdata_relative)[[size_factor]] <- as.factor(sample_data(psdata_relative)[[size_factor]])
      n_sizes <- length(unique(sample_data(psdata_relative)[[size_factor]]))
      sizeset <<- seq(2, 2 + 1.2 * (n_sizes - 1), by = 1.2)
    }
    if (!is.null(alpha_factor)) {
      sample_data(psdata_relative)[[alpha_factor]] <- as.factor(sample_data(psdata_relative)[[alpha_factor]])
      alphaset <<- scale_alpha_continuous(range = c(0.3, 1))
    }

    na_types <- unique(sample_data(psdata_relative)$na_type)

    if (length(na_types) == 1) {
      # Create beta-diversity plots for a single na_type
      plot_Jac <- base_beta_plot(psdata_relative, ordination_method, "jaccard", "Jaccard\n(binary presence only)",
                                 color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")
      plot_BC <- base_beta_plot(psdata_relative, ordination_method, "bray", "Bray-Curtis\n(presence + abundance)",
                                color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")
      plot_uu <- base_beta_plot(psdata_relative, ordination_method, "uunifrac", "Unweighted UniFrac\n(lineage presence only)",
                                color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")
      plot_wu <- base_beta_plot(psdata_relative, ordination_method, "wunifrac", "Weighted UniFrac\n(lineage presence and abundance)",
                                color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")

      legend <- get_legend(plot_Jac + theme(legend.position = "right"))
      combined_plot_relative <- cowplot::plot_grid(plot_Jac, plot_BC, plot_uu, plot_wu,
                                                   ncol = 2, labels = c("A", "B", "C", "D"))
      combined_plot_relative <- cowplot::plot_grid(combined_plot_relative, legend, ncol = 2,
                                                   rel_widths = c(3, 0.8))

      print(combined_plot_relative)

      # Save the relative beta diversity plot (using the provided 'level' in the filename)
      figure_file_path <- paste0(asv_folder, project_name, "_beta_diversity_relative_", ordination_method, "_asv_level.pdf")
      ggsave(filename = figure_file_path, plot = combined_plot_relative, width = 16, height = 7)
      log_message(paste("Relative beta diversity plot saved:", figure_file_path), log_file)

      # Generate and save absolute beta diversity plots if norm_method is provided
      if (!is.null(norm_method)) {
        plot_man <- base_beta_plot(psdata_absolute, ordination_method, "manhattan", "Manhattan\n(PCoA)",
                                   color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "right")
        print(plot_man)
        figure_file_path <- paste0(asv_folder, project_name, "_beta_diversity_absolute_", ordination_method, "_asv_level.pdf")
        ggsave(filename = figure_file_path, plot = plot_man, width = 14, height = 7)
        log_message(paste("Absolute beta diversity plot saved:", figure_file_path), log_file)
      }
      return(combined_plot_relative)
    } else if (length(na_types) == 2) {
      # Process DNA and RNA separately
      psdata_relative_dna <- subset_samples(psdata_relative, na_type == "dna")
      plot_Jac_dna <- base_beta_plot(psdata_relative_dna, ordination_method, "jaccard", "Jaccard\n(binary presence only)",
                                     color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")
      plot_BC_dna <- base_beta_plot(psdata_relative_dna, ordination_method, "bray", "Bray-Curtis\n(presence + abundance)",
                                    color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")
      plot_uu_dna <- base_beta_plot(psdata_relative_dna, ordination_method, "uunifrac", "Unweighted UniFrac\n(lineage presence only)",
                                    color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")
      plot_wu_dna <- base_beta_plot(psdata_relative_dna, ordination_method, "wunifrac", "Weighted UniFrac\n(lineage presence and abundance)",
                                    color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")

      legend_dna <- get_legend(plot_Jac_dna + theme(legend.position = "right"))
      combined_plot_relative_dna <- cowplot::plot_grid(plot_Jac_dna, plot_BC_dna, plot_uu_dna, plot_wu_dna,
                                                       ncol = 2, labels = c("A", "B", "C", "D"))
      combined_plot_relative_dna <- cowplot::plot_grid(combined_plot_relative_dna, legend_dna, ncol = 2,
                                                       rel_widths = c(3, 0.8))

      print(combined_plot_relative_dna)

      figure_file_path <- paste0(asv_folder, project_name, "_beta_diversity_relative_", ordination_method, "_asv_level_dna.pdf")
      ggsave(filename = figure_file_path, plot = combined_plot_relative_dna, width = 16, height = 7)
      log_message(paste("Relative beta diversity DNA plot saved:", figure_file_path), log_file)

      if (!is.null(norm_method)) {
        psdata_absolute_dna <- subset_samples(psdata_absolute, na_type == "dna")
        plot_man_dna <- base_beta_plot(psdata_absolute_dna, ordination_method, "manhattan", "Manhattan\n(PCoA)",
                                       color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "right")
        print(plot_man_dna)

        figure_file_path <- paste0(asv_folder, project_name, "_beta_diversity_absolute_", ordination_method, "_asv_level_dna.pdf")
        ggsave(filename = figure_file_path, plot = plot_man_dna, width = 14, height = 7)
        log_message(paste("Absolute beta diversity DNA plot saved:", figure_file_path), log_file)
      }

      psdata_relative_rna <- subset_samples(psdata_relative, na_type == "rna")
      plot_Jac_rna <- base_beta_plot(psdata_relative_rna, ordination_method, "jaccard", "Jaccard\n(binary presence only)",
                                     color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")
      plot_BC_rna <- base_beta_plot(psdata_relative_rna, ordination_method, "bray", "Bray-Curtis\n(presence + abundance)",
                                    color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")
      plot_uu_rna <- base_beta_plot(psdata_relative_rna, ordination_method, "uunifrac", "Unweighted UniFrac\n(lineage presence only)",
                                    color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")
      plot_wu_rna <- base_beta_plot(psdata_relative_rna, ordination_method, "wunifrac", "Weighted UniFrac\n(lineage presence and abundance)",
                                    color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")

      legend_rna <- get_legend(plot_Jac_rna + theme(legend.position = "right"))
      combined_plot_relative_rna <- cowplot::plot_grid(plot_Jac_rna, plot_BC_rna, plot_uu_rna, plot_wu_rna,
                                                       ncol = 2, labels = c("A", "B", "C", "D"))
      combined_plot_relative_rna <- cowplot::plot_grid(combined_plot_relative_rna, legend_rna, ncol = 2,
                                                       rel_widths = c(3, 0.8))

      print(combined_plot_relative_rna)

      figure_file_path <- paste0(asv_folder, project_name, "_beta_diversity_relative_", ordination_method, "_asv_level_rna.pdf")
      ggsave(filename = figure_file_path, plot = combined_plot_relative_rna, width = 16, height = 7)
      log_message(paste("Relative beta diversity RNA plot saved:", figure_file_path), log_file)

      if (!is.null(norm_method)) {
        psdata_absolute_rna <- subset_samples(psdata_absolute, na_type == "rna")
        plot_man_rna <- base_beta_plot(psdata_absolute_rna, ordination_method, "manhattan", "Manhattan\n(PCoA)",
                                       color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "right")
        print(plot_man_rna)

        figure_file_path <- paste0(asv_folder, project_name, "_beta_diversity_absolute_", ordination_method, "_asv_level_rna.pdf")
        ggsave(filename = figure_file_path, plot = plot_man_rna, width = 14, height = 7)
        log_message(paste("Absolute beta diversity RNA plot saved:", figure_file_path), log_file)
      }
    }
  } else {
    for (tax in taxrank) {
      log_message(paste("Processing taxonomic level:", tax), log_file)

      if (is.null(norm_method)) {
        psdata_relative <- physeq
        psdata_absolute <- NULL
      } else if (norm_method == "fcm") {
        psdata_relative <- physeq[[paste0("psdata_copy_number_corrected_", tax)]]
        psdata_absolute <- physeq[[paste0("psdata_fcm_norm_rarefied_", tax)]]
      } else if (norm_method == "qpcr") {
        psdata_relative <- physeq[[paste0("psdata_copy_number_corrected_", tax)]]
        psdata_absolute <- physeq[[paste0("psdata_qpcr_norm_rarefied_", tax)]]
      }

      tax_folder <- paste0(beta_div_folder, tax, "/")
      if (!dir.exists(tax_folder)) {
        dir.create(tax_folder, recursive = TRUE)
      }

      psdata_relative <- transform_sample_counts(psdata_relative, function(x) x / sum(x) * 100)

      if (!is.null(norm_method)) {
        if (!is.null(color_factor))
          sample_data(psdata_absolute)[[color_factor]] <- as.factor(sample_data(psdata_absolute)[[color_factor]])
        if (!is.null(shape_factor))
          sample_data(psdata_absolute)[[shape_factor]] <- as.factor(sample_data(psdata_absolute)[[shape_factor]])
        if (!is.null(size_factor))
          sample_data(psdata_absolute)[[size_factor]] <- as.factor(sample_data(psdata_absolute)[[size_factor]])
        if (!is.null(alpha_factor))
          sample_data(psdata_absolute)[[alpha_factor]] <- as.factor(sample_data(psdata_absolute)[[alpha_factor]])
      }

      if (!is.null(color_factor)) {
        sample_data(psdata_relative)[[color_factor]] <- as.factor(sample_data(psdata_relative)[[color_factor]])
        unique_colors <- levels(sample_data(psdata_relative)[[color_factor]])
        n_colors <- length(unique_colors)
        if (color_continuous == FALSE) {
          colorset <<- scales::hue_pal()(n_colors)
        }
      }
      if (!is.null(shape_factor)) {
        sample_data(psdata_relative)[[shape_factor]] <- as.factor(sample_data(psdata_relative)[[shape_factor]])
        n_shapes <- length(unique(sample_data(psdata_relative)[[shape_factor]]))
        shapeset <<- seq_len(n_shapes)
      }
      if (!is.null(size_factor)) {
        sample_data(psdata_relative)[[size_factor]] <- as.factor(sample_data(psdata_relative)[[size_factor]])
        n_sizes <- length(unique(sample_data(psdata_relative)[[size_factor]]))
        sizeset <<- seq(2, 2 + 1.2 * (n_sizes - 1), by = 1.2)
      }
      if (!is.null(alpha_factor)) {
        sample_data(psdata_relative)[[alpha_factor]] <- as.factor(sample_data(psdata_relative)[[alpha_factor]])
        alphaset <<- scale_alpha_continuous(range = c(0.3, 1))
      }

      na_types <- unique(sample_data(psdata_relative)$na_type)

      if (length(na_types) == 1) {
        plot_Jac <- base_beta_plot(psdata_relative, ordination_method, "jaccard", "Jaccard\n(binary presence only)",
                                   color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")
        plot_BC <- base_beta_plot(psdata_relative, ordination_method, "bray", "Bray-Curtis\n(presence + abundance)",
                                  color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")
        plot_uu <- base_beta_plot(psdata_relative, ordination_method, "uunifrac", "Unweighted UniFrac\n(lineage presence only)",
                                  color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")
        plot_wu <- base_beta_plot(psdata_relative, ordination_method, "wunifrac", "Weighted UniFrac\n(lineage presence and abundance)",
                                  color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")

        legend <- get_legend(plot_Jac + theme(legend.position = "right"))
        combined_plot_relative <- cowplot::plot_grid(plot_Jac, plot_BC, plot_uu, plot_wu,
                                                     ncol = 2, labels = c("A", "B", "C", "D"))
        combined_plot_relative <- cowplot::plot_grid(combined_plot_relative, legend, ncol = 2,
                                                     rel_widths = c(3, 0.8))

        print(combined_plot_relative)

        figure_file_path <- paste0(tax_folder, project_name, "_beta_diversity_relative_", ordination_method, "_", tax, "_level.pdf")
        ggsave(filename = figure_file_path, plot = combined_plot_relative, width = 16, height = 7)
        log_message(paste("Relative beta diversity plot saved:", figure_file_path), log_file)

        if (!is.null(norm_method)) {
          plot_man <- base_beta_plot(psdata_absolute, ordination_method, "manhattan", "Manhattan\n(PCoA)",
                                     color_factor, shape_factor, size_factor, alpha_factor) +
            theme(legend.position = "right")
          print(plot_man)
          figure_file_path <- paste0(tax_folder, project_name, "_beta_diversity_absolute_", ordination_method, "_", tax, "_level.pdf")
          ggsave(filename = figure_file_path, plot = plot_man, width = 14, height = 7)
          log_message(paste("Absolute beta diversity plot saved:", figure_file_path), log_file)
        }
      } else if (length(na_types) == 2) {
        psdata_relative_dna <- subset_samples(psdata_relative, na_type == "dna")
        plot_Jac_dna <- base_beta_plot(psdata_relative_dna, ordination_method, "jaccard", "Jaccard\n(binary presence only)",
                                       color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")
        plot_BC_dna <- base_beta_plot(psdata_relative_dna, ordination_method, "bray", "Bray-Curtis\n(presence + abundance)",
                                      color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")
        plot_uu_dna <- base_beta_plot(psdata_relative_dna, ordination_method, "uunifrac", "Unweighted UniFrac\n(lineage presence only)",
                                      color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")
        plot_wu_dna <- base_beta_plot(psdata_relative_dna, ordination_method, "wunifrac", "Weighted UniFrac\n(lineage presence and abundance)",
                                      color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")

        legend_dna <- get_legend(plot_Jac_dna + theme(legend.position = "right"))
        combined_plot_relative_dna <- cowplot::plot_grid(plot_Jac_dna, plot_BC_dna, plot_uu_dna, plot_wu_dna,
                                                         ncol = 2, labels = c("A", "B", "C", "D"))
        combined_plot_relative_dna <- cowplot::plot_grid(combined_plot_relative_dna, legend_dna, ncol = 2,
                                                         rel_widths = c(3, 0.8))

        print(combined_plot_relative_dna)

        figure_file_path <- paste0(tax_folder, project_name, "_beta_diversity_relative_", ordination_method, "_", tax, "_level_dna.pdf")
        ggsave(filename = figure_file_path, plot = combined_plot_relative_dna, width = 16, height = 7)
        log_message(paste("Relative beta diversity DNA plot saved:", figure_file_path), log_file)

        if (!is.null(norm_method)) {
          psdata_absolute_dna <- subset_samples(psdata_absolute, na_type == "dna")
          plot_man_dna <- base_beta_plot(psdata_absolute_dna, ordination_method, "manhattan", "Manhattan\n(PCoA)",
                                         color_factor, shape_factor, size_factor, alpha_factor) +
            theme(legend.position = "right")
          print(plot_man_dna)
          plots_out[[tax]] <- plot_man_dna
          figure_file_path <- paste0(tax_folder, project_name, "_beta_diversity_absolute_", ordination_method, "_", tax, "_level_dna.pdf")
          ggsave(filename = figure_file_path, plot = plot_man_dna, width = 14, height = 7)
          log_message(paste("Absolute beta diversity DNA plot saved:", figure_file_path), log_file)
        }

        psdata_relative_rna <- subset_samples(psdata_relative, na_type == "rna")
        plot_Jac_rna <- base_beta_plot(psdata_relative_rna, ordination_method, "jaccard", "Jaccard\n(binary presence only)",
                                       color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")
        plot_BC_rna <- base_beta_plot(psdata_relative_rna, ordination_method, "bray", "Bray-Curtis\n(presence + abundance)",
                                      color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")
        plot_uu_rna <- base_beta_plot(psdata_relative_rna, ordination_method, "uunifrac", "Unweighted UniFrac\n(lineage presence only)",
                                      color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")
        plot_wu_rna <- base_beta_plot(psdata_relative_rna, ordination_method, "wunifrac", "Weighted UniFrac\n(lineage presence and abundance)",
                                      color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")

        legend_rna <- get_legend(plot_Jac_rna + theme(legend.position = "right"))
        combined_plot_relative_rna <- cowplot::plot_grid(plot_Jac_rna, plot_BC_rna, plot_uu_rna, plot_wu_rna,
                                                         ncol = 2, labels = c("A", "B", "C", "D"))
        combined_plot_relative_rna <- cowplot::plot_grid(combined_plot_relative_rna, legend_rna, ncol = 2,
                                                         rel_widths = c(3, 0.8))

        print(combined_plot_relative_rna)

        figure_file_path <- paste0(tax_folder, project_name, "_beta_diversity_relative_", ordination_method, "_", tax, "_level_rna.pdf")
        ggsave(filename = figure_file_path, plot = combined_plot_relative_rna, width = 16, height = 7)
        log_message(paste("Relative beta diversity RNA plot saved:", figure_file_path), log_file)

        if (!is.null(norm_method)) {
          psdata_absolute_rna <- subset_samples(psdata_absolute, na_type == "rna")
          plot_man_rna <- base_beta_plot(psdata_absolute_rna, ordination_method, "manhattan", "Manhattan\n(PCoA)",
                                         color_factor, shape_factor, size_factor, alpha_factor) +
            theme(legend.position = "right")

          print(plot_man_rna)

          figure_file_path <- paste0(tax_folder, project_name, "_beta_diversity_absolute_", ordination_method, "_", tax, "_level_rna.pdf")
          ggsave(filename = figure_file_path, plot = plot_man_rna, width = 14, height = 7)
          log_message(paste("Absolute beta diversity RNA plot saved:", figure_file_path), log_file)
        }
      }
    } # end for loop
  }
}
