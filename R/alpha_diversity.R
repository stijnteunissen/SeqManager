#' Generate Alpha Diversity Plots
#'
#' This function calculates alpha diversity metrics from a phyloseq object and
#' generates alpha diversity plots at various taxonomic levels or at the ASV
#' level. The alpha diversity measures (Observed, Chao1, Shannon, Simpson) are
#' computed using the \code{estimate_richness} function from the phyloseq
#' package. Depending on the input parameters, the function can handle
#' normalized data (using flow cytometry or qPCR methods) and can separate plots
#' by DNA and RNA types.
#'
#' @param physeq A phyloseq object containing microbial community data.
#' @param norm_method A character string specifying the normalization method for the data. Options include:
#'   \itemize{
#'     \item \code{"fcm"}: Use flow cytometry-normalized data.
#'     \item \code{"qpcr"}: Use qPCR-normalized data.
#'     \item \code{NULL}: Use only copy number corrected data (default).
#'   }
#' @param taxrank A character vector specifying the taxonomic levels for which
#'   alpha diversity is to be calculated. If the first element is (taxrank = `asv`),
#'   ASV-level data is processed. Otherwise, the function
#'   processes data for each taxonomic level provided default is taxrank = c('Phylum', 'Class', 'Order', 'Family', 'Tax_label').
#' @param date_factor An optional character string indicating the name of the
#'   date column in the sample metadata. If provided, the column is converted to
#'   a Date object ("%d/%m/%Y") and used to order the data.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the appropriate data object from \code{physeq} based on the chosen \code{norm_method} and taxonomic level.
#'   \item Estimates alpha diversity metrics (Observed, Chao1, Shannon, Simpson) using \code{estimate_richness}.
#'   \item Merges the alpha diversity estimates with sample metadata.
#'   \item Optionally orders the data by a date factor if provided.
#'   \item Appends dummy data to the dataset for visualization purposes.
#'   \item Exports the combined alpha diversity data as CSV files.
#'   \item Generates bar plots for the diversity metrics, optionally separating DNA and RNA data if both are present.
#'   \item Saves the resulting plots as PDF files in the project's figures folder.
#' }
#'
#' @return A combined ggplot object containing the generated alpha diversity plots.
#'
#' @examples
#' \dontrun{
#'   # Generate alpha diversity plots at the ASV level without normalization
#'   alpha_plot <- alpha_diversity(physeq = my_physeq, taxrank = "asv")
#'
#'   # Generate alpha diversity plots at the Phylum level using flow cytometry-normalized data
#'   alpha_plot <- alpha_diversity(physeq = my_physeq, norm_method = "fcm", taxrank = "Phylum")
#'
#'   # Generate alpha diversity plots with a specified date factor for ordering samples
#'   alpha_plot <- alpha_diversity(physeq = my_physeq, date_factor = "Sample_Date")
#' }
#'
#' @export
alpha_diversity = function(physeq = physeq,
                           norm_method = NULL,
                           taxrank = c("Phylum", "Class", "Order", "Family", "Tax_label"),
                           date_factor = NULL) {

  base_alpha_plot = function(alpha_data, x_value, y_value, x_label, y_label) {

    if (!is.null(date_factor)) {
      alpha_data = alpha_data %>%
        mutate(!!date_factor := as.Date(.data[[date_factor]], format = "%d/%m/%Y")) %>%
        arrange(.data[[date_factor]])

      alpha_data = alpha_data %>%
        mutate(grouping_factor = do.call(paste, c(across(all_of(present_factors)), sep = "_")))

      plot = ggplot(alpha_data, aes(x = !!sym(x_value), y = !!sym(y_value), group = grouping_factor)) +
        #geom_jitter(aes(color = .data[[date_factor]]), size = 2, width = 0.2, show.legend = FALSE) +
        #scale_color_date(low = "lightblue", high = "darkgreen") +
        geom_col(fill = "steelblue", color = "steelblue", show.legend = FALSE) +
        theme_classic() +
        labs(x = x_label, y = y_label) +
        theme(
          legend.position = "bottom",
          legend.text = element_markdown(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(face = "bold"),
          #strip.text = element_text(face = "bold", angle = 90, vjust = 0.5, hjust = 0),
          strip.background = element_blank(),
          ggh4x.facet.nestline = element_line(colour = "black")
        ) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +  # Correct placement outside theme with 5%
        expand_limits(y = c(min(alpha_data[[y_value]]) - 1, max(alpha_data[[y_value]]) + 1))  # Correct placement outside theme

    } else {
      alpha_data = alpha_data %>%
        mutate(grouping_factor = do.call(paste, c(across(all_of(present_factors)), sep = "_")))

      plot = ggplot(alpha_data, aes(x = !!sym(x_value), y = !!sym(y_value), group = grouping_factor)) +
        #geom_jitter(aes(color = grouping_factor), size = 2, width = 0.2, show.legend = FALSE) +
        geom_col(fill = "steelblue", color = "steelblue", show.legend = FALSE) +
        #scale_color_manual(values = colorset) +
        theme_classic() +
        labs(x = x_label, y = y_label) +
        theme(
          legend.position = "bottom",
          legend.text = element_markdown(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(face = "bold"),
          #strip.text = element_text(angle = 90, vjust = 0.5, hjust = 0),
          strip.background = element_blank(),
          ggh4x.facet.nestline = element_line(colour = "black")) +
        scale_y_continuous(expand = c(0, 0))
    }
    return(plot)
  }

  facet_add = function(present_factors) {
    if (!is.null(present_factors) && length(present_factors) > 0) {
      return(facet_nested(cols = vars(!!!syms(present_factors)), scales = "free", space = "free", nest_line = element_line(linetype = 1)))
    } else {
      return(NULL)
    }
  }

  project_name = projects
  project_folder = paste0(base_path, project_name)
  figure_folder = paste0(project_folder, "/figures/")
  output_folder_csv_files = paste0(project_folder, "/output_data/csv_files/")
  output_folder_rds_files = paste0(project_folder, "/output_data/rds_files/")

  if (tolower(taxrank[1]) == "asv") {
    log_message("Processing ASV-level alpha diversity", log_file)

    if (is.null(norm_method)) {
      psdata = physeq[["psdata_asv_copy_number_corrected"]]
    } else if (norm_method == "fcm") {
      psdata = physeq[["psdata_asv_fcm_norm_rarefied"]]
    } else if (norm_method == "qpcr") {
      psdata = physeq[["psdata_asv_qpcr_norm_rarefied"]]
    }

    alpha_div_folder = paste0(figure_folder, "Alpha_diversity/")
    if(!dir.exists(alpha_div_folder)){dir.create(alpha_div_folder)}

    asv_folder = paste0(alpha_div_folder, "ASV/")
    if(!dir.exists(asv_folder)){dir.create(asv_folder)}

    variable_columns = intersect(present_variable_factors, colnames(sample_data(psdata)))
    factor_columns = unique(c(variable_columns))
    present_factors = if (length(factor_columns) > 0) factor_columns else NULL

    alpha_data = estimate_richness(psdata, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
    alpha_data = alpha_data %>% rownames_to_column(var = "sampleid")
    metadata = sample_data(psdata) %>% data.frame() %>% as_tibble()
    alpha_data_full = inner_join(metadata, alpha_data, by = "sampleid")

    if ("Sample_Date" %in% present_factors) {
      alpha_data_full = alpha_data_full %>%
        mutate(Sample_Date = as.Date(Sample_Date, format = "%d/%m/%Y")) %>%
        arrange(Sample_Date)
    }

    # # adding dummy data to the dataset.
    # dummy_sample_names <- sprintf("DUMMY_%03d", 1:12)
    # dummy_row <- tibble(
    #   sampleid = dummy_sample_names,
    #   timepoint = rep(1, 12),
    #   Chao1 = 0,
    #   Shannon = 0,
    #   Observed = 0,
    #   Simpson = 0,
    #   na_type = "dna",
    #   soil_type = rep(c("L1", "L5", "L7"), each = 4),
    #   treatment = rep(c("Non-Pesticide", "Pesticide", "Pesticide", "Pesticide"), times = 3),
    #   replica   = rep(c(4, 1, 2, 3), times = 3)
    # )
    # alpha_data_full <- bind_rows(alpha_data_full, dummy_row)

    alpha_div_csv_folder = paste0(output_folder_csv_files, "Alpha_diversity/")
    if(!dir.exists(alpha_div_csv_folder)){dir.create(alpha_div_csv_folder)}

    asv_csv_folder = paste0(alpha_div_csv_folder, "ASV/")
    if(!dir.exists(asv_csv_folder)){dir.create(asv_csv_folder)}

    alpha_data_full_csv = alpha_data_full %>% mutate(Observed = round(Observed, 2),
                                                     Chao1 = round(Chao1, 2),
                                                     Shannon = round(Shannon, 2),
                                                     Simpson = round(Simpson, 2))

    output_file_path = paste0(asv_csv_folder, project_name, "_alpha_diversity_asv_level.csv")
    write.csv(alpha_data_full_csv, file = output_file_path, row.names = FALSE)
    log_message(paste("Alpha diversity asv level saved as .csv object in", output_file_path), log_file)

    na_types = unique(alpha_data_full$na_type)

    if (length(na_types) == 1) {
      chao1_plot = base_alpha_plot(alpha_data_full, "sampleid", "Chao1", x_label = "Sample", y_label = "Chao1 Index") +
        facet_add(present_factors)

      shannon_plot = base_alpha_plot(alpha_data_full, "sampleid", "Shannon", x_label = "Sample", y_label = "Shannon Index") +
        facet_add(present_factors)

      combined_plot = plot_grid(chao1_plot + theme(legend.position = "none"), shannon_plot + theme(legend.position = "none"),
                                align = "hv", labels = c("A", "B"), nrow = 1)

    } else if (length(na_types) == 2) {
      alpha_data_full_dna = alpha_data_full %>% filter(na_type == "dna")

      chao1_plot_dna =
        base_alpha_plot(alpha_data_full_dna, "sampleid", "Chao1", x_label = "Sample", y_label = "Chao1 Index") +
        facet_add(present_factors)

      shannon_plot_dna =
        base_alpha_plot(alpha_data_full_dna, "sampleid", "Shannon", x_label = "Sample", y_label = "Shannon Index") +
        facet_add(present_factors)

      alpha_data_full_rna = alpha_data_full %>% filter(na_type == "rna")

      chao1_plot_rna =
        base_alpha_plot(alpha_data_full_rna, "sampleid", "Chao1", x_label = "Sample", y_label = "Chao1 Index") +
        facet_add(present_factors)

      shannon_plot_rna =
        base_alpha_plot(alpha_data_full_rna, "sampleid", "Shannon", x_label = "Sample", y_label = "Shannon Index") +
        facet_add(present_factors)

      separator_line = ggdraw() +
        draw_line(x = c(0.25, 0.75), y = c(0.5, 0.5), size = 0.5, color = "black") +
        theme_void()

      dna_label = ggdraw() + draw_label("DNA", fontface = "bold", size = 14, hjust = 0.5)
      rna_label = ggdraw() + draw_label("RNA", fontface = "bold", size = 14, hjust = 0.5)

      combined_plot_dna = plot_grid(chao1_plot_dna, shannon_plot_dna,
                                    ncol = 2, labels = c("A", "B"))

      combined_plot_rna = plot_grid(chao1_plot_rna, shannon_plot_rna,
                                    ncol = 2, labels = c("C", "D"))

      combined_plot = plot_grid(
        plot_grid(dna_label, separator_line, combined_plot_dna, ncol = 1, rel_heights = c(0.1, 0.05, 1)),
        plot_grid(rna_label, separator_line, combined_plot_rna, ncol = 1, rel_heights = c(0.1, 0.05, 1)),
        ncol = 1)
    }

    print(combined_plot)

    figure_file_path = paste0(asv_folder, project_name, "_alpha_diversity_asv_level.pdf")
    ggsave(filename = figure_file_path, plot = combined_plot, width = 30, height = 8)
    log_message(paste("alpha diversity asv level saved as .pdf object in", figure_file_path), log_file)
  } else {

    for (tax in taxrank) {
      log_message(paste("Processing taxonomic level:", tax), log_file)

      if (is.null(norm_method)) {
        psdata = physeq
      } else if (norm_method == "fcm") {
        psdata = physeq[[paste0("psdata_fcm_norm_rarefied_", tax)]]
      } else if (norm_method == "qpcr") {
        psdata = physeq[[paste0("psdata_qpcr_norm_rarefied_", tax)]]
      }

      alpha_div_folder = paste0(figure_folder, "Alpha_diversity/")
      if(!dir.exists(alpha_div_folder)){dir.create(alpha_div_folder)}

      tax_folder = paste0(alpha_div_folder, tax, "/")
      if(!dir.exists(tax_folder)){dir.create(tax_folder)}

      variable_columns = intersect(present_variable_factors, colnames(sample_data(psdata)))
      factor_columns = unique(c(variable_columns))
      present_factors = if (length(factor_columns) > 0) factor_columns else NULL

      alpha_data = estimate_richness(psdata, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
      alpha_data = alpha_data %>% rownames_to_column(var = "sampleid")
      metadata = sample_data(psdata) %>% data.frame() %>% as_tibble()
      alpha_data_full = inner_join(metadata, alpha_data, by = "sampleid")

      if ("Sample_Date" %in% present_factors) {
        alpha_data_full = alpha_data_full %>%
          mutate(Sample_Date = as.Date(Sample_Date, format = "%d/%m/%Y")) %>%
          arrange(Sample_Date)
      }

      # # adding dummy data to the dataset.
      # dummy_sample_names <- sprintf("DUMMY_%03d", 1:12)
      # dummy_row <- tibble(
      #   sampleid = dummy_sample_names,
      #   timepoint = rep(1, 12),
      #   Chao1 = 0,
      #   Shannon = 0,
      #   Observed = 0,
      #   Simpson = 0,
      #   na_type = "dna",
      #   soil_type = rep(c("L1", "L5", "L7"), each = 4),
      #   treatment = rep(c("Non-Pesticide", "Pesticide", "Pesticide", "Pesticide"), times = 3),
      #   replica   = rep(c(4, 1, 2, 3), times = 3)
      # )
      # alpha_data_full <- bind_rows(alpha_data_full, dummy_row)

      alpha_div_csv_folder = paste0(output_folder_csv_files, "Alpha_diversity/")
      if(!dir.exists(alpha_div_csv_folder)){dir.create(alpha_div_csv_folder)}

      tax_csv_folder = paste0(alpha_div_csv_folder, tax, "/")
      if(!dir.exists(tax_csv_folder)){dir.create(tax_csv_folder)}

      alpha_data_full_csv = alpha_data_full %>% mutate(Observed = round(Observed, 2),
                                                       Chao1 = round(Chao1, 2),
                                                       Shannon = round(Shannon, 2),
                                                       Simpson = round(Simpson, 2))

      output_file_path = paste0(tax_csv_folder, project_name, "_alpha_diversity_", tax, "_level.csv")
      write.csv(alpha_data_full_csv, file = output_file_path, row.names = FALSE)
      log_message(paste("Alpha diversity", tax, "level saved as .csv object in", output_file_path), log_file)

      na_types = unique(alpha_data_full$na_type)

      if (length(na_types) == 1) {
        chao1_plot = base_alpha_plot(alpha_data_full, "sampleid", "Chao1", x_label = "Sample", y_label = "Chao1 Index") +
          facet_add(present_factors)

        shannon_plot = base_alpha_plot(alpha_data_full, "sampleid", "Shannon", x_label = "Sample", y_label = "Shannon Index") +
          facet_add(present_factors)

        combined_plot = plot_grid(chao1_plot + theme(legend.position = "none"), shannon_plot + theme(legend.position = "none"),
                                  align = "hv", labels = c("A", "B"), nrow = 1)

      } else if (length(na_types) == 2) {
        alpha_data_full_dna = alpha_data_full %>% filter(na_type == "dna")

        chao1_plot_dna =
          base_alpha_plot(alpha_data_full_dna, "sampleid", "Chao1", x_label = "Sample", y_label = "Chao1 Index") +
          facet_add(present_factors)

        shannon_plot_dna =
          base_alpha_plot(alpha_data_full_dna, "sampleid", "Shannon", x_label = "Sample", y_label = "Shannon Index") +
          facet_add(present_factors)

        alpha_data_full_rna = alpha_data_full %>% filter(na_type == "rna")

        chao1_plot_rna =
          base_alpha_plot(alpha_data_full_rna, "sampleid", "Chao1", x_label = "Sample", y_label = "Chao1 Index") +
          facet_add(present_factors)

        shannon_plot_rna =
          base_alpha_plot(alpha_data_full_rna, "sampleid", "Shannon", x_label = "Sample", y_label = "Shannon Index") +
          facet_add(present_factors)

        separator_line = ggdraw() +
          draw_line(x = c(0.25, 0.75), y = c(0.5, 0.5), size = 0.5, color = "black") +
          theme_void()

        dna_label = ggdraw() + draw_label("DNA", fontface = "bold", size = 14, hjust = 0.5)
        rna_label = ggdraw() + draw_label("RNA", fontface = "bold", size = 14, hjust = 0.5)

        combined_plot_dna = plot_grid(chao1_plot_dna, shannon_plot_dna,
                                      ncol = 2, labels = c("A", "B"))

        combined_plot_rna = plot_grid(chao1_plot_rna, shannon_plot_rna,
                                      ncol = 2, labels = c("C", "D"))

        combined_plot = plot_grid(
          plot_grid(dna_label, separator_line, combined_plot_dna, ncol = 1, rel_heights = c(0.1, 0.05, 1)),
          plot_grid(rna_label, separator_line, combined_plot_rna, ncol = 1, rel_heights = c(0.1, 0.05, 1)),
          ncol = 1)
      }

      print(combined_plot)

      figure_file_path = paste0(tax_folder, project_name, "_alpha_diversity_", tax, "_level.pdf")
      ggsave(filename = figure_file_path, plot = combined_plot, width = 30, height = 8)
      log_message(paste("alpha diversity", tax, "level saved as .pdf object in", figure_file_path), log_file)
    }
  }
  return(combined_plot)
}
