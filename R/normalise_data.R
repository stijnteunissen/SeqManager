#' Normalize Phyloseq Data
#'
#' This function applies normalization to a `phyloseq` object, allowing the
#' conversion of data to absolute values based on either flow cytometry (FCM)
#' data or qPCR data. The funciton accounts for DNA and RNA normalization,
#' considering extra steps involved in prodessing RNA. Additionally, it
#' incorporates Anna16 correction to adjust ASV counts based on predicted copy
#' numbers. The normalized data can later be used for relative abundance
#' calculations.
#'
#' @inheritParams remove_mock
#'
#' @param norm_method A character string specifying the normalization method to use. Options are:
#'   \itemize{
#'     \item `"fcm"`: Normalize based on flow cytometry data, converting abundances to cell concentrations (cells/mL or gram sample).
#'     \item `"qpcr"`: Normalize based on qPCR data, converting abundances to cell equivalents (cells/ml or gram sample).
#'     \item `NULL`: Apply only Anna16 correction without additional normalization.
#'   }
#'
#' @details
#' The function performs the following steps based on the chosen normalization method:
#'
#' 1. **Anna16 Correction (applied in all cases):**
#'    - Correct ASV abundances by dividing each count by the predicted copy number. The predicted 16S rRNA gene copy numbers are
#'    based on the method described in ["Deep learning for predicting 16S rRNA gene copy number"](https://dx.doi.org/10.1038/s41598-024-64658-5).
#'    - This prediction enables the calculation of cell equivalents by adjusting for variablility in 16S rRNA gene copy numbers across taxa.
#'
#' 2. **FCM Normalization (`norm_method = "fcm"`):**
#'    - Combine sample metadata (e.g., `cells_per_ml`, `sample_volume_ml`) with corrected abundances.
#'    - Calculate absolute abundances as `cells_per_ml_sample` and adjust ASV abundances accordingly.
#'
#' 3. **qPCR Normalization (`norm_method = "qpcr"`):**
#'    - Use qPCR-derived `sq_calc_mean` and predicted copy numbers to calculate absolute abundances.
#'    - Normalize abundances to represent cell equivalents, ensuring compatibility with downstream analysis.
#'
#' ### DNA vs RNA Normalization
#' The interpretation of normalized data depends on the type of nucleic acid (DNA or RNA):
#'   - **DNA:** Normalized abundances typically represent **cells per ml (or gram)**. This assumes one genome copy per cell.
#'   - **RNA:** Normalized abundances often represent **copies per cell equivalent per ml (or gram)**.
#'              RNA reflects transcriptional activity and can vary widely depending on the condition of the cells.
#'
#' When working with RNA data, ensure that metadata contains information about RNA extraction efficiency or other factors influencing RNA copy numbers.
#'
#' @references
#' Miao J, Chen T, Misir M, Lin Y. Deep learning for predicting 16S rRNA gene
#' copy number. Sci Rep. 2024 Jun 20;14(1):14282. doi:
#' [10.1038/s41598-024-64658-5](https://dx.doi.org/10.1038/s41598-024-64658-5).
#' PMID: 38902329; PMCID: PMC11190246.
#'
#' @return
#' A list containing one or more `phyloseq` objects:
#'   \itemize{
#'     \item `psdata_anna16_corrected`: Phyloseq object after Anna16 correction.
#'     \item `psdata_fcm_norm`: (if `norm_method = "fcm"`) Phyloseq object with FCM-normalized cell concentrations.
#'     \item `psdata_qpcr_norm`: (if `norm_method = "qpcr"`) Phyloseq object with qPCR-normalized cell equivalents.
#'   }
#' The normalized objects are also saved as `.rds` files in the specified output directory.
#'
#' @examples
#' \dontrun{
#' # Apply Anna16 correction only
#' result <- normalize_data(physeq = physeq, norm_method = NULL)
#'
#' # Normalize using flow cytometry (FCM) data
#' result <- normalize_data(physeq = physeq, norm_method = "fcm")
#'
#' # Normalize using qPCR data
#' result <- normalize_data(physeq = physeq, norm_method = "qpcr")
#' }
#'
#' @export
normalise_data = function(physeq = without_mock_physeq,
                          norm_method = NULL) {

  psdata = physeq
  project_name = projects

  project_folder = paste0(base_path, project_name)
  destination_folder = paste0(project_folder, "/input_data/")
  output_folder_rds_files_before = paste0(project_folder, "/output_data/rds_files/Before_cleaning_rds_files/")
  output_folder_rds_files_after = paste0(project_folder, "/output_data/rds_files/After_cleaning_rds_files/")
  output_asv_rds_files = paste0(output_folder_rds_files_after, "ASV/")
  if(!dir.exists(output_asv_rds_files)){dir.create(output_asv_rds_files)}

  # anna16 correctie
  df_psdata = data.frame(otu_table(psdata))
  df_psdata$OTU = rownames(df_psdata)
  pstibble =
    as_tibble(df_psdata) %>%
    select(OTU, everything())

  anna16_file = list.files(destination_folder, pattern = "dna-sequences\\.csv$", full.names = TRUE)
  anna16_data = read_csv(anna16_file, show_col_types = FALSE) %>% rename(index = "OTU")

  joined_pstibble =
    pstibble %>%
    inner_join(., anna16_data, by = "OTU")

  corrected_joined_pstibble =
    joined_pstibble %>%
    rowwise() %>%
    mutate(across(c(everything(), -OTU, -predicted_copy_number), ~ . / predicted_copy_number)) %>%
    mutate(across(c(everything(), -OTU, -predicted_copy_number), ~ ceiling(.))) %>%
    select(-predicted_copy_number)

  otu_corrected = otu_table(data.frame(corrected_joined_pstibble[, -1]), taxa_are_rows = TRUE)
  taxa_names(otu_corrected) = corrected_joined_pstibble$OTU

  psdata_anna16_corrected = psdata
  otu_table(psdata_anna16_corrected) = otu_corrected

  # output_file_path = paste0(output_folder_rds_files, project_name, "_phyloseq_asv_level_anna16_corrected_counts.rds")
  # saveRDS(psdata_anna16_corrected, file = output_file_path)
  # log_message(paste("Anna16 corrected phyloseq data saved as .rds object in", output_file_path), log_file)

  # if (is.null(norm_method)) {
  #   log_message("No normalization method specified for absolute data. Only anna16 correction applied.", log_file)
  #   return(list(psdata_anna16_corrected = psdata_anna16_corrected))
  # }

  # FCM normalisatie
  if (!is.null(norm_method) && norm_method == "fcm") {

    df_psdata_fcm = data.frame(otu_table(psdata_anna16_corrected))
    df_psdata_fcm$OTU = rownames(df_psdata_fcm)

    psdata_fcm_long =
      df_psdata_fcm %>%
      pivot_longer(cols = -OTU,
                   names_to = "SampleID",
                   values_to = "Abundance")

    df_sample_data_fcm = data.frame(sample_data(psdata_anna16_corrected))
    df_sample_data_fcm$SampleID = rownames(df_sample_data_fcm)

    joined_pstibble_fcm =
      psdata_fcm_long %>%
      inner_join(., df_sample_data_fcm, by = "SampleID")

    joined_pstibble_fcm_norm =
      joined_pstibble_fcm %>%
      rowwise() %>%
      #mutate(cells_per_ml = cells_per_ml / sample_volume_ml) %>%
      mutate(norm_abund = ceiling(Abundance * cells_per_ml)) %>%
      select(OTU, norm_abund, SampleID)

    fcm_norm_wide =
      joined_pstibble_fcm_norm %>%
      pivot_wider(names_from = SampleID,
                  values_from = norm_abund)

    otu_fcm_norm = otu_table(data.frame(fcm_norm_wide[, -1]), taxa_are_rows = TRUE)
    taxa_names(otu_fcm_norm) = fcm_norm_wide$OTU

    psdata_fcm_norm = psdata_anna16_corrected
    otu_table(psdata_fcm_norm) = otu_fcm_norm

    # return(list(psdata_anna16_corrected = psdata_anna16_corrected, psdata_fcm_norm = psdata_fcm_norm))

    # qpcr normalisatie
  } else if (!is.null(norm_method) && norm_method == "qpcr") {

    df_psdata_qpcr = data.frame(otu_table(psdata))
    df_psdata_qpcr$OTU = rownames(df_psdata_qpcr)

    psdata_qpcr_long =
      df_psdata_qpcr %>%
      pivot_longer(cols = -OTU,
                   names_to = "SampleID",
                   values_to = "Abundance")

    df_sample_data_qpcr = data.frame(sample_data(psdata))
    df_sample_data_qpcr$SampleID = rownames(df_sample_data_qpcr)

    joined_pstibble_qpcr =
      psdata_qpcr_long %>%
      inner_join(., df_sample_data_qpcr, by = "SampleID")

    # Merge `predicted_copy_number` from Anna16 data
    joined_pstibble_qpcr_anna16 =
      joined_pstibble_qpcr %>%
      inner_join(., anna16_data %>% select(OTU, predicted_copy_number), by = "OTU")

    # cell equivlants abundance in plaats van relative
    # abdance is cell equivlants
    # aangeven wat de eenheiden zijn vor relatief als het relatief als is of absoluut en de juiste eenheid kan per gram of ml

    results_list = list()

    if ("dna" %in% joined_pstibble_qpcr_anna16$na_type) {
      joined_pstibble_qpcr_norm_dna =
        joined_pstibble_qpcr_anna16 %>%
        filter(na_type == "dna") %>%
        group_by(SampleID) %>%
        mutate(a = (sq_mean / insert_volume) * dilution_factor,
               b = a * elution_volume_ul,
               sq_calc_mean = b / dry_sample_mass_g) %>%
        ungroup()
      results_list$dna = joined_pstibble_qpcr_norm_dna

      # Add sq_calc_mean to sample_data
      df = data.frame(sample_data(psdata))
      df = df %>% mutate(SampleID = rownames(df))
      df_2 = df %>% left_join(joined_pstibble_qpcr_norm_dna %>% select(SampleID, sq_calc_mean) %>% distinct(), by = "SampleID")

      original_sample_names <- sample_names(psdata)
      rownames(df_2) <- original_sample_names

      sample_data(psdata) = sample_data(df_2)

    }

    if ("rna" %in% joined_pstibble_qpcr_anna16$na_type) {
      joined_pstibble_qpcr_norm_rna =
        joined_pstibble_qpcr_anna16 %>%
        filter(na_type == "rna") %>%
        group_by(SampleID) %>%
        mutate(a = (sq_mean / insert_volume) * dilution_factor,
               b = a * Final_volume_of_cDNA,
               c = b * (RNA_volume_used_Dnase_treatment / Dnase_Volume_used_cDNA_synthesis),
               d = c * (elution_volume_ul / RNA_volume_used_Dnase_treatment),
               sq_calc_mean = d / dry_sample_mass_g) %>%
        ungroup()
      results_list$rna = joined_pstibble_qpcr_norm_rna

      # Add sq_calc_mean to sample_data
      df = data.frame(sample_data(psdata))
      df = df %>% mutate(SampleID = rownames(df))
      df_2 = df %>% left_join(joined_pstibble_qpcr_norm_rna %>% select(SampleID, sq_calc_mean) %>% distinct(), by = "SampleID")

      original_sample_names <- sample_names(psdata)
      rownames(df_2) <- original_sample_names

      sample_data(psdata) = sample_data(df_2)
    }

    if ("dna" %in% joined_pstibble_qpcr_anna16$na_type & "rna" %in% joined_pstibble_qpcr_anna16$na_type) {
      df = data.frame(sample_data(psdata))
      df_2 = df %>%
        mutate(sq_calc_mean = coalesce(sq_calc_mean.x, sq_calc_mean.y)) %>%
        select(-sq_calc_mean.x, -sq_calc_mean.y)
      sample_data(psdata) = sample_data(df_2)
    }

    if (!is.null(results_list$dna) && !is.null(results_list$rna)) {
      joined_pstibble_combined =
        bind_rows(results_list$dna, results_list$rna) %>%
        group_by(SampleID) %>%
        mutate(relative_abundance = Abundance / sum(Abundance),
               absolute_abundance_qpcr = relative_abundance * sq_calc_mean,
               norm_abund = ceiling(absolute_abundance_qpcr / predicted_copy_number)) %>%
        ungroup() %>%
        select(OTU, norm_abund, SampleID)
    } else if (!is.null(results_list$dna)) {
      joined_pstibble_combined =
        results_list$dna %>%
        group_by(SampleID) %>%
        mutate(relative_abundance = Abundance / sum(Abundance),
               absolute_abundance_qpcr = relative_abundance * sq_calc_mean,
               norm_abund = ceiling(absolute_abundance_qpcr / predicted_copy_number)) %>%
        ungroup() %>%
        select(OTU, norm_abund, SampleID)
    } else if (!is.null(results_list$rna)) {
      joined_pstibble_combined =
        results_list$rna %>%
        group_by(SampleID) %>%
        mutate(relative_abundance = Abundance / sum(Abundance),
               absolute_abundance_qpcr = relative_abundance * sq_calc_mean,
               norm_abund = ceiling(absolute_abundance_qpcr / predicted_copy_number)) %>%
        ungroup() %>%
        select(OTU, norm_abund, SampleID)
    }

    qpcr_norm_wide =
      joined_pstibble_combined %>%
      pivot_wider(names_from = SampleID,
                  values_from = norm_abund)

    otu_qpcr_norm = otu_table(data.frame(qpcr_norm_wide[, -1]), taxa_are_rows = TRUE)
    taxa_names(otu_qpcr_norm) = qpcr_norm_wide$OTU

    psdata_qpcr_norm = psdata
    otu_table(psdata_qpcr_norm) = otu_qpcr_norm

  }

  # modify tax table
  modify_tax_table = function(psdata) {
    tax_table_df = as.data.frame(tax_table(psdata))
    otu_names = taxa_names(psdata)

    tax_table_df =
      tax_table_df %>%
      rowwise() %>%
      mutate(Species = case_when(
        grepl("\\d", Genus) ~ {
          first_non_numeric <- case_when(
            !grepl("\\d", Family) ~ paste(Family, Genus, sep = " "),
            !grepl("\\d", Order) ~ paste(Order, Genus, sep = " "),
            !grepl("\\d", Class) ~ paste(Class, Genus, sep = " "),
            !grepl("\\d", Phylum) ~ paste(Phylum, Genus, sep = " "),
            !grepl("\\d", Kingdom) ~ paste(Kingdom, Genus, sep = " "),
            TRUE ~ NA_character_
          )
          first_non_numeric
        },
        TRUE ~ Genus)) %>%
      ungroup()

    tax_table_df =
      tax_table_df %>%
      dplyr::rename(Tax_label = Species)

    tax_table_matrix = as.matrix(tax_table_df)
    rownames(tax_table_matrix) = otu_names
    tax_table(psdata) = tax_table_matrix

    return(psdata)
  }

  psdata_anna16_corrected = modify_tax_table(psdata_anna16_corrected)

  output_file_path = paste0(output_asv_rds_files, project_name, "_phyloseq_asv_level_anna16_corrected_counts.rds")
  saveRDS(psdata_anna16_corrected, file = output_file_path)
  log_message(paste("phyloseq data anna16 corrected counts asv level saved as .rds object in", output_file_path), log_file)

  if (is.null(norm_method)) {
    log_message("No normalization method specified for absolute data. Only anna16 correction applied.", log_file)
    return(list(psdata_asv_anna16_corrected = psdata_anna16_corrected))
  }

  if (norm_method == "fcm") {
    psdata_fcm_norm = modify_tax_table(psdata_fcm_norm)

    output_file_path = paste0(output_folder_rds_files_before, project_name, "_phyloseq_asv_level_fcm_normalised_cell_concentration.rds")
    saveRDS(psdata_fcm_norm, file = output_file_path)
    log_message(paste("Phyloseq data fcm normalised cell concentration (cells per ml/gram sample) asv level saved as .rds object in", output_file_path), log_file)

    return(list(psdata_asv_anna16_corrected = psdata_anna16_corrected, psdata_asv_fcm_norm = psdata_fcm_norm))

  } else if (norm_method == "qpcr") {
    psdata_qpcr_norm = modify_tax_table(psdata_qpcr_norm)

    output_file_path = paste0(output_folder_rds_files_before, project_name, "_phyloseq_asv_level_qpcr_normalised_celL_concentration.rds")
    saveRDS(psdata_qpcr_norm, file = output_file_path)
    log_message(paste("phyloseq qpcr normalised cell concentration (cells per ml/gram sample) asv level saved as .rds object in", output_file_path), log_file)

    return(list(psdata_asv_anna16_corrected = psdata_anna16_corrected, psdata_asv_qpcr_norm = psdata_qpcr_norm))

  } else {
    log_message("Error: Invalid normalization method specified. Use 'fcm' or 'qpcr'.")
    stop("Error: Invalid normalization method specified. Use 'fcm' or 'qpcr'.")
  }

}
