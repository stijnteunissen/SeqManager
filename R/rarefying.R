#' Rarefying Phyloseq Data
#'
#' This function performs rarefaction on a normalized `phyloseq` object. Rarefaction is applied either based on
#' flow cytometry (FCM) data or qPCR data, depending on the specified normalization method. The rarefied data is
#' saved as `.rds` files for downstream analysis.
#'
#' @param physeq A list containing the `phyloseq` objects. Required components:
#'               \itemize{
#'               \item `psdata_asv_copy_number_corrected`: With copy number corrected phyloseq data.
#'               \item `psdata_asv_fcm_norm`: Required if `norm_method` = "fcm".
#'               \item `psdata_asv_qpcr_norm`: Required if `norm_method` = "qpcr".
#'               }
#'
#' @param norm_method A string specifying the normalization method. Acceptable values:
#'                    \itemize{
#'                    \item `NULL`: Returns the `psdata_asv_copy_number_corrected` object without modifications.
#'                    \item `"fcm"`: Rarefied data using FCM normalized.
#'                    \item `"qpcr"`: Rarefied data using qPCR normalized.
#'                    }
#'
#' @details
#' - For `"fcm"` normalization:
#'   - Rarefies only the `fcm` normalized data, while the `copy number corrected` data remains unchanged.
#'   - Rarefies based on the calculated sampling depth, derived from the ratio of total reads per sample to the
#'     estimated cell counts.
#'   - Uses a custom function (`avgrarefy`) to perform multiple iterations of rarefaction and averages the results.
#' - For `"qpcr"` normalization:
#'   - Rarefies only the `qpcr` normalized data, while the `copy number corrected` data remains unchanged.
#'   - Rarefies based on the calculated sampling depth, derived from the ratio of total reads per sample to
#'     predicted 16S rRNA gene copy numbers.
#'   - Uses the same `avgrarefy` function for averaging rarefied counts.
#'
#' - The input values are subject to a scaling limit of 1e7. If input values exceed this limit due to the prior normalization,
#'   all values are scaled down by a calculated scaling factor. The rarefied counts are rescaled back to their original scale
#'   after rarefaction. However, due to this scaling and the rounding of scaled values, slight variations in the rarefied counts
#'   may occur.
#'
#' The function uses parallel processing to improve the efficiency of rarefaction. It saves the rarefied `phyloseq`
#' object as an `.rds` file in the appropriate output folder.
#'
#' @return The rarefied `phyloseq` object is returned and also saved as an `.rds` file. The file name and location
#' depend on the specified normalization method:
#' - For `"fcm"`: `"project_name_phyloseq_asv_level_fcm_normalised_cell_concentration_rarefied.rds"`
#' - For `"qpcr"`: `"project_name_phyloseq_asv_level_qpcr_normalised_cell_concentration_rarefied.rds"`
#'
#' @examples
#' # Rarefy using FCM normalization
#' rarefied_physeq <- rarefying(physeq = normalised_physeq, norm_method = "fcm")
#'
#' # Rarefy using qPCR normalization
#' rarefied_physeq <- rarefying(physeq = normalised_physeq, norm_method = "qpcr")
#'
#' @note Ensure that the `phyloseq` object is properly normalized before applying this function.
#' Missing or invalid `rarefy_to` values will result in warnings and skipped samples.
#'
#' @export
rarefying = function(physeq = physeq,
                     norm_method = NULL,
                     iteration = 100) {

  project_name = projects

  project_folder = paste0(base_path, project_name)
  figure_folder = paste0(project_folder, "/figures/")
  output_folder_csv_files = paste0(project_folder, "/output_data/csv_files/")
  output_folder_rds_files = paste0(project_folder, "/output_data/rds_files/After_cleaning_rds_files/")
  output_asv_rds_files = paste0(output_folder_rds_files, "ASV/")
  if(!dir.exists(output_asv_rds_files)){dir.create(output_asv_rds_files)}

  avgrarefy = function(x, sample, iterations = iteration, seed = 711) {
    set.seed(seed)
    cl = makeCluster(detectCores() -4) # all cores minus two
    clusterEvalQ(cl, library(vegan)) # load vegan on all worker nodes

    tablist = parLapply(cl, seq_len(iterations), function(i) {
      rarefied_sample = suppressWarnings(rrarefy(x, sample = sample))
      return(rarefied_sample)
    })
    stopCluster(cl)

    afunc = array(unlist(tablist), c(dim(tablist[[1]]), iterations))
    output = apply(afunc, 1:2, mean)

    return(round(output, 0))
  }

  if (is.null(norm_method)) {
    psdata = physeq[["psdata_asv_copy_number_corrected"]]

    log_message(paste("Message: Normalization method is NULL. Returning unchanged phyloseq object."), log_file)

    return(psdata_asv_copy_number_corrected = psdata)

  } else if (norm_method == "fcm") {
    psdata = physeq[["psdata_asv_copy_number_corrected"]]
    psdata_fcm = physeq[["psdata_asv_fcm_norm"]]

    psdata_fcm <- prune_samples(sample_sums(psdata_fcm) > 0, psdata_fcm)        # Remove samples with zero counts
    psdata_fcm <- prune_taxa(rowSums(otu_table(psdata_fcm)) > 0, psdata_fcm)  # Remove taxa with zero counts across all samples

    # convert phyloseq to data frame and calculate cells per ml in the sample
    sample_data = data.frame(sample_data(psdata_fcm))

    # convert otu table to matrix
    ps_matrix = as(t(otu_table(psdata_fcm)), "matrix")

    # Determine the scaling factor based on the maxium value in the otu matrix
    max_sample_sum = max(sample_sums(psdata_fcm))
    limit = 1e7
    scaling_factor = 10 ^ ceiling(log10(max_sample_sum / limit))

    # scale down to OTU matrix and cells per ml
    scaled_ps_matrix = ceiling(ps_matrix / scaling_factor)
    sample_data$cells_per_ml = sample_data$cells_per_ml / scaling_factor

    # calculate sample sizes and rarefy depths
    sample_size = rowSums(scaled_ps_matrix) # total reads per sample
    cell_count_table = sample_data$cells_per_ml # extract cell counts
    sampling_depths = sample_size / cell_count_table # sampling depths (total reads divided by cell count)
    minimum_sampling_depth = min(sampling_depths) # minimum sampling depth across all samples
    rarefy_to = round(cell_count_table * minimum_sampling_depth, digits = 0) # number of reads to rarefy for each sample

    #psdata_phyloseq = otu_table(psdata_fcm, taxa_are_rows = FALSE)
    psdata_phyloseq = t(scaled_ps_matrix)

    # rarefy each sample based on the calculated rarefying targets (rarfy_to)
    rarefied_matrix = matrix(nrow = nrow(psdata_phyloseq),
                             ncol = ncol(psdata_phyloseq),
                             dimnames = list(rownames(psdata_phyloseq), colnames(psdata_phyloseq)))

    for (i in seq_len(ncol(psdata_phyloseq))) {
      sample_name = colnames(psdata_phyloseq)[i]
      sample_counts = psdata_phyloseq[, sample_name, drop = FALSE]

      if (!is.na(rarefy_to[i]) && rarefy_to[i] > 0) {
        rarefied_sample = avgrarefy(x = sample_counts, sample = rarefy_to[i], iterations = iteration, seed = 711)

        # Store rarefied counts in the corresponding column
        rarefied_matrix[, sample_name] <- as.numeric(rarefied_sample)
      } else {
        warning(paste("Skipping sample", sample_name, "due to invalid rarefy_to value"))
      }
    }

    # rescale the rarefied matrix back to the original scale
    rarefied_matrix = rarefied_matrix * scaling_factor

    rarefied_matrix_t = t(rarefied_matrix)
    colnames(rarefied_matrix_t) = taxa_names(psdata_fcm)
    rownames(rarefied_matrix_t) = sample_names(psdata_fcm)
    otu_rare = otu_table(rarefied_matrix_t, taxa_are_rows = FALSE)
    psdata_rarefied = psdata_fcm
    otu_table(psdata_rarefied) = otu_rare

    assign(paste0(project_name, "_rarefied_physeq"), psdata_rarefied, envir = .GlobalEnv)

    output_file_path = paste0(output_asv_rds_files, project_name, "_phyloseq_asv_level_fcm_normalised_cell_concentration_rarefied.rds")
    saveRDS(psdata_rarefied, file = output_file_path)
    log_message(paste("Phyloseq data fcm normalised cell concentration (cells per ml/gram sample) asv level rarefied saved as .rds object in", output_file_path), log_file)

    return(list(psdata_asv_copy_number_corrected = psdata, psdata_asv_fcm_norm_rarefied = psdata_rarefied))

  } else if (norm_method == "qpcr") {

    psdata_ccn = physeq[["psdata_asv_copy_number_corrected"]]
    psdata_qpcr = physeq[["psdata_asv_qpcr_norm"]]

    psdata_qpcr <- prune_samples(sample_sums(psdata_qpcr) > 0, psdata_qpcr)        # Remove samples with zero counts
    psdata_qpcr <- prune_taxa(rowSums(otu_table(psdata_qpcr)) > 0, psdata_qpcr)

    # convert phyloseq to data frame
    sample_data = data.frame(sample_data(psdata_qpcr))

    # convert otu table to matrix
    ps_matrix = as(t(otu_table(psdata_qpcr)), "matrix")

    # Determine the scaling factor based on the maxium value in the otu matrix
    max_sample_sum = max(sample_sums(psdata_qpcr))
    limit = 1e7
    scaling_factor = 10 ^ ceiling(log10(max_sample_sum / limit))

    # scale down to OTU matrix and cells per ml
    scaled_ps_matrix = ceiling(ps_matrix / scaling_factor)
    sample_data$sq_calc_mean = sample_data$sq_calc_mean / scaling_factor

    sample_size = rowSums(scaled_ps_matrix) # total reads per sample
    copy_count_table = sample_data$sq_calc_mean # extract copy counts???
    sampling_dephts = sample_size / copy_count_table # sampling dephts (total reads divided by copy counts??)
    minimum_sampling_depth = min(sampling_dephts) # minimum sampling depht across all samples
    rarefy_to = round(copy_count_table * minimum_sampling_depth, digits = 0) # number of read tot rarefy for each sample

    #psdata_phyloseq = otu_table(psdata_qpcr, taxa_are_rows = FALSE)
    psdata_phyloseq = t(scaled_ps_matrix)

    # rarefy each sample based on the calculated rarefying targets (rarfy_to)
    rarefied_matrix = matrix(nrow = nrow(psdata_phyloseq),
                             ncol = ncol(psdata_phyloseq),
                             dimnames = list(rownames(psdata_phyloseq), colnames(psdata_phyloseq)))

    for (i in seq_len(ncol(psdata_phyloseq))) {
      sample_name = colnames(psdata_phyloseq)[i]
      sample_counts = psdata_phyloseq[, sample_name, drop = FALSE]

      if (!is.na(rarefy_to[i]) && rarefy_to[i] > 0) {
        rarefied_sample = avgrarefy(x = sample_counts, sample = rarefy_to[i], iterations = iteration, seed = 711)

        # Store rarefied counts in the corresponding column
        rarefied_matrix[, sample_name] <- as.numeric(rarefied_sample)
      } else {
        warning(paste("Skipping sample", sample_name, "due to invalid rarefy_to value"))
      }
    }

    # rescale the rarefied matrix back to the original scale
    rarefied_matrix = rarefied_matrix * scaling_factor

    rarefied_matrix_t = t(rarefied_matrix)
    colnames(rarefied_matrix_t) = taxa_names(psdata_qpcr)
    rownames(rarefied_matrix_t) = sample_names(psdata_qpcr)
    otu_rare = otu_table(rarefied_matrix_t, taxa_are_rows = FALSE)
    psdata_rarefied = psdata_qpcr
    otu_table(psdata_rarefied) = otu_rare

    assign(paste0(project_name, "_rarefied_physeq"), psdata_rarefied, envir = .GlobalEnv)

    output_file_path = paste0(output_asv_rds_files, project_name, "_phyloseq_asv_level_qpcr_normalised_celL_concentration_rarefied.rds")
    saveRDS(psdata_rarefied, file = output_file_path)
    log_message(paste("phyloseq qpcr normalised cell concentration (cells per ml/gram sample) asv level rarefied saved as .rds object in", output_file_path), log_file)

    return(list(psdata_asv_copy_number_corrected = psdata_ccn, psdata_asv_qpcr_norm_rarefied = psdata_rarefied))

  } else {
    log_message(paste("Error: Invalid normalization method specified. Use 'fcm' or 'qpcr'."), log_file)
    stop("Error: Invalid normalization method specified. Use 'fcm' or 'qpcr'.")
  }
}
