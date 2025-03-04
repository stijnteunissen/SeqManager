#' Aggregate Phyloseq Data by Taxonomic Level
#'
#' This function aggregates a phyloseq object to higher taxonomic levels (e.g.,
#' Phylum, Class, Order, Family, or Tax_label) using the \code{tax_glom}
#' function. Depending on the specified normalization method, the function
#' processes and saves both copy number corrected data and normalized data (using
#' flow cytometry or qPCR) as RDS files in dedicated directories.
#'
#' @param physeq A phyloseq object containing ASV-level data. For example, \code{rarefied_asv_physeq}.
#' @param norm_method A character string specifying the normalization method. Options are:
#'   \itemize{
#'     \item \code{NULL} (default): Process copy number corrected counts only.
#'     \item \code{"fcm"}: Process both copy number corrected counts and flow cytometry-normalized data.
#'     \item \code{"qpcr"}: Process both copy number corrected counts and qPCR-normalized data.
#'   }
#' @param taxrank A character vector indicating the taxonomic levels at which to aggregate the data. The default is \code{c("Phylum", "Class", "Order", "Family", "Tax_label")}.
#'
#' @details
#' For each taxonomic level specified in \code{taxrank}, the function:
#' \enumerate{
#'   \item Logs a message indicating the taxonomic level being processed.
#'   \item Creates a subfolder (if it does not already exist) under the project's \code{After_cleaning_rds_files} directory.
#'   \item Aggregates the data at the specified taxonomic level using \code{tax_glom}.
#'   \item Saves the aggregated phyloseq object as an RDS file. The filename and folder structure reflect the taxonomic level and the normalization method used.
#'   \item If a normalization method is provided (\code{"fcm"} or \code{"qpcr"}), both the copy number corrected data and the corresponding normalized data are processed and saved.
#' }
#'
#' @return A list of aggregated phyloseq objects. The list elements are named
#'   according to the data type and taxonomic level (e.g.,
#'   \code{"psdata_copy_number_corrected_Phylum"},
#'   \code{"psdata_fcm_norm_rarefied_Class"}, etc.).
#'
#' @examples
#' \dontrun{
#'   # Aggregate data at default taxonomic levels using only copy number corrected counts
#'   result <- group_tax(physeq = rarefied_asv_physeq)
#'
#'   # Aggregate data using flow cytometry normalization
#'   result <- group_tax(physeq = rarefied_asv_physeq, norm_method = "fcm")
#'
#'   # Aggregate data using qPCR normalization
#'   result <- group_tax(physeq = rarefied_asv_physeq, norm_method = "qpcr")
#' }
#'
#' @export
group_tax = function(physeq = rarefied_asv_physeq,
                     norm_method = NULL,
                     taxrank = c("Phylum", "Class", "Order", "Family", "Tax_label")) {

  project_name = projects
  project_folder = paste0(base_path, project_name)
  output_folder_rds_files = paste0(project_folder, "/output_data/rds_files/After_cleaning_rds_files/")

  results = list()

  for (tax in taxrank) {
    log_message(paste("Processing taxonomic level:", tax), log_file)
    tax_folder = file.path(output_folder_rds_files, paste0(tax))
    if(!dir.exists(tax_folder)){dir.create(tax_folder)}

    if (is.null(norm_method)) {
      psdata = physeq
      psdata_tax = tax_glom(psdata, taxrank = tax)

      output_file_path = paste0(tax_folder, "/", project_name, "_phyloseq_", tax, "_level_copy_number_corrected_counts.rds")
      saveRDS(psdata_tax, file = output_file_path)
      log_message(paste("phyloseq data copy number corrected counts", tax, "level saved as .rds object in", output_file_path), log_file)

      results[[paste0("psdata_copy_number_corrected_", tax)]] <- psdata_tax
    } else if (norm_method == "fcm") {
      psdata = physeq[["psdata_asv_copy_number_corrected"]]
      psdata_tax = tax_glom(psdata, taxrank = tax)

      output_file_path = paste0(tax_folder, "/", project_name, "_phyloseq_", tax, "_level_copy_number_corrected_counts.rds")
      saveRDS(psdata_tax, file = output_file_path)
      log_message(paste("phyloseq data copy number corrected counts", tax, "level saved as .rds object in", output_file_path), log_file)

      psdata_fcm = physeq[["psdata_asv_fcm_norm_rarefied"]]
      psdata_tax_fcm = tax_glom(psdata_fcm, taxrank = tax)

      output_file_path = paste0(tax_folder, "/", project_name, "_phyloseq_", tax, "_level_fcm_normalised_cell_concentration_rarefied.rds")
      saveRDS(psdata_tax_fcm, file = output_file_path)
      log_message(paste("Phyloseq data fcm normalised cell concentration (cells per ml/gram sample)", tax, " level rarefied saved as .rds object in", output_file_path), log_file)

      results[[paste0("psdata_copy_number_corrected_", tax)]] <- psdata_tax
      results[[paste0("psdata_fcm_norm_rarefied_", tax)]] <- psdata_tax_fcm
    } else if (norm_method == "qpcr") {
      psdata = physeq[["psdata_asv_copy_number_corrected"]]
      psdata_tax = tax_glom(psdata, taxrank = tax)

      output_file_path = paste0(tax_folder, "/", project_name, "_phyloseq_", tax, "_level_copy_number_corrected_counts.rds")
      saveRDS(psdata_tax, file = output_file_path)
      log_message(paste("phyloseq data copy_number corrected counts", tax, " level saved as .rds object in", output_file_path), log_file)

      psdata_qpcr = physeq[["psdata_asv_qpcr_norm_rarefied"]]
      psdata_tax_qpcr = tax_glom(psdata_qpcr, taxrank = tax)

      output_file_path = paste0(tax_folder, "/", project_name, "_phyloseq_", tax, "_level_qpcr_normalised_cell_concentration_rarefied.rds")
      saveRDS(psdata_tax_qpcr, file = output_file_path)
      log_message(paste("phyloseq qpcr normalised cell concentration (cells per ml/gram sample)", tax, " level rarefied saved as .rds object in", output_file_path), log_file)

      results[[paste0("psdata_copy_number_corrected_", tax)]] <- psdata_tax
      results[[paste0("psdata_qpcr_norm_rarefied_", tax)]] <- psdata_tax_qpcr
    }
  }

  return(results)
}
