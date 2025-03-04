#' Convert Phyloseq Data to Tibble and Export by Taxonomic Level
#'
#' This function converts a phyloseq object into a tibble for each specified taxonomic level by using the \code{psmelt} function. The resulting tibbles are then saved as RDS files in a designated output folder. Depending on the provided normalization method, the function processes:
#' \itemize{
#'   \item \code{NULL} (default): Exports only the copy number corrected counts.
#'   \item \code{"fcm"}: Exports both the copy number corrected counts and flow cytometry-normalized counts.
#'   \item \code{"qpcr"}: Exports both the copy number corrected counts and qPCR-normalized counts.
#' }
#'
#' @param physeq A phyloseq object containing genus-level data. The default is \code{rarefied_genus_physeq}.
#' @param norm_method A character string specifying the normalization method to use. Options include:
#'   \itemize{
#'     \item \code{NULL} (default): Process and export copy number corrected counts.
#'     \item \code{"fcm"}: Process and export both copy number corrected counts and flow cytometry-normalized counts.
#'     \item \code{"qpcr"}: Process and export both copy number corrected counts and qPCR-normalized counts.
#'   }
#' @param taxrank A character vector indicating the taxonomic levels to process. The default is \code{c("Phylum", "Class", "Order", "Family", "Tax_label")}.
#'
#' @details
#' For each taxonomic level specified in \code{taxrank}, the function performs the following steps:
#' \enumerate{
#'   \item Logs a message indicating the taxonomic level currently being processed.
#'   \item Creates a folder for that taxonomic level under the project's \code{After_cleaning_rds_files} directory if it does not already exist.
#'   \item Converts the phyloseq data (or a subset thereof) into a tibble via \code{psmelt()} and \code{as_tibble()}.
#'   \item Saves the resulting tibble as an RDS file with a filename that includes the project name, taxonomic level, and normalization method.
#'   \item Stores each exported tibble in a list with names that reflect the data type (e.g., \code{"psmelt_copy_number_corrected_Phylum"}, \code{"psmelt_fcm_norm_rarefied_Class"}, etc.).
#' }
#'
#' @return A list containing the psmelt tibble data for each taxonomic level. The list names indicate both the processing method and the taxonomic level.
#'
#' @examples
#' \dontrun{
#'   # Export data without normalization (copy number corrected counts only)
#'   result <- psdata_to_tibble(physeq = rarefied_genus_physeq)
#'
#'   # Export data with flow cytometry normalization
#'   result <- psdata_to_tibble(physeq = rarefied_genus_physeq, norm_method = "fcm")
#'
#'   # Export data with qPCR normalization
#'   result <- psdata_to_tibble(physeq = rarefied_genus_physeq, norm_method = "qpcr")
#' }
#'
#' @export
psdata_to_tibble = function(physeq = rarefied_genus_physeq,
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
      psdata_tibble = psdata %>% psmelt() %>% as_tibble()

      output_file_path = paste0(tax_folder, "/", project_name, "_psmelt_", tax, "_level_copy_number_corrected_counts.rds")
      saveRDS(psdata_tibble, file = output_file_path)
      log_message(paste("psmelt data copy number corrected counts", tax, "level saved as .rds object in", output_file_path), log_file)

      results[[paste0("psmelt_copy_number_corrected_", tax)]] <- psdata_tibble
    } else if (norm_method == "fcm") {
      psdata = physeq[[paste0("psdata_copy_number_corrected_", tax)]]
      psdata_tibble = psdata %>% psmelt() %>% as_tibble()

      output_file_path = paste0(tax_folder, "/", project_name, "_psmelt_", tax, "_level_copy_number_corrected_counts.rds")
      saveRDS(psdata_tibble, file = output_file_path)
      log_message(paste("psmelt data copy number corrected counts", tax, "level saved as .rds object in", output_file_path), log_file)

      psdata_fcm = physeq[[paste0("psdata_fcm_norm_rarefied_", tax)]]
      psdata_tibble_fcm = psdata_fcm %>% psmelt() %>% as_tibble()

      output_file_path = paste0(tax_folder, "/", project_name, "_psmelt_", tax, "_level_fcm_normalised_cell_concentration_rarefied.rds")
      saveRDS(psdata_tibble_fcm, file = output_file_path)
      log_message(paste("psmelt data fcm normalised cell concentration (cells per ml/gram sample)", tax, "level rarefied saved as .rds object in", output_file_path), log_file)

      results[[paste0("psmelt_copy_number_corrected_", tax)]] <- psdata_tibble
      results[[paste0("psmelt_fcm_norm_rarefied_", tax)]] <- psdata_tibble_fcm
    } else if (norm_method == "qpcr") {
      psdata = physeq[[paste0("psdata_copy_number_corrected_", tax)]]
      psdata_tibble = psdata %>% psmelt() %>% as_tibble()

      output_file_path = paste0(tax_folder, "/", project_name, "_psmelt_", tax, "_level_copy_number_corrected_counts.rds")
      saveRDS(psdata_tibble, file = output_file_path)
      log_message(paste("psmelt data copy number corrected counts", tax, "level saved as .rds object in", output_file_path), log_file)

      psdata_qpcr = physeq[[paste0("psdata_qpcr_norm_rarefied_", tax)]]
      psdata_tibble_qpcr = psdata_qpcr %>% psmelt() %>% as_tibble()

      output_file_path = paste0(tax_folder, "/", project_name, "_psmelt_", tax, "_level_qpcr_normalised_celL_concentration_rarefied.rds")
      saveRDS(psdata_tibble_qpcr, file = output_file_path)
      log_message(paste("psmelt qpcr normalised cell concentration (cells per ml/gram sample)", tax, "level rarefied saved as .rds object in", output_file_path), log_file)

      results[[paste0("psmelt_copy_number_corrected_", tax)]] <- psdata_tibble
      results[[paste0("psmelt_qpcr_norm_rarefied_", tax)]] <- psdata_tibble_qpcr
    }
  }
  return(results)
}
