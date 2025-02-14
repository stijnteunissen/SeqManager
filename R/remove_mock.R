#' Remove Mock ASVs from Phyloseq Object
#'
#' This function removes mock ASVs (amplicon sequence variants) from a `phyloseq` object.
#' Mock ASVs are typically used for decontamination purposes and are identified by their genus names.
#' The function also filters the dataset to retain only samples with `sample_type == "sample`.
#' The user can choose whether to filter out the mock ASVs from the dataset by setting the `mock` parameter.
#'
#' @inheritParams decontam
#'
#' @param mock_genera A vector of genera representing the taxa that make up the mock community.
#'                    These are the taxa used to identify mock ASVs that can be removed from the `phyloseq` object.
#'
#' @param mock A logical value (default `TRUE`) that determines whether or not to filter out the mock ASVs.
#'             \itemize{
#'             \item `TRUE`: Remove mock ASVs from the `phyloseq` object and retain only the samples for downstream analysis.
#'             \item `FALSE`: Retain mock ASVs but filter the dataset to include only samples with `sample_type == "sample"`.
#'             Use this option if no mock community is present in the dataset.
#'             }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item If `mock = FALSE`, the function filters the dataset to retain only samples with `sample_type == "sample"`, without removing mock ASVs.
#'         This option is suitable for datasets where no mock community is included.
#'   \item If `mock = TRUE`, the function:
#'     \itemize{
#'       \item Identifies mock ASVs based on the provided `mock_genera`.
#'       \item Removes the mock ASVs from the dataset.
#'       \item Retains only samples with `sample_type == "sample"`.
#'     }
#'   \item The filtered `phyloseq` object is saved as an RDS file with the name `project_name_phyloseq_asv_level_without_mock.rds`.
#' }
#'
#' @return
#' A `phyloseq` object with filtered samples.
#' If `mock = TRUE`, mock ASVs are removed. If `mock = FALSE`, mock ASVs are retained, but the dataset includes only `sample_type == "sample"`.
#' The object is also saved as an RDS file for future use.
#'
#' @examples
#' \dontrun{
#' # Remove mock ASVs from the phyloseq object
#' physeq_no_mock <- remove_mock(physeq = physeq, mock_genera = c("Mock_Genus1", "Mock_Genus2"), mock = TRUE)
#'
#' # Retain mock ASVs but filter to only 'sample' types (or when no mock community is present)
#' physeq_no_filter <- remove_mock(physeq = physeq, mock_genera = c("Mock_Genus1", "Mock_Genus2"), mock = FALSE)
#' }
#'
#' @export
remove_mock = function(physeq = decontam_physeq,
                       mock_genera = mock_genera,
                       mock = TRUE) {

  psdata = physeq
  project_name = projects

  project_folder = paste0(base_path, project_name)
  output_folder_rds_files = paste0(project_folder, "/output_data/rds_files/Before_cleaning_rds_files/")

  if (mock == FALSE) {
    physeq_filtered = psdata
    message = paste0("message: mock parameter is FALSE. mock ASVs not filtered")
    log_message(message, log_file)
  } else if (mock == TRUE) {

    # if (!is.null(extra_genera)) {
    #   mock_genera = c(mock_genera, extra_genera)
    # } else
    #   mock_genera = mock_genera

    # extract mock OTU mock_genera
    mock_ASVs =
      psdata %>%
      subset_taxa(Genus %in% mock_genera) %>%
      taxa_names()

    # Filter phyloseq to remove mock samples
    physeq_no_mock =
      psdata %>%
      subset_samples(., sample_or_control == "sample") %>%
      prune_taxa(taxa_sums(.) > 0, .)

    physeq_filtered =
      physeq_no_mock %>%
      prune_taxa(!taxa_names(.) %in% mock_ASVs, .)
  }

  output_file_path = paste0(output_folder_rds_files, project_name, "_phyloseq_asv_level_without_mock.rds")
  saveRDS(physeq_filtered, file = output_file_path)
  log_message(paste("Phyloseq object without mock saved as .rds object in", output_file_path), log_file)

  return(physeq_filtered)
}
