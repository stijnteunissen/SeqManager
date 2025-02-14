#' Create Project Folders and Process Files
#'
#' This function creates the necessary directory structure for specified projects,
#' checks for required and optional files in the source directory, and copies
#' relevant files to the appropriate locations.
#'
#' @param projecten A character vector containing the names of projects. Each project will have its own folder structure created.
#' @return No return value. The function creates directories and copies files as a side effect.
#' @details
#' The function performs the following steps for each project:
#' \itemize{
#' \item Creates directories if they do not already exist.
#' \item Checks for required files (`table.qza`, `rooted-tree.qza`, `classifier.qza`, `metadata.tsv`).
#' \item Logs errors if required files are missing and stops execution.
#' \item Checks for optional files (`dna-sequences.csv`, `qPCR.csv`, `fcm.csv`, `metadata_extra.tsv`).
#' \item Logs warnings if optional files are missing.
#' \item Copies relevant files to the `input_R_data` directory.
#' }
#'
#' @export
#' @examples
#' # Example usage:
#' projecten <- c("Project1", "Project2")
#' creating_folders(projecten)
creating_folders <- function(projecten) {

  for (project in projecten) {

    project_folder <- paste0(base_path, project)

    # Create the required directories if they don't exist
    dir.create(paste0(project_folder, "/input_R_data"), recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(project_folder, "/output_R_data/csv_files"), recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(project_folder, "/output_R_data/rds_files"), recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(project_folder, "/messages"), recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(project_folder, "/figures"), recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(project_folder, "/scripts"), recursive = TRUE, showWarnings = FALSE)

    source_folder <- paste0(project_folder, "/output")
    destination_folder <- paste0(project_folder, "/input_R_data")

    # Messages file for warnings and errors
    messages_file <- paste0(project_folder, "/messages/", project, "_warnings.txt")

    # Remove the messages file if it exists (start fresh)
    if (file.exists(messages_file)) {
      file.remove(messages_file)
    }

    files <- list.files(source_folder, full.names = TRUE)

    # Check for required files - stop if missing
    required_files <- c("table.*\\.qza$", "rooted-tree.*\\.qza$", "classifier.*\\.qza", "metadata\\.tsv$")
    for (file_pattern in required_files) {
      if (!any(grepl(file_pattern, files))) {
        error_message <- paste("Error:", file_pattern, "does not exist in", source_folder, "for project:", project, "\n")
        cat(error_message, file = messages_file, append = TRUE)
        stop(error_message)
      }
    }

    # Check for optional files - only warning if missing
    optional_files <- c("dna-sequences.*\\.csv$", "qPCR.*\\.csv$", "fcm.*\\.csv$", "metadata_extra\\.tsv$")
    for (file_pattern in optional_files) {
      if (!any(grepl(file_pattern, files))) {
        warning_message <- paste("Warning:", file_pattern, "does not exist in", source_folder, "for project:", project, "\n")
        cat(warning_message, file = messages_file, append = TRUE)
      }
    }

    # Copy the relevant files
    files_for_phyloseq_object <- list.files(source_folder, pattern = "table.*\\.qza$|rooted-tree.*\\.qza$|classifier.*\\.qza|metadata\\.tsv$|metadata_extra\\.tsv$|dna-sequences.*\\.csv$|fcm.*\\.csv$|qPCR.*\\.csv$", full.names = TRUE)

    file.copy(files_for_phyloseq_object, destination_folder, overwrite = TRUE)
  }
}
