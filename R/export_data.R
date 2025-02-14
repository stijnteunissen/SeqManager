#' Export Project Data to a Zip Archive
#'
#' This function exports project-related figures and output data by creating a dedicated export folder within the project directory. It organizes the export folder into subdirectories for figures, CSV files, and RDS files, copies the relevant files from their original locations, and finally compresses the export folder into a zip archive.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Constructs paths for the project folder, figures, and output data using global variables (e.g., \code{projects} and \code{base_path}).
#'   \item Creates an export folder named \code{export_data_<project_name>} within the project folder.
#'   \item Creates subdirectories within the export folder for figures, CSV files, and RDS files.
#'   \item Copies figures from the project's figure folder to the export folder.
#'   \item Copies CSV files from the output CSV folder to the corresponding subfolder in the export directory.
#'   \item Copies RDS files from the \code{After_cleaning_rds_files} subfolder (if present) to the export folder.
#'   \item Zips the contents of the export folder into a zip archive containing only the relative file paths.
#' }
#'
#' @return This function does not return a value. It performs file operations and creates a zip archive in the project folder.
#'
#' @examples
#' \dontrun{
#'   # Export project data to a zip archive
#'   export_data()
#' }
#'
#' @export
export_data <- function() {

  project_name = projects
  project_folder = paste0(base_path, project_name)
  figure_folder = paste0(project_folder, "/figures/")
  output_folder = paste0(project_folder, "/output_data/")
  output_folder_csv_files = paste0(project_folder, "/output_data/csv_files/")
  output_folder_rds_files = paste0(project_folder, "/output_data/rds_files/")

  # create export folder
  if(!dir.exists(paste0(project_folder, "/export_data_", project_name))){dir.create(paste0(project_folder, "/export_data_", project_name))}

  # path to export folder
  export_folder = paste0(project_folder, "/export_data_", project_name)

  # create sub folders in export folder and path to sub folders
  if(!dir.exists(paste0(export_folder, "/figures"))){dir.create(paste0(export_folder, "/figures"))}
  if(!dir.exists(paste0(export_folder, "/output_data"))){dir.create(paste0(export_folder, "/output_data"))}
  if(!dir.exists(paste0(export_folder, "/output_data/csv_files"))){dir.create(paste0(export_folder, "/output_data/csv_files"))}
  if(!dir.exists(paste0(export_folder, "/output_data/rds_files"))){dir.create(paste0(export_folder, "/output_data/rds_files"))}

  export_figures = paste0(export_folder, "/figures/")
  export_output_folder = paste0(export_folder, "/output_data/")
  export_csv_output_folder = paste0(export_folder, "/output_data/csv_files/")
  export_rds_output_folder = paste0(export_folder, "/output_data/rds_files/")

  # copy figures naar export folder
  if (dir.exists(figure_folder)) {
    file.copy(from = list.files(figure_folder, full.names = TRUE),
              to = export_figures, recursive = TRUE, overwrite = TRUE)
  }

  # Copy the entire CSV output folder to the export folder
  if (dir.exists(output_folder_csv_files)) {
    file.copy(from = list.files(output_folder_csv_files, full.names = TRUE),
              to = export_csv_output_folder,
              recursive = TRUE,
              overwrite = TRUE)
  }

  # Copy only the "After_cleaning_rds_files" folder from the RDS output folder.
  after_cleaning_folder <- file.path(output_folder_rds_files, "After_cleaning_rds_files")
  if (dir.exists(after_cleaning_folder)) {
    file.copy(from = list.files(after_cleaning_folder, full.names = TRUE),
              to = export_rds_output_folder,
              recursive = TRUE,
              overwrite = TRUE)
  }

  # # Define RDS files and copy them
  # rds_files <- c(
  #   paste0(output_folder_rds_files, project_name, "_phyloseq_asv_level_anna16_corrected_counts.rds"),
  #   paste0(output_folder_rds_files, project_name, "_phyloseq_asv_level_fcm_normalised_cell_concentration_rarefied.rds"),
  #   paste0(output_folder_rds_files, project_name, "_phyloseq_asv_level_qpcr_normalised_celL_concentration_rarefied.rds")
  # )
  #
  # # Copy RDS files if they exist
  # for (rds_file in rds_files) {
  #   if (file.exists(rds_file)) {
  #     file.copy(rds_file, export_rds_output_folder, overwrite = TRUE)
  #   }
  # }

  # # zip export folder
  # zip_file <- paste0(export_folder, ".zip")
  # files_to_zip <- list.files(export_folder, full.names = TRUE, recursive = TRUE)
  # zip(zipfile = zip_file, files = files_to_zip)

  # Zip the export folder but include only the relative paths.
  # Save the current working directory.
  old_wd <- getwd()
  # Change working directory to the export folder so that only its contents are zipped.
  setwd(export_folder)
  files_to_zip <- list.files(".", recursive = TRUE)
  zip_file <- paste0(export_folder, ".zip")
  zip(zipfile = zip_file, files = files_to_zip)
  # Reset working directory.
  setwd(old_wd)

}
