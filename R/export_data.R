#' Export Project Data to a Zip Archive
#'
#' This function exports project-related figures and output data by creating a dedicated export folder
#' within the project directory. It organizes the export folder into subdirectories for figures, CSV files,
#' and RDS files, copies the relevant files from their original locations, and finally compresses the export
#' folder into a zip archive containing only relative file paths.
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
  # Check that the required global variables exist
  if (!exists("base_path") || !exists("projects")) {
    stop("Ensure that 'base_path' and 'projects' are defined.")
  }

  project_name <- projects
  project_folder <- paste0(base_path, project_name)
  figure_folder <- file.path(project_folder, "figures")
  output_folder <- file.path(project_folder, "output_data")
  output_folder_csv_files <- file.path(output_folder, "csv_files")
  output_folder_rds_files <- file.path(output_folder, "rds_files")

  # Create the export folder
  export_folder <- file.path(project_folder, paste0("export_data_", project_name))
  if (!dir.exists(export_folder)) {
    dir.create(export_folder, recursive = TRUE)
  }

  # Create subdirectories within the export folder
  export_figures <- file.path(export_folder, "figures")
  export_csv_output_folder <- file.path(export_folder, "output_data", "csv_files")
  export_rds_output_folder <- file.path(export_folder, "output_data", "rds_files")

  if (!dir.exists(export_figures)) {
    dir.create(export_figures, recursive = TRUE)
  }
  if (!dir.exists(export_csv_output_folder)) {
    dir.create(export_csv_output_folder, recursive = TRUE)
  }
  if (!dir.exists(export_rds_output_folder)) {
    dir.create(export_rds_output_folder, recursive = TRUE)
  }

  # Copy figures to the export folder
  if (dir.exists(figure_folder)) {
    file.copy(from = list.files(figure_folder, full.names = TRUE),
              to = export_figures, recursive = TRUE, overwrite = TRUE)
  }

  # Copy CSV files to the export folder
  if (dir.exists(output_folder_csv_files)) {
    file.copy(from = list.files(output_folder_csv_files, full.names = TRUE),
              to = export_csv_output_folder, recursive = TRUE, overwrite = TRUE)
  }

  # Copy RDS files from the "After_cleaning_rds_files" folder
  after_cleaning_folder <- file.path(output_folder_rds_files, "After_cleaning_rds_files")
  if (dir.exists(after_cleaning_folder)) {
    file.copy(from = list.files(after_cleaning_folder, full.names = TRUE),
              to = export_rds_output_folder, recursive = TRUE, overwrite = TRUE)
  }


  # # Zip the export folder contents using relative paths.
  # # Save the current working directory.
  # old_wd <- getwd()
  # on.exit(setwd(old_wd))  # Ensure the old working directory is restored even if an error occurs.
  #
  # # Change working directory to the export folder so that only its contents are zipped.
  # setwd(export_folder)
  #
  # files_to_zip <- list.files(".", recursive = TRUE)
  # zip_file <- paste0(export_folder, "_data.zip")
  #
  # # Check that the zip_file path is not empty.
  # if (nchar(zip_file) == 0) {
  #   stop("The zip file path is empty!")
  # }
  #
  # # Warn if no files were found to zip.
  # if (length(files_to_zip) == 0) {
  #   warning("No files found in the export folder to zip.")
  # }
  #
  # # Create the zip archive.
  # utils::zip(zipfile = zip_file, files = files_to_zip)
  #
  # message("Export completed! Zip file created: ", zip_file)
}
