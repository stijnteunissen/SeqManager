#' Create Project Folder Structure
#'
#' This function creates a folder structure for projects, ensuring that
#' all necessary directories exist and that specific files required for downstream
#' analysis are copied from the `qiime2_output` folder into the `input_data` folder.
#'
#' @param projects A character vector containing the name of the project folders.
#'                 Each project folder must already exist within the specified
#'                 base path. The `projects` argument represents the name of the
#'                 project folder, which must contain a subfolder named `qiime2_output`
#'                 with the output files generated from QIIME2. This folder structure
#'                 is necessary for setting up downstream analysis.
#'                 `projects` is typically reference in a configuration file.
#' @param base_path A character string indicating the base directory where
#'                  project folder is located. This path is used to locate
#'                  the project folder and create the required subfolders.
#'                  `base_path` is typically referenced in a configuration file.
#'
#' @details
#' This function is designed to facilitate the setup of downstream analyses by:
#' \itemize{
#' \item Creating a consistent directory structure for each project, including
#'       subfolders like `input_data`, `output_data`, `figures`, and `messages`.
#' \item Copying essential files from the `qiime2_output` folder into the `input_data` folder.
#'       These essential files include:
#'       \itemize{
#'       \item `table.qza`: The feature table from QIIME2.
#'       \item `rooted-tree.qza`: The phylogenetic tree used for diversity analysis.
#'       \item `classifier.qza`: The classifier used for taxonomy assignment.
#'       \item `metadata.tsv`: The QIIME2 sample metadata file required for analyses.
#'       \item `dna-sequences.csv`: Contains 16S rRNA gene copy number predictions.
#'       \item `metadata_extra.tsv`: Any additional metadata provided for the samples.
#'       }
#' \item Checking for optional files that enhance analyses, such as:
#'       \itemize{
#'       \item `qPCR.csv`: Contains quantitative PCR data.
#'       \item `fcm.csv`: Contains flow cytometry data.
#'       }
#' \item Logging warnings for missing optional files and errors for missing required files.
#' }
#'
#' The function ensures that downstream analyses, which rely on specific input files
#' (e.g., the `table.qza`, `rooted-tree.qza`, etc.), have access to these files in the
#' correct directory structure. If any required files are missing from the
#' `qiime2_output` folder, the function stops execution and logs an error message.
#'
#' @note
#' Each project folder must already exist within the `base_path` directory and must
#' contain a subfolder named `qiime2_output`, which holds the outputs of QIIME2 analysis.
#' The function sets up the folder structure for downstream analysis within this project folder.
#'
#' @return None. This function is called for its side effects.
#'
#' @examples
#' \dontrun{
#' # Define the base path and log file location
#' base_path <- "path/to/projects"
#' log_file <- "path/to/logfile.txt"
#'
#' # Create folder structures for two projects
#' create_folders(projects)
#' }
#'
#' @export
create_folders = function(projects) {

  for (project in projects) {

    project_folder = paste0(base_path, project)

    # Create the required directories if they don't exist
    if(!dir.exists(paste0(project_folder, "/input_data"))){dir.create(paste0(project_folder, "/input_data"))}
    if(!dir.exists(paste0(project_folder, "/output_data"))){dir.create(paste0(project_folder, "/output_data"))}
    if(!dir.exists(paste0(project_folder, "/output_data/csv_files"))){dir.create(paste0(project_folder, "/output_data/csv_files"))}
    if(!dir.exists(paste0(project_folder, "/output_data/rds_files"))){dir.create(paste0(project_folder, "/output_data/rds_files"))}
    if(!dir.exists(paste0(project_folder, "/output_data/rds_files/Before_cleaning_rds_files"))){dir.create(paste0(project_folder, "/output_data/rds_files/Before_cleaning_rds_files"))}
    if(!dir.exists(paste0(project_folder, "/output_data/rds_files/After_cleaning_rds_files"))){dir.create(paste0(project_folder, "/output_data/rds_files/After_cleaning_rds_files"))}
    if(!dir.exists(paste0(project_folder, "/messages"))){dir.create(paste0(project_folder, "/messages"))}
    if(!dir.exists(paste0(project_folder, "/figures"))){dir.create(paste0(project_folder, "/figures"))}

    source_folder = paste0(project_folder, "/qiime2_output")
    destination_folder = paste0(project_folder, "/input_data")

    files = list.files(source_folder, full.names = TRUE)

    # Check for required files - stop if missing
    required_files = c("table.*\\.qza$", "rooted-tree.*\\.qza$", "classifier.*\\.qza", "metadata\\.tsv$")
    for (file_pattern in required_files) {
      if (!any(grepl(file_pattern, files))) {  # Use 'files' instead of 'file'
        error_message = paste("Error:", file_pattern, "does not exist in", source_folder, "for project:", project, "\n")
        log_message(error_message, log_file)
        stop(error_message)
      }
    }

    # Check for optional files - only warning if missing
    optional_files = c("dna-sequences.*\\.csv$", "qPCR.*\\.csv$", "fcm.*\\.csv$", "metadata_extra\\.tsv$")
    for (file_pattern in optional_files) {
      if (!any(grepl(file_pattern, files))) {  # Use 'files' instead of 'file'
        warning_message = paste("Warning:", file_pattern, "does not exist in", source_folder, "for project:", project, "\n")
        log_message(warning_message, log_file)
      }
    }

    # Copy the relevant files
    files_for_phyloseq_object =
      list.files(source_folder,
                 pattern = "table.*\\.qza$|rooted-tree.*\\.qza$|classifier.*\\.qza|metadata\\.tsv$|metadata_extra\\.tsv$|dna-sequences.*\\.csv$|fcm.*\\.csv$|qPCR.*\\.csv$",
                 full.names = TRUE)

    file.copy(files_for_phyloseq_object, destination_folder, overwrite = TRUE)
  }
}
