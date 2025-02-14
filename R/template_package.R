#' Function Title (Template)
#'
#' Short description of what the function does.
#'
#' @param projects A character vector containing the names of the projects.
#'                 Specify what this parameter represents and how it should be structured.
#' @param base_path (Optional) A string specifying the base path where project directories are located.
#'                  If not provided, a default path is assumed.
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#' \item Step 1: Explain the first operation performed by the function.
#' \item Step 2: Provide details about the next step in the function.
#' \item Step 3: Continue describing all relevant steps in the process.
#' }
#'
#' Key assumptions or specific requirements:
#' \itemize{
#' \item Ensure that required files exist in the expected directories.
#' \item Metadata must contain the `SampleID` column.
#' \item Describe any other key considerations here.
#' }
#'
#' @return
#' Describe the output of the function. For example:
#' A data frame containing processed and merged metadata, or a logical value indicating success.
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic usage
#' unified_metadata <- unify_metadata(projects = c("Project1", "Project2"))
#'
#' # Example 2: Using a custom base path
#' unified_metadata <- unify_metadata(projects = "Project1", base_path = "/custom/path/")
#' }
#'
#' @seealso
#' Other related functions or links to external resources.
#'
#' @export
