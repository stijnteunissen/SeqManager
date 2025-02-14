#' Resolve Multichotomies in Phyloseq Object's Tree
#'
#' This function processes a `phyloseq` object to resolve multichotomies (non-binary nodes) in the phylogenetic tree
#' using the `ape` package. It ensures the tree is binary and saves the updated `phyloseq` object as an RDS file.
#'
#' @inheritParams tax_clean
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Checks if the `phyloseq` object contains a binary phylogenetic tree.
#'   \item Resolves polychotomous nodes (if present) using `ape::multi2di()`.
#'   \item Ensures that the tree is binary after resolution; raises an error if unresolved.
#'   \item Merges the resolved tree back with the `otu_table`, `sample_data`, and `tax_table` in the `phyloseq` object.
#'   \item Saves the updated `phyloseq` object with the resolved tree as an RDS file.
#' }
#'
#' @return A `phyloseq` object with a binary phylogenetic tree.
#'
#' @examples
#' \dontrun{
#' # Resolve tree and save the updated phyloseq object
#' resolved_physeq <- resolve_tree(physeq = cleaned_physeq)
#' }
#'
#' @export
# Define function to process physeq objects for resolving multichotomies using Ape
# Define function to process physeq objects for resolving multichotomies using Ape
resolve_tree = function(physeq = cleaned_physeq) {

  psdata = physeq
  project_name = projects

  project_folder = paste0(base_path, projects)
  output_folder_csv_files = paste0(project_folder, "/output_data/csv_files/")
  output_folder_rds_files = paste0(project_folder, "/output_data/rds_files/Before_cleaning_rds_files/")

  if (!inherits(psdata, "phyloseq")) {
    stop("Error: psdata is not a phyloseq object.")
  }

  # Check if the tree is binary
  if (!ape::is.binary(phy_tree(psdata))) {
    # Resolve polychotomous nodes
    phy_tree_resolved <- ape::multi2di(phy_tree(psdata))
    # Check if resolved
    if (!ape::is.binary(phy_tree_resolved)) {
      stop("Error: Unable to resolve polychotomous nodes.")
    }
    # Update tree
    tree2 <- phy_tree_resolved
  } else {
    # Use the original tree if it's already binary
    tree2 <- phyloseq::phy_tree(psdata)
  }

  # Merge new tree with sample_data and otu_table
  new_tree <- phyloseq::merge_phyloseq(
    phyloseq::otu_table(psdata),
    phyloseq::sample_data(psdata),
    phyloseq::tax_table(psdata),
    tree2)

  # Add sample IDs to the sample_data
  phyloseq::sample_data(new_tree)$sampleid <- phyloseq::sample_names(new_tree)

  output_file_path = paste0(output_folder_rds_files, project_name, "_phyloseq_resolved_tree.rds")

  saveRDS(new_tree, file = output_file_path)

  log_message(paste("Phyloseq object saved as .rds object in", output_file_path), log_file)

  return(new_tree)
}
