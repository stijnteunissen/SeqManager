#' Clean Taxonomy Table in Phyloseq Object
#'
#' This function cleans and filters the taxonomy table within a `phyloseq` object.
#' It removes unclassified or ambiguous names, and replaces these and missing taxon names at genus level with placeholders derived from higher
#' taxonomic ranks. Optionally  this function filters specific taxa (e.g., chloroplasts, mitochondria, and ASVs unclassified at Kingdom and Phylum levels).
#' The cleaned taxonomy enables taxonomic agglomeration (`tax_glom`) to genus level without merging unclassified taxa from diverse ancestry.
#'
#' @param physeq A `phyloseq` object containing the microbiome data.
#'               This is the input object that the function processes.
#'
#' @param tax_filter A logical value (default `TRUE`) indicating whether to apply additional
#'                   filtering to remove unwanted taxa (e.g., chloroplasts, mitochondria).
#'                   \itemize{
#'                   \item `TRUE`: Apply filtering to remove taxa such as chloroplasts, mitochondria, and unclassified taxa.
#'                   \item `FALSE`: Skip filtering; the taxonomy table will only be cleaned but no taxa will be removed.
#'                   }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#' \item Cleans the taxonomy table:
#'   \itemize{
#'   \item Replaces ambiguous or placeholder taxa names (e.g., "uncultured organism") with `NA`.
#'   \item Assigns missing taxon names based on their higher-level ranks (e.g., "Phylum of Kingdom").
#'   \item Cleans the taxonomic ranks to remove overly nested placeholders.
#'   \item Ensures consistent naming for `Kingdom` (e.g., replacing "d__Bacteria" with "Bacteria").
#'   }
#' \item If `tax_filter = TRUE`, filters out unwanted taxa such as:
#'   \itemize{
#'   \item Chloroplasts (at `Class` or `Order` level).
#'   \item Mitochondria (at `Family` level).
#'   \item Unassigned taxa or taxa classified as `Eukaryota`.
#'   }
#' \item Saves the cleaned `phyloseq` object as an RDS file.
#' \item Logs the number of ASVs (amplicon sequence variants) removed during cleaning.
#' }
#'
#' The cleaned `phyloseq` object is returned for further analysis.
#'
#' @return A cleaned `phyloseq` object with ambiguous and unwanted taxa removed or replaced.
#'
#' @examples
#' \dontrun{
#' # Clean taxonomy table and apply filtering
#' physeq_cleaned <- tax_clean(physeq = physeq, tax_filter = TRUE)
#'
#' # Clean taxonomy table without filtering
#' physeq_cleaned_no_filter <- tax_clean(physeq = physeq, tax_filter = FALSE)
#' }
#'
#' @export
tax_clean = function(physeq = physeq,
                     tax_filter = TRUE) {

  psdata = physeq
  project_name = projects

  project_folder = paste0(base_path, projects)
  output_folder_csv_files = paste0(project_folder, "/output_data/csv_files/")
  output_folder_rds_files = paste0(project_folder, "/output_data/rds_files/Before_cleaning_rds_files/")


  # save uncleaned psdata
  output_file_path = paste0(output_folder_rds_files, project_name, "_phyloseq_uncleaned.rds")
  saveRDS(psdata, file = output_file_path)
  log_message(paste("Uncleaned phyloseq object saved as .rds object in", output_file_path), log_file)

  # ASV count start
  psdata_in = psdata

  # specify NA taxon name tags to last known taxon names
  tax.clean <- data.frame(phyloseq::tax_table(psdata))

  tax.clean2 =
    tax.clean %>%
    mutate_if(is.factor, as.character) %>%
    mutate(across(everything(),
                  ~ str_replace_all(.x,  "Ambiguous_taxa|metagenome|uncultured archeaon|uncultured bacterium|uncultured prokaryote|uncultured soil bacterium|uncultured rumen bacterium|uncultured compost bacterium|uncultured organism|uncultured$",
                                    replacement = NA_character_
                  ))) %>%
    replace(is.na(.), NA_character_) %>%
    mutate(Phylum = if_else(is.na(Phylum), paste0("Phylum of ", Kingdom), Phylum),
           Class = if_else(is.na(Class), paste0("Class of ", Phylum), Class),
           Order = if_else(is.na(Order), paste0("Order of ", Class), Order),
           Family = if_else(is.na(Family), paste0("Family of ", Order), Family),
           Genus = if_else(is.na(Genus), paste0("Genus of ", Family), Genus),
           Species = if_else(is.na(Species), paste0("Species of ", Genus), Species)
    ) %>%
    mutate(across(.cols = Kingdom:Species, .fns = ~if_else(str_detect(.,'\\bof\\b.*\\bof\\b'), paste0(word(., 1)," ", word(., 2)," ", word(., -1)), .))) %>%
    mutate(Kingdom = if_else(Kingdom == "d__Bacteria", "Bacteria", Kingdom),
           Kingdom = if_else(Kingdom == "d__Archaea", "Archaea", Kingdom))

  # put cleaned tax_table into phyloseq object
  phyloseq::tax_table(psdata) <- phyloseq::tax_table(as.matrix(tax.clean2))

  # apply taxa filter if tax_filter is TRUE
  # remove unwanted taxa such as Mitochondria, Chloroplasts, Unclassified Kingdom, Eukaryota, etc.
  if (tax_filter == TRUE) {
    psdata <-
      psdata %>% subset_taxa(
        Class != "Chloroplast" &
          Order != "Chloroplast" &
          Family != "Mitochondria" &
          Kingdom != "d__Eukaryota" &
          Kingdom != "Unassigned" &
          Phylum != "Phylum of d__Bacteria" &
          Phylum != "Phylum of d__Archaea" &
          Phylum != "Phylum of Bacteria")
  }

  # count ASVs after cleaning
  psdata_out = psdata

  # difference
  removed_ASV_count = ntaxa(psdata_in) - ntaxa(psdata_out)

  #save psdata after cleaning as RDS object
  output_file_path = paste0(output_folder_rds_files, project_name, "_phyloseq_cleaned.rds")
  saveRDS(psdata, file = output_file_path)
  log_message(paste("taxonomy table cleaned: ",  removed_ASV_count, " ASVs removed.  Phyloseq object saved as .rds object in", output_file_path), log_file)

  return(psdata)
}
