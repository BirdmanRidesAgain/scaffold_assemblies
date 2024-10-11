###################################################
### INFORMATION ###
###################################################
# get_contact_map.R
# KCollier 28 May 2024

# This is a command-line-usable script that outputs a synteny/contact map for two aligned genomes.
# It is generally assumed the output comes from RagTag.
# tidyverse is a dependency

###################################################
### HELPER FUNCTIONS
###################################################
print_help <- function() {
  sink(stdout(), type = "message")
  message("\nget_contact_map.R, Keiler Collier, 28 May 2023, v.0.2.0)
  
  Suggested usage: get_contact_map.R --paf <input_paf> --prefix <output_prefix>
          
  Options: 
          --raw: If set, all alignments are reported in the dotplot and quality filtering is not performed. Any micro/macro bounding remains.
                 It is not recommended to set this, as secondary, short and low-mapping-quality alignments can confound interpretation.
                 If --raw is unset (default behavior), secondary mappings, alignments under 500kbp and alignments with a mapq > 40 are removed.
         
          --target_lab <'Target'>: The label of the genome you are mapping the query sequences to. Defaults to 'Target'. 
                 More generally, this is the name of your 'truth set'. Reads in Target that don't have an alignment with Query will not be included in the dotplot.
          
          --query_lab <'Query'>: The label of the genome that you are attempting to find matches in. Defaults to 'Query'.
                 This is your genome-of-interest. All reads in Query will be displayed, regardless of matches with Target.
          
          --min_len <0>: Takes an integer value and sets the minimum alignment length for an alignment to be mapped. If unset, no micro/macro analyses will be done.
                Defaults to 0. Suggested boundary between avian macro/microchromosomes is 10 000 000bp. (10000000)
          
          --contig_plots: If set, 'n' number of contig-wise alignment plots will be written, where 'n' is the number of high-quality plots in 'Query'.
                IT IS NOT RECOMMENDED TO RUN THIS WITH --raw.
          --output_paf: This option outputs a csv of all filtered/plotted alignments. Turned on by default.
          ")
  sink(NULL, type = "message")
  quit()
}

check_installed_libs <- function(libs) {
  message("Checking for uninstalled libraries...")
  new_libs <- libs[!(libs %in% installed.packages()[,"Package"])]
  if(length(new_libs)) {
    message("Missing packages found. Installing.")
    install.packages(new_libs)
  }
  message("")
}

lib_load <- function(libs) {
  # Function to smooth out library loading process
  message("Loading libraries:")
  for (i in 1:length(libs)) { 
    suppressPackageStartupMessages(library(libs[i], character.only = TRUE))
    message(paste("\t",toString(libs[i]),"successfully loaded"))
  }
  message("All packages loaded\n")
}

print_params <- function(args) {
  message("User-defined parameters:")
  
  for (i in 1:length(args)) {
    name <- names(args[i])
    if (!name == "") {
      message(stringr::str_c("\t",name, ": ", toString(args[i])))
    }
  }
  message("")
}