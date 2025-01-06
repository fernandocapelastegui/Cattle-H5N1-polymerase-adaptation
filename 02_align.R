# Load necessary libraries
library(tools)  # for file path manipulation


run_mafft <- function(input_fasta, output_fasta) {
  # Construct the MAFFT command. --auto lets MAFFT decide the best algorithm
  mafft_command <- paste("/usr/local/bin/mafft --auto --reorder", shQuote(input_fasta))
  
  # Run the command and capture the aligned output
  tryCatch({
    # Run the MAFFT command
    result <- system(mafft_command, intern = TRUE, wait = TRUE)
    
    # Write the output to the specified output file
    writeLines(result, output_fasta)
    
    # Check if the command ran successfully
    if (file.exists(output_fasta)) {
      cat(sprintf("MAFFT alignment completed successfully. Output saved to %s\n", output_fasta))
      now <- Sys.time()
      cat(format(now, "%H:%M:%S"), "\n")  # Print current time
    } else {
      cat(sprintf("Error running MAFFT. Output file not created: %s\n", output_fasta))
    }
  }, error = function(e) {
    cat(sprintf("Error while running MAFFT: %s\n", e$message))
  })
}


raw_folder_path <- paste0("/Users/capelastegui.f/git/bovine_tree_figure/Data/01_raw_data/",tdate,"/")
output_folder_path <- paste0("/Users/capelastegui.f/git/bovine_tree_figure/Data/02_aligned_data/", tdate)

# Make sure the output directory exists
if (!dir.exists(output_folder_path)) {
  dir.create(output_folder_path, recursive = TRUE)
}

# List all FASTA files in the raw folder and process them
file_list <- list.files(raw_folder_path, pattern = "_concat.fasta", full.names = TRUE)

for (input_fasta in file_list) {
  # Remove the original extension and add '_aligned.fasta'
  base_name <- file_path_sans_ext(basename(input_fasta))  # Get filename without extension
  output_fasta <- file.path(output_folder_path, paste0(base_name, "_aligned.fasta"))
  
  # Run alignment function
  run_mafft(input_fasta, output_fasta)
}
