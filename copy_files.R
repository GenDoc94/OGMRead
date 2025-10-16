#Just to copy files from Dropbox to our directory
library(stringr)

# Define source folders and corresponding destination folders
folders <- list(
        list(source = "C:/Users/Juan/Dropbox/DATOS OGM/ARCHIVOS VARIANTES",
             dest   = "C:/Users/Juan/Desktop/R/OGMRead/variants_files",
             check_duplicates = TRUE),   # Check numeric prefix
        list(source = "C:/Users/Juan/Dropbox/DATOS OGM/ARCHIVOS JSON",
             dest   = "C:/Users/Juan/Desktop/R/OGMRead/json_files",
             check_duplicates = FALSE),  # Skip duplicate check
        list(source = "C:/Users/Juan/Dropbox/DATOS OGM/ARCHIVOS ANEUPLOIDÍAS",
             dest   = "C:/Users/Juan/Desktop/R/OGMRead/aneuploidy_files",
             check_duplicates = TRUE)    # Check numeric prefix
)

# Loop over each folder
for (f in folders) {
        
        # Create destination folder if it doesn't exist
        if (!dir.exists(f$dest)) dir.create(f$dest, recursive = TRUE)
        
        # List files in the source folder
        files_source <- list.files(f$source, full.names = TRUE)
        
        # Clean filenames to avoid Windows errors
        safe_files <- gsub("[:*?<>|]", "_", basename(files_source))
        
        
        # Full path of destination files
        files_dest <- file.path(f$dest, safe_files)
        
        # Copy files, overwriting existing ones
        file.copy(from = files_source, to = files_dest, overwrite = TRUE)
        
        # Check numeric prefixes for duplicates if required
        if (f$check_duplicates) {
                prefix_numbers <- str_extract(safe_files, "^\\d+")
                duplicated_numbers <- prefix_numbers[duplicated(prefix_numbers)]
                
                if (length(duplicated_numbers) == 0) {
                        message("OK (", f$dest, "): All files have a unique starting number.")
                } else {
                        message("WARNING (", f$dest, "): Duplicate numeric prefixes: ", 
                                paste(unique(duplicated_numbers), collapse = ", "))
                        
                        for (num in unique(duplicated_numbers)) {
                                dup_files <- safe_files[prefix_numbers == num]
                                message("Files with prefix ", num, ": ", paste(dup_files, collapse = "; "))
                        }
                }
        }
}

# Remove all variables created until now
rm(list = ls())