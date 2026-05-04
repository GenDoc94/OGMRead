copy_ogm_files <- function() {
        library(stringr)
        
        # 1. Definición de carpetas (Variables LOCALES a la función)
        folders <- list(
                list(source = "C:/Users/Juan/Dropbox/DATOS OGM/ARCHIVOS VARIANTES",
                     dest   = "C:/Users/Juan/Desktop/R/OGMRead/files/variants_files",
                     check_duplicates = TRUE),
                list(source = "C:/Users/Juan/Dropbox/DATOS OGM/ARCHIVOS JSON",
                     dest   = "C:/Users/Juan/Desktop/R/OGMRead/files/json_files",
                     check_duplicates = FALSE),
                list(source = "C:/Users/Juan/Dropbox/DATOS OGM/ARCHIVOS ANEUPLOIDÍAS",
                     dest   = "C:/Users/Juan/Desktop/R/OGMRead/files/aneuploidy_files",
                     check_duplicates = TRUE),
                list(source = "C:/Users/Juan/Dropbox/DATOS OGM/BASES CLÍNICAS",
                     dest   = "C:/Users/Juan/Desktop/R/OGMRead/files/clinical_files",
                     check_duplicates = FALSE)
        )
        
        # 2. Bucle de ejecución
        for (f in folders) {
                if (!dir.exists(f$dest)) dir.create(f$dest, recursive = TRUE)
                
                files_source <- list.files(f$source, full.names = TRUE)
                safe_files <- gsub("[:*?<>|]", "_", basename(files_source))
                files_dest <- file.path(f$dest, safe_files)
                
                file.copy(from = files_source, to = files_dest, overwrite = TRUE)
                
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
        # Al llegar aquí, R limpia automáticamente 'folders', 'f', 'safe_files', etc.
}

# 3. Llamada a la función para que se ejecute al hacer el source
copy_ogm_files()

# 4. Borramos solo la función para no dejar rastro, sin tocar el start_time
rm(copy_ogm_files)