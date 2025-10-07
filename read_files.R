library(tidyverse)

#List all files that start with number
archives <- list.files(
        path = "variants_files",
        pattern = "^[0-9]+_.*\\.txt$",  #start with number
        full.names = TRUE
)

#Read all and add column Id
all_data <- archives %>%
        set_names() %>%  #keep the name of the archive
        map_dfr(~ {
                read_tsv(.x, na = c("null", "-", "NA"), comment = "") %>%
                        rename_with(~ sub("^#", "", .x), .cols = 1) %>%
                        mutate(Id = as.integer(str_extract(basename(.x), "^[0-9]+")), .before = 1)
        }) %>% mutate(
                across(c(Classification, Type, Zygosity), as.factor)
        )

#Nested base
n_data <- all_data |> group_by(Id) |> nest()


#Filtering P/LP/VUS variants
n_f_data <- n_data %>%
        mutate(
                data = map(data, ~ filter(.x, Classification %in% c("Uncertain significance", 
                                                                    "Likely pathogenic", 
                                                                    "Pathogenic")))
        )