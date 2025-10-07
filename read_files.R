library(tidyverse)
library(readxl)

#List all files that start with number
archives <- list.files(
        path = "variants_files",
        pattern = "^[0-9]+_.*\\.txt$",  #start with number
        full.names = TRUE
)

#Read all and add column Id
all_variants <- archives %>%
        set_names() %>%  #keep the name of the archive
        map_dfr(~ {
                read_tsv(.x, na = c("null", "-", "NA"), comment = "") %>%
                        rename_with(~ sub("^#", "", .x), .cols = 1) %>%
                        mutate(Id = as.integer(str_extract(basename(.x), "^[0-9]+")), .before = 1)
        }) %>% mutate(
                across(c(Classification, Type, Zygosity), as.factor)
        )


ddx <- read_tsv("bionapp_files/DDx.txt",
                col_names = TRUE,
                locale = locale(encoding = "Latin1")) %>%
        mutate(label=Dx) %>% select(-Dx) #rename doesnt work

dmuestra <- read_tsv("bionapp_files/DMuestra.txt",
                     col_names = TRUE,
                     locale = locale(encoding = "Latin1"))

all_metadata <- read_tsv("bionapp_files/Muestras.txt",
                         col_names = TRUE,
                         locale = locale(encoding = "Latin1")) %>% 
        mutate(Id=NumBN, .before=1) %>% select(-NumBN) %>% #rename doesnt work
        left_join(ddx, by = c("Dx" = "Cod")) %>%
        mutate(Dx = label) %>%
        select(-label) %>% mutate(Dx = as.factor(Dx)) %>%
        left_join(dmuestra, by = c("Muestra" = "Cod")) %>%
        mutate(Muestra = TipoMuestra) %>%
        select(-TipoMuestra) %>% mutate(Muestra = as.factor(Muestra))

rm(archives, ddx, dmuestra)

#Nested base
n_variant <- all_variants |> group_by(Id) |> nest()

#Mixing data with metadata (just keeping data)
n_variant_complete <- left_join(n_variant,all_metadata, by="Id")

#Filtering P/LP/VUS variants
n_variant_filter <- n_variant %>%
        mutate(
                data = map(data, ~ filter(.x, Classification %in% c("Uncertain significance", 
                                                                    "Likely pathogenic", 
                                                                    "Pathogenic")))
        )