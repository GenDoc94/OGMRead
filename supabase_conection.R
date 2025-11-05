library(httr2)
library(jsonlite)
library(tidyverse)

supabase_url <- "https://gavggbrvhzddfnfbkjun.supabase.co"
api_key <- Sys.getenv("SUPABASE_KEY")

#DDx
res <- request(paste0(supabase_url, "/rest/v1/DDx")) |>
        req_headers(apikey = api_key, Authorization = paste("Bearer", api_key)) |>
        req_perform()
ddx <- fromJSON(resp_body_string(res)) %>% mutate(label=Dx) %>% select(-Dx)

#DMuestra
res <- request(paste0(supabase_url, "/rest/v1/DMuestra")) |>
        req_headers(apikey = api_key, Authorization = paste("Bearer", api_key)) |>
        req_perform()
dmuestra <- fromJSON(resp_body_string(res))

#Muestras
res <- request(paste0(supabase_url, "/rest/v1/Muestras")) |>
        req_headers(apikey = api_key, Authorization = paste("Bearer", api_key)) |>
        req_perform()
metadata <- fromJSON(resp_body_string(res))  %>% 
        mutate(Id=NumBN, .before=1) %>% select(-NumBN) %>% #rename doesnt work
        left_join(ddx, by = c("Dx" = "Cod")) %>%
        mutate(Dx = label) %>%
        select(-label) %>% mutate(Dx = as.factor(Dx)) %>%
        left_join(dmuestra, by = c("Muestra" = "Cod")) %>%
        mutate(Muestra = TipoMuestra) %>%
        select(-TipoMuestra) %>% mutate(Muestra = as.factor(Muestra))

rm(api_key, ddx, dmuestra, res, supabase_url)