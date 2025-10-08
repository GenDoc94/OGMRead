library(circlize)
library(purrr)

case74 <- cb_filter |> select(Id, data) |> filter(Id == 74) |> unnest()

cytobands <- read.cytoband(species = "hg38")

deletions <- case74 |> 
        filter(Type == "deletion") |> 
        select(Id, RefcontigID1, RefStartPos, RefEndPos) |>
        rename(chr = RefcontigID1) |>
        mutate(chr = paste("chr", chr, sep = "")) |>
        ungroup() |> select(-Id)


circos.clear()

circos.par(start.degree = 90)  # Cromosoma 1 arriba

circos.initializeWithIdeogram(
        cytoband = cytobands,
        species = "hg38",
        plotType = c("ideogram", "axis", "labels"),
        sort.chr = paste0("chr", 1:22)
)

circos.trackPlotRegion(
        ylim = c(0, 1),
        track.height = 0.05,
        bg.border = "grey40",
        panel.fun = function(x, y) {
                chr <- get.cell.meta.data("sector.index")
                
                # Seleccionar deleciones del cromosoma actual
                d <- deletions[deletions$chr == chr, ]
                
                if (nrow(d) > 0) {
                        # Dibujar puntos en el centro de cada deleción
                        mid_pos <- (d$RefStartPos + d$RefEndPos) / 2
                        circos.points(mid_pos, rep(0.5, nrow(d)), 
                                      col = "darkred", 
                                      pch = 16,        # tipo de punto (16 = círculo sólido)
                                      cex = 0.8)       # tamaño del punto
                }
        }
)