library(circlize)
library(purrr)

# insertion -> #2e868b V
# deletion -> #f3794b V
# inversión -> #6fa9db
# traslocation -> #eb008b
# duplication -> #9999ff V
# CNV gain -> #9999ff V
# CNV loss -> #f3794b V

case74 <- cb_filter |> select(Id, data) |> filter(Id == 74) |> unnest() #trans, del, gain
case61 <- cb_filter |> select(Id, data) |> filter(Id == 61) |> unnest() #ins, del
case83 <- cb_filter |> select(Id, data) |> filter(Id == 83) |> unnest() #trans, del, dup, gain, loss
case66 <- cb_filter |> select(Id, data) |> filter(Id == 66) |> unnest() #trans, inv, dup, gain, loss


case <- case66

cytobands <- read.cytoband(species = "hg38")

deletions <- case |> 
        filter(Type == "deletion") |> 
        select(Id, RefcontigID1, RefStartPos, RefEndPos) |>
        rename(chr = RefcontigID1) |>
        mutate(chr = paste("chr", chr, sep = "")) |>
        ungroup() |> select(-Id) |>
        mutate(chr = case_when(
                chr == "chr23" ~ "chrX",
                chr == "chr24" ~ "chrY",
                TRUE ~ chr))

insertions <- case |> 
        filter(Type == "insertion") |> 
        select(Id, RefcontigID1, RefStartPos, RefEndPos) |>
        rename(chr = RefcontigID1) |>
        mutate(chr = paste("chr", chr, sep = "")) |>
        ungroup() |> select(-Id) |>
        mutate(chr = case_when(
                chr == "chr23" ~ "chrX",
                chr == "chr24" ~ "chrY",
                TRUE ~ chr))

duplications <- case |> 
        filter(Type == "duplication") |> 
        select(Id, RefcontigID1, RefStartPos, RefEndPos) |>
        rename(chr = RefcontigID1) |>
        mutate(chr = paste("chr", chr, sep = "")) |>
        ungroup() |> select(-Id) |>
        mutate(chr = case_when(
                chr == "chr23" ~ "chrX",
                chr == "chr24" ~ "chrY",
                TRUE ~ chr))

gains <- case |> 
        filter(Type == "gain") |> 
        select(Id, RefcontigID1, RefStartPos, RefEndPos, fractionalCopyNumber) |>
        rename(chr = RefcontigID1) |>
        mutate(chr = paste("chr", chr, sep = "")) |>
        ungroup() |> select(-Id) |>
        mutate(chr = case_when(
                chr == "chr23" ~ "chrX",
                chr == "chr24" ~ "chrY",
                TRUE ~ chr))

losses <- case |> 
        filter(Type == "loss") |> 
        select(Id, RefcontigID1, RefStartPos, RefEndPos, fractionalCopyNumber) |>
        rename(chr = RefcontigID1) |>
        mutate(chr = paste("chr", chr, sep = "")) |>
        ungroup() |> select(-Id) |>
        mutate(chr = case_when(
                chr == "chr23" ~ "chrX",
                chr == "chr24" ~ "chrY",
                TRUE ~ chr))


circos.clear()

circos.par(start.degree = 90)  # Chr1 at 12am.

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
                                      col = "#f3794b", 
                                      pch = 16,        # tipo de punto (16 = círculo sólido)
                                      cex = 0.8)       # tamaño del punto
                }
                
                # Seleccionar deleciones del cromosoma actual
                i <- insertions[insertions$chr == chr, ]
                
                if (nrow(i) > 0) {
                        # Dibujar puntos en el centro de cada deleción
                        mid_pos <- (i$RefStartPos + i$RefEndPos) / 2
                        circos.points(mid_pos, rep(0.5, nrow(i)), 
                                      col = "#2e868b", 
                                      pch = 16,        # tipo de punto (16 = círculo sólido)
                                      cex = 0.8)       # tamaño del punto
                }
                
                # Seleccionar deleciones del cromosoma actual
                dp <- duplications[duplications$chr == chr, ]
                
                if (nrow(dp) > 0) {
                        # Dibujar puntos en el centro de cada deleción
                        mid_pos <- (dp$RefStartPos + dp$RefEndPos) / 2
                        circos.points(mid_pos, rep(0.5, nrow(dp)), 
                                      col = "#9999ff", 
                                      pch = 16,        # tipo de punto (16 = círculo sólido)
                                      cex = 0.8)       # tamaño del punto
                }
                
                
        }
)




circos.trackPlotRegion(
        ylim = c(min(losses$fractionalCopyNumber, gains$fractionalCopyNumber, 2) - 0.5,
                 max(losses$fractionalCopyNumber, gains$fractionalCopyNumber, 2) + 0.5),
        track.height = 0.1,
        bg.border = "grey40",
        panel.fun = function(x, y) {
                current_chr <- get.cell.meta.data("sector.index")
                if (!grepl("^chr", current_chr)) {
                        current_chr <- paste0("chr", current_chr)
                }
                
                xlim <- get.cell.meta.data("xlim")
                
                # Línea base punteada en 2
                circos.lines(x = xlim, y = rep(2, 2), col = "black", lty = 2, lwd = 1)
                
                # Filtrar ganancias y pérdidas del cromosoma actual
                d_gains <- gains %>% filter(chr == current_chr)
                d_losses <- losses %>% filter(chr == current_chr)
                
                # Función auxiliar para dibujar escalones
                draw_cnvs <- function(df, color) {
                        if (nrow(df) > 0) {
                                for (i in seq_len(nrow(df))) {
                                        xs <- c(df$RefStartPos[i], df$RefStartPos[i], df$RefEndPos[i], df$RefEndPos[i])
                                        ys <- c(2, df$fractionalCopyNumber[i], df$fractionalCopyNumber[i], 2)
                                        circos.lines(xs, ys, col = color, lwd = 2)
                                }
                        }
                }
                
                draw_cnvs(d_gains, "#9999ff")
                draw_cnvs(d_losses, "#f3794b")
        }
)
