library(circlize)
library(purrr)
library(stringr)

# insertion -> #2e868b V
# deletion -> #f3794b V
# inversión -> #6fa9db
# traslocation -> #eb008b V
# duplication -> #9999ff V
# CNV gain -> #9999ff V
# CNV loss -> #f3794b V
# aneuploidy gain -> #9999ff (=duplication) V
# aneuploidy loss -> #f3794b (=deletion) V

case74 <- cb_filter |> select(Id, data) |> filter(Id == 74) |> unnest() #trans, del, gain
case61 <- cb_filter |> select(Id, data) |> filter(Id == 61) |> unnest() #ins, del
case83 <- cb_filter |> select(Id, data) |> filter(Id == 83) |> unnest() #trans, del, dup, gain, loss
case66 <- cb_filter |> select(Id, data) |> filter(Id == 66) |> unnest() #trans, inv, dup, gain, loss
case90 <- cb_filter |> select(Id, data) |> filter(Id == 90) |> unnest()

case <- case74

cytobands <- read.cytoband(species = "hg38")

aneuplos <- case |>
        filter(Type == "an-gain" | Type == "an-loss") |>
        select(Id, RefcontigID1, Type) |>
        rename(chr = RefcontigID1) |>
        mutate(chr = paste("chr", chr, sep = "")) |>
        ungroup() |> select(-Id) |>
        mutate(chr = case_when(
                chr == "chr23" ~ "chrX",
                chr == "chr24" ~ "chrY",
                TRUE ~ chr))

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
        filter(str_starts(Type, "duplication")) |> 
        select(Id, RefcontigID1, RefStartPos, RefEndPos) |>
        rename(chr = RefcontigID1) |>
        mutate(chr = paste("chr", chr, sep = "")) |>
        ungroup() |> select(-Id) |>
        mutate(chr = case_when(
                chr == "chr23" ~ "chrX",
                chr == "chr24" ~ "chrY",
                TRUE ~ chr))

inversions <- case |> 
        filter(str_starts(Type, "inversion")) |> 
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

translocations <- case |> 
        filter(str_starts(Type, "translocation")) |> 
        select(Id, RefcontigID1, RefcontigID2, RefStartPos, RefEndPos) |>
        rename(chr1 = RefcontigID1, chr2 = RefcontigID2) |>
        mutate(chr1 = paste("chr", chr1, sep = ""), chr2 = paste("chr", chr2, sep = "")) |>
        ungroup() |> select(-Id) |>
        mutate(chr1 = case_when(
                chr1 == "chr23" ~ "chrX",
                chr1 == "chr24" ~ "chrY",
                TRUE ~ chr1)) |>
        mutate(chr2 = case_when(
                chr2 == "chr23" ~ "chrX",
                chr2 == "chr24" ~ "chrY",
                TRUE ~ chr2))





circos.clear()

circos.par(start.degree = 90)  # Chr1 at 12am.

circos.initializeWithIdeogram(
        cytoband = cytobands,
        species = "hg38",
        plotType = c("ideogram", "axis", "labels"),
        sort.chr = paste0("chr", 1:22)
)

#SNVs
circos.trackPlotRegion(
        ylim = c(0, 1),
        track.height = 0.05,
        bg.border = "grey40",
        panel.fun = function(x, y) {
                chr <- get.cell.meta.data("sector.index")
                
                #DELECTIONS
                d <- deletions[deletions$chr == chr, ]
                
                if (nrow(d) > 0) {
                        # Dibujar puntos en el centro de cada deleción
                        mid_pos <- (d$RefStartPos + d$RefEndPos) / 2
                        circos.points(mid_pos, rep(0.5, nrow(d)), 
                                      col = "#f3794b", 
                                      pch = 16,        # tipo de punto (16 = círculo sólido)
                                      cex = 0.8)       # tamaño del punto
                }
                
                #INSERTIONS
                i <- insertions[insertions$chr == chr, ]
                
                if (nrow(i) > 0) {
                        # Dibujar puntos en el centro de cada deleción
                        mid_pos <- (i$RefStartPos + i$RefEndPos) / 2
                        circos.points(mid_pos, rep(0.5, nrow(i)), 
                                      col = "#2e868b", 
                                      pch = 16,        # tipo de punto (16 = círculo sólido)
                                      cex = 0.8)       # tamaño del punto
                }
                
                #DUPLICATIONS
                dp <- duplications[duplications$chr == chr, ]
                
                if (nrow(dp) > 0) {
                        # Dibujar puntos en el centro de cada deleción
                        mid_pos <- (dp$RefStartPos + dp$RefEndPos) / 2
                        circos.points(mid_pos, rep(0.5, nrow(dp)), 
                                      col = "#9999ff", 
                                      pch = 16,        # tipo de punto (16 = círculo sólido)
                                      cex = 0.8)       # tamaño del punto
                }
                
                #INVERSIONS
                inv <- inversions[inversions$chr == chr, ]
                
                if (nrow(inv) > 0) {
                        # Dibujar puntos en el centro de cada deleción
                        mid_pos <- (inv$RefStartPos + inv$RefEndPos) / 2
                        circos.points(mid_pos, rep(0.5, nrow(inv)), 
                                      col = "#6fa9db", 
                                      pch = 16,        # tipo de punto (16 = círculo sólido)
                                      cex = 0.8)       # tamaño del punto
                }
                
                
        }
)

#CNVs
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
                
                #ANEUPLOIDIES
                d_aneu <- aneuplos %>% filter(chr == current_chr)
                if (nrow(d_aneu) > 0) {
                        y_min <- min(c(losses$fractionalCopyNumber, gains$fractionalCopyNumber, 2)) - 0.4
                        for (i in seq_len(nrow(d_aneu))) {
                                color <- ifelse(d_aneu$Type[i] == "an-gain", "#9999ff", "#f3794b")
                                circos.lines(
                                        x = xlim, 
                                        y = rep(y_min, 2), 
                                        col = color, 
                                        lwd = 3
                                )
                        }
                }
        }
)

#TRANSLOCATIONS
for(i in seq_len(nrow(translocations))) {
        chr_a <- translocations$chr1[i]
        chr_b <- translocations$chr2[i]
        circos.link(
                sector.index1 = chr_a, point1 = translocations$RefStartPos[i],
                sector.index2 = chr_b, point2 = translocations$RefEndPos[i],
                col = "#eb008b", lwd = 1.5
        )
}

# Definir colores y etiquetas
legend_colors <- c(
        "#f3794b",   # deletion / CNV loss / aneuploidy loss
        "#2e868b",   # insertion
        "#6fa9db",   # inversion
        "#eb008b",   # translocation
        "#9999ff",   # duplication / CNV gain / aneuploidy gain
        "#9999ff",   # gain (CNV)
        "#f3794b"    # loss (CNV)
)

legend_labels <- c(
        "Deletion / Loss / Aneuploidy Loss",
        "Insertion",
        "Inversion",
        "Translocation",
        "Duplication / CNV Gain / Aneuploidy Gain",
        "Gain (CNV)",
        "Loss (CNV)"
)

# Añadir la leyenda a la derecha del gráfico
legend(
        "right",                 # posición
        legend = legend_labels,  # etiquetas
        col = legend_colors,     # colores
        pch = 16,                # símbolo sólido
        pt.cex = 1.2,            # tamaño del símbolo
        bty = "n",               # sin borde
        xpd = NA,                # permite dibujar fuera de la ventana
        cex = 0.8                # tamaño del texto
)