# =========================================================
# Visualización genómica integrada con circlize (CNV, deleciones, inversiones, translocaciones)
# =========================================================

# -------------------------------
# 0. Cargar librerías
# -------------------------------
library(circlize)
library(purrr)

# -------------------------------
# 1. Descargar bandas cromosómicas
# -------------------------------
cytobands <- read.cytoband(species = "hg38")

# -------------------------------
# 2. Preparar deleciones
# -------------------------------
deletions <- data.frame(
        chr = c("chr1", "chr1", "chr2"),
        start = c(1e6, 5e7, 1e6),
        end   = c(2e6, 6e7, 2e6)
)

# -------------------------------
# 3. Preparar translocaciones
# -------------------------------
translocations <- data.frame(
        chr1 = c("chr1", "chr3"),
        start1 = c(1e7, 5e7),
        end1   = c(2e7, 6e7),
        chr2 = c("chr2", "chr5"),
        start2 = c(3e7, 7e7),
        end2   = c(4e7, 8e7)
)

# -------------------------------
# 4. Preparar inversiones
# -------------------------------
inversions <- data.frame(
        chr = c("chr1", "chr2"),
        pos = c(1.5e7, 7e7)
)

# -------------------------------
# 5. Preparar variaciones en el número de copias (CNV)
# -------------------------------
# Segmentos CN (ganancias/pérdidas)
cnv_segments <- data.frame(
        chr = c("chr1", "chr1", "chr2", "chr3"),
        start = c(1e7, 5e7, 1e7, 3e7),
        end   = c(2e7, 6e7, 2e7, 4e7),
        CN    = c(1, 3, 2, 4)
)

# Valores continuos (log2 ratio)
cnv_log2 <- data.frame(
        chr = rep("chr1", 10),
        start = seq(1e6, 10e6, by = 1e6),
        end = seq(1e6, 10e6, by = 1e6) + 1e5,
        log2ratio = c(-1, -0.5, 0, 0.3, 1, 0.8, 0.2, -0.2, -0.6, -1)
)

# -------------------------------
# 6. Limpiar gráfico previo
# -------------------------------
circos.clear()

# -------------------------------
# 7. Configurar posición inicial
# -------------------------------
circos.par(start.degree = 90)  # Cromosoma 1 arriba

# -------------------------------
# 8. Inicializar ideograma
# -------------------------------
circos.initializeWithIdeogram(
        cytoband = cytobands,
        species = "hg38",
        plotType = c("ideogram", "axis", "labels"),
        sort.chr = paste0("chr", 1:22)
)

# -------------------------------
# 9. Añadir track CN segmentado (ganancias/pérdidas)
# -------------------------------
circos.trackPlotRegion(
        ylim = c(0, 5),
        track.height = 0.08,
        bg.border = "grey40",
        bg.col = NA,
        panel.fun = function(x, y) {
                chr <- get.cell.meta.data("sector.index")
                d <- cnv_segments[cnv_segments$chr == chr, ]
                if (nrow(d) > 0) {
                        for (i in seq_len(nrow(d))) {
                                col <- if (d$CN[i] < 2) "#1E90FF66" else if (d$CN[i] > 2) "#FF000066" else "#B0B0B066"
                                circos.rect(d$start[i], 0, d$end[i], d$CN[i], col = col, border = "grey40")
                        }
                }
        }
)

# -------------------------------
# 10. Añadir track con log2 ratios
# -------------------------------
circos.genomicTrack(
        cnv_log2,
        ylim = c(-2, 2),
        track.height = 0.06,
        bg.border = "grey40",
        panel.fun = function(region, value, ...) {
                circos.genomicPoints(region, value, col = "darkorange", pch = 16, cex = 0.5)
                circos.genomicLines(region, value, col = "#FF8C0077", lwd = 1)
        }
)

# -------------------------------
# 11. Añadir deleciones e inversiones
# -------------------------------
circos.trackPlotRegion(
        ylim = c(0, 1),
        track.height = 0.05,
        bg.border = "grey40",
        panel.fun = function(x, y) {
                chr <- get.cell.meta.data("sector.index")
                
                # Deleciones
                d <- deletions[deletions$chr == chr, ]
                if (nrow(d) > 0) {
                        circos.rect(d$start, 0, d$end, 1, col = "#FF0000AA", border = "grey40")
                }
                
                # Inversiones
                inv <- inversions[inversions$chr == chr, ]
                if (nrow(inv) > 0) {
                        circos.points(x = inv$pos, y = rep(0.5, nrow(inv)), col = "#00FF00BB", pch = 16, cex = 1)
                }
        }
)

# -------------------------------
# 12. Añadir translocaciones
# -------------------------------
walk(1:nrow(translocations), function(i) {
        circos.link(
                sector.index1 = translocations$chr1[i],
                point1 = c(translocations$start1[i], translocations$end1[i]),
                sector.index2 = translocations$chr2[i],
                point2 = c(translocations$start2[i], translocations$end2[i]),
                col = "#4169E180",
                border = "grey40"
        )
})

# -------------------------------
# 13. Añadir leyenda
# -------------------------------
legend(
        "topright",
        legend = c(
                "CNV ganancia (CN>2)",
                "CNV pérdida (CN<2)",
                "CNV log2 ratio",
                "Deleciones",
                "Inversiones",
                "Translocaciones"
        ),
        fill = c("#FF000066", "#1E90FF66", "darkorange", "#FF0000AA", "#00FF00BB", "#4169E180"),
        border = "grey40",
        bty = "n"
)


#library(karyoploteR)
#
# Suponiendo genome hg38
#kp <- plotKaryotype(genome = "hg38")
#
# Añadir deleciones
#kpRect(kp, chr = delec$chr, x0 = delec$start, x1 = delec$end, y0 = 0, y1 = 0.1, col = "red")