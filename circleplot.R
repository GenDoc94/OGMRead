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
# 5. Limpiar gráfico previo
# -------------------------------
circos.clear()

# -------------------------------
# 6. Configurar posición inicial
# -------------------------------
circos.par(start.degree = 90)  # Cromosoma 1 arriba

# -------------------------------
# 7. Inicializar ideograma
# -------------------------------
circos.initializeWithIdeogram(
        cytoband = cytobands,
        species = "hg38",
        plotType = c("ideogram", "axis", "labels"),
        sort.chr = paste0("chr", 1:22)
)

# -------------------------------
# 8. Añadir deleciones e inversiones en el mismo track
# -------------------------------
circos.trackPlotRegion(
        ylim = c(0, 1),
        panel.fun = function(x, y) {
                chr <- get.cell.meta.data("sector.index")
                
                # Deleciones
                d <- deletions[deletions$chr == chr, ]
                if(nrow(d) > 0){
                        circos.rect(d$start, 0, d$end, 1, col = "red", border = NA)
                }
                
                # Inversiones
                inv <- inversions[inversions$chr == chr, ]
                if(nrow(inv) > 0){
                        circos.points(x = inv$pos, y = rep(0.5, nrow(inv)), col = "green", pch = 16, cex = 1)
                }
        }
)

# -------------------------------
# 9. Añadir translocaciones
# -------------------------------
walk(1:nrow(translocations), function(i) {
        circos.link(
                sector.index1 = translocations$chr1[i],
                point1 = c(translocations$start1[i], translocations$end1[i]),
                sector.index2 = translocations$chr2[i],
                point2 = c(translocations$start2[i], translocations$end2[i]),
                col = "blue", border = NA
        )
})

# -------------------------------
# 10. Añadir leyenda
# -------------------------------
legend(
        "topright",
        legend = c("Deleciones", "Translocaciones", "Inversiones"),
        fill = c("red", "blue", "green"),
        border = NA,
        bty = "n"
)



#library(karyoploteR)
#
# Suponiendo genome hg38
#kp <- plotKaryotype(genome = "hg38")
#
# Añadir deleciones
#kpRect(kp, chr = delec$chr, x0 = delec$start, x1 = delec$end, y0 = 0, y1 = 0.1, col = "red")