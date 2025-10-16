ids <- unique(complete_base$Id)
#if you want to filter:
ids <- complete_base |> filter(Dx == "MM") |> pull(Id) |> unique()

pdf("MM_graphs.pdf", width = 10, height = 14)
on.exit(dev.off()) #just if there is an error...

# Recorremos los IDs de 8 en 8
for (i in seq(1, length(ids), by = 8)) {
        subset_ids <- ids[i:min(i + 7, length(ids))]
        
        # Definimos una cuadrícula de 4 filas x 2 columnas
        par(mfrow = c(4, 2), mar = c(1, 1, 2, 1))  # márgenes pequeños + espacio para título
        
        # Dibujamos los gráficos
        for (id in subset_ids) {
                circleplot(id, cb_filter)
                title(main = paste("Paciente", id), line = -1)  # título más cercano al gráfico
        }
}

dev.off()
