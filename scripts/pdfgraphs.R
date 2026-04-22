# pdfgraphs.R

# Function to generate PDF plots
make_pdf_graphs <- function(filter_dx = NULL, filename = "graphs.pdf") {
        
        # Select IDs
        ids <- unique(complete_base$Id)
        
        if (!is.null(filter_dx)) {
                ids <- complete_base |> filter(Dx == filter_dx) |> pull(Id) |> unique()
        }
        
        pdf(filename, width = 10, height = 14)
        on.exit(dev.off())  # close PDF even if error occurs
        
        # Loop through IDs in chunks of 8
        for (i in seq(1, length(ids), by = 8)) {
                subset_ids <- ids[i:min(i + 7, length(ids))]
                
                # Set up 4x2 plotting grid
                par(mfrow = c(4, 2), mar = c(1, 1, 2, 1))
                
                # Plot each ID
                for (id in subset_ids) {
                        circleplot(id, cb_filter)
                        title(main = paste("Paciente", id), line = -1)
                }
        }
        
        message("PDF saved as ", filename)
}
