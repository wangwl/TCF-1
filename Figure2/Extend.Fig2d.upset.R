data=read.csv("data/boundary_enrich.csv")

ecdfplot <- function(mydata, x) {
  plot <- (ggplot(data = mydata, aes_string(x = x, colour = "color")) + 
             stat_ecdf() + scale_color_identity() + theme_classic())
}

pdf("FigureS2D.pdf",  useDingbats = FALSE, width=10, height=6)
upset(data, 
      queries = list(list(query = intersects, params = list("DP"), color='blue', active = T, query.name="A"),
                     list(query = intersects, params = list("DN3"), color='blue', active = T, query.name="A"),
                     list(query = intersects, params = list("DN4"), color='blue', active = T, query.name="A"),
                     list(query = intersects, params = list("DP", "DN4"), color='#63C5DA', active = T, query.name="B"),
                     list(query = intersects, params = list("DN3", "DN4"), color='#63C5DA', active = T, query.name="B"),
                     list(query = intersects, params = list("DP", "DN3"), color='#63C5DA', active = T, query.name="B"),
                     list(query = intersects, params = list("DP", "DN3", "DN4"), color='#63C5DA', active = T, query.name="B"),
                     list(query = intersects, params = list("ETP"), color='green', active = T, query.name="C"),
                     list(query = intersects, params = list("DN2"), color='green', active = T, query.name="C"),
                     list(query = intersects, params = list("CLP"), color='green', active = T, query.name="C"),
                     list(query = intersects, params = list("DN2", "ETP"), color='#C7EA46', active = T, query.name="D"),
                     list(query = intersects, params = list("CLP", "ETP"), color='#C7EA46', active = T, query.name="D"),
                     list(query = intersects, params = list("CLP", "DN2"), color='#C7EA46', active = T, query.name="D"),
                     list(query = intersects, params = list("CLP", "ETP", "DN2"), color = "yellow", active = T, query.name="E"),
                     list(query = intersects, params = list("CLP", "ETP", "DN2", "DN3"), color='orange', active = T, query.name="F"), 
                     list(query = intersects, params = list("CLP", "ETP", "DN2", "DN3", "DN4"), color = "red", active = T, query.name="G"), 
                     list(query = intersects, params = list("CLP", "ETP", "DN2", "DN3", "DN4", "DP"), color='purple', active = T, query.name="H")),
      attribute.plots = list(gridrows = 120, plots = list(list(plot = ecdfplot, x = "CLPPC1", queries = T),
                                                          list(plot = ecdfplot, x = "DPPC1", queries = T),
                                                          list(plot = ecdfplot, x = "CTCF.noTCF1", queries = T),
                                                          list(plot = ecdfplot, x = "CTCF.TCF1", queries = T),
                                                          list(plot = ecdfplot, x = "noCTCF.TCF1", queries = T),
                                                          list(plot = ecdfplot, x = "noCTCF.noTCF1", queries = T)), ncols = 3), 
      main.bar.color = "black", sets = c("CLP", "ETP", "DN2", "DN3", "DN4", "DP"), nintersects=100, group.by = "degree", keep.order = T, query.legend = 'top')

dev.off()
