# install and load libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
pkgs <- c(
    "igraph",
    "RCy3"
)
BiocManager::install(pkgs, ask = FALSE, update = FALSE)
invisible(lapply(pkgs, library, character.only = TRUE))

# load GRN
grn <- readRDS("networks/estrogen_grn_igraph.rds")

# visualize in Cytoscape
cytoscapePing()
copyVisualStyle("default", "GENIE3")
setEdgeTargetArrowShapeDefault("Arrow", style.name = "GENIE3")

createNetworkFromIgraph(grn, title = paste("GRN of Estrogen-Associated Genes"), collection = "GRN Collection")
layoutNetwork("force-directed")
setVisualStyle("GENIE3")
