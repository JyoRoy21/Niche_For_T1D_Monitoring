BiocManager::install("escape")
suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(SingleCellExperiment))
BiocManager::install("scran")
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))

#1.Loading Processed Single-Cell Data----

Tcells_w18<-readRDS("./annotated_TCells_T1D_Timepoints_v1.rds")
DimPlot(Tcells_w18)
sce.Tcells_w18 <- as.SingleCellExperiment(Tcells_w18, assay = "RNA")
DefaultAssay(Tcells_w18) = "RNA"

#2.Getting Gene Sets----
GS.hallmark <- getGeneSets(species="Mus musculus",library = "H")
#GS.hallmark <- msigdbr(species="Mus musculus", category="H") # Large df w/ categories
gene.sets <- list(Bcells = c("Abca1 "," Abcb8","Acaa2","Igh1","Igh2"),
                  Myeloid = c("Spi1","Fcer1g","Csf1r"),
                  Tcells = c("Cd3e", "Cd3d", "Cd3g", "Cd7","Cd8a"))

enrichment.scores <- escape.matrix(Tcells_w18, 
                                   gene.sets = GS.hallmark , 
                                   method = "ssGSEA",
                                   groups = 1000, 
                                   min.size = 5)

ggplot(data = as.data.frame(enrichment.scores), 
       mapping = aes(enrichment.scores[,1], enrichment.scores[,2])) + 
  geom_point() + 
  theme_classic() + 
  theme(axis.title = element_blank())
DefaultAssay(Tcells_w18) <- "RNA"
Layers(Tcells_w18)
Tcells_w18 <- runEscape(Tcells_w18, 
                        method = "ssGSEA",
                        gene.sets = GS.hallmark, 
                        groups = 1000, 
                        min.size = 5,
                        new.assay.name = "escape.ssGSEA",
                        use.layer = "counts"
)

sce.Tcells_w18 <- runEscape(sce.Tcells_w18, 
                      method = "UCell",
                      gene.sets = GS.hallmark, 
                      groups = 1000, 
                      min.size = 5,
                      new.assay.name = "escape.UCell")

heatmapEnrichment(sce.Tcells_w18, 
                  group.by = "ident",
                  assay = "escape.UCell",
                  scale = TRUE,
                  cluster.rows = TRUE,
                  cluster.columns = TRUE)

FeaturePlot(Tcells_w18, "Bcells") + 
  scale_color_gradientn(colors = colorblind_vector) + 
  theme(plot.title = element_blank())

Tcells_w18 <- performNormalization(sc.data = Tcells_w18, 
                                   assay = "escape.ssGSEA", 
                                   gene.sets = GS.hallmark)
heatmapEnrichment(Tcells_w18, 
                  assay = "escape.ssGSEA",
                  palette = "Spectral") 
