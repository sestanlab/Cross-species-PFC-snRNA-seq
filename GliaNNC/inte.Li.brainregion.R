source("../scripts/pfc.fun.R")

library(future)
plan("multiprocess", workers = 6)
options(future.globals.maxSize = 30*1000*1024^2)


load(file = "./load_files/OMIX931-99-01.RData")


## Do the clustering analysis
fetal <- Microglia.Brain.Region
rm(Microglia.Brain.Region)


subfetal <- fetal[, fetal@meta.data$state != "Neuron"]
subfetal$cluster <- as.character(subfetal$seurat_clusters)
subfetal$inte.batch <- "fetal"
subfetal$stage <- "fetal"


hsb <- readRDS("../SpecclsMars/load_files/Immune_insp_Human_v3.rds")
hsb <- hsb[, !as.character(hsb@meta.data$fig1cluster) %in% c("B EBF1 IGKC", "T SKAP1 CD247", "Myeloid LSP1 LYZ")]
hsb@meta.data$inte.batch <- hsb@meta.data$samplename
hsb$stage <- "adult"


## Highly variable genes
hvg_fetal <- FindVariableFeatures(subfetal, nfeatures = 1500) %>%
                VariableFeatures()
hvg_hsb <- SplitObject(hsb, split.by = "samplename") %>%
				lapply(., function(x) FindVariableFeatures(x, nfeatures = 1500)) %>%
                SelectIntegrationFeatures(., nfeatures = 1000)
hvg <- union(hvg_hsb, hvg_fetal)


## Use the shared genes
sh_genes <- intersect(rownames(subfetal), rownames(hsb))
hvg <- intersect(hvg, sh_genes)
seu.list <- merge(x = hsb[sh_genes, ], y = subfetal[sh_genes, ]) %>%
            SplitObject(., split.by = "inte.batch")


file_name <- paste0("Inte_Li_fetalMicro")


## Do the integration
source("./preprocess.fun.R")
seu <- Integratelist.seurat(obj.list = seu.list, hvg = hvg, file_name = file_name, input_dir = inputdir, inte.dims = 1:20, cluster.dims = 1:20, do.cluster = FALSE)










