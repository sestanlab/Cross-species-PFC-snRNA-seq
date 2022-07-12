source("../scripts/pfc.fun.R")


library(future)
plan("multiprocess", workers = 6)
options(future.globals.maxSize = 30*1000*1024^2)


load(file = paste0(inputdir, "Geirsdottir_human_micro.Rdata"))
## all_genes, seu

geir <- seu
rm(seu)
geir$cluster <- as.character(geir$seurat_clusters)
geir$inte.batch <- "Geirsdottir"
geir$stage <- "Geirsdottir"



hsb <- readRDS("../SpecclsMars/load_files/Immune_insp_Human_v3.rds")
hsb <- hsb[, !as.character(hsb@meta.data$fig1cluster) %in% c("B EBF1 IGKC", "T SKAP1 CD247", "Myeloid LSP1 LYZ")]
hsb@meta.data$inte.batch <- hsb@meta.data$samplename
hsb$stage <- "adult"




## Highly variable genes
hvg <- SplitObject(hsb, split.by = "samplename") %>%
				lapply(., function(x) FindVariableFeatures(x, nfeatures = 1500)) %>%
                SelectIntegrationFeatures(., nfeatures = 1000)



## Use the shared genes
sh_genes <- intersect(all_genes, rownames(hsb))
hvg <- intersect(hvg, sh_genes)
seu.list <- merge(x = hsb[sh_genes, ], y = geir[intersect(sh_genes, rownames(geir)), ]) %>%
            SplitObject(., split.by = "inte.batch")


file_name <- paste0("Inte_Geirsdottir_Micro_v2")


## Do the integration
source("./preprocess.fun.R")
seu <- Integratelist.seurat(obj.list = seu.list, hvg = hvg, file_name = file_name, input_dir = inputdir, inte.dims = 1:20, cluster.dims = 1:20, do.cluster = TRUE)





