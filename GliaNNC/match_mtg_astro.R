source("../scripts/pfc.fun.R")


## Subset to Astrocytes
seu <- readRDS(file = paste0(dataDir, "Human_MTG_org_05072021.rds"))
mastro <- subset(seu, cluster %in% c("Astro L1-2 FGFR3 GFAP", "Astro L1-6 FGFR3 SLC14A1"))
mhvg <- FindVariableFeatures(mastro, nfeatures = 500) %>%
            VariableFeatures()


hastro <- readRDS(file = paste0("../recluster_v2/load_files/Cls_bysp_detail_Astro_Human.rds"))
DefaultAssay(hastro) <- "RNA"
dfchvg <- SplitObject(hastro, split.by = "samplename") %>%
            lapply(., function(x) FindVariableFeatures(x, nfeatures = 1200)) %>%
            SelectIntegrationFeatures(., nfeatures = 1000)


Idents(mastro) <- "cluster"
mavgs <- log(AverageExpression(mastro)$RNA + 1)


Idents(hastro) <- "fig1cluster"
havgs <- log(AverageExpression(hastro)$RNA + 1)



sh_genes <- intersect(rownames(mavgs), rownames(havgs))
hvg <- union(mhvg, dfchvg) %>%
            intersect(., sh_genes)


mat <- cor(mavgs[hvg, ], havgs[hvg, ], method = "p")


mord <- c("Astro L1-2 FGFR3 GFAP", "Astro L1-6 FGFR3 SLC14A1")
hord <- c("iAstro GFAP FABP7", "pAstro AQP4 SLC1A2", "fAstro GFAP AQP1", "rAstro AQP4 OSMR")
pdf(paste0(outputdir, "MTG-dlPFC.human.match.pdf"), width = 7, height = 5)
pheatmap::pheatmap(mat[mord, hord], cluster_rows =FALSE, cluster_cols = FALSE, color = colorRampPalette(viridis::viridis(3))(30), border_color = NA, show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 5, fontsize_row = 8)
dev.off()




















