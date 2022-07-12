source("../scripts/pfc.fun.R")


## Subset to Astrocytes
seu <- readRDS(file = paste0("../MF4_FOXP2/motor_cortex/", "Motor.", "Human", ".full.rds"))
mastro <- subset(seu, cluster_label %in% c("Astro L1 FGFR3 SERPINI2", "Astro L1-6 FGFR3 AQP1", "Astro L1-6 FGFR3 PLCG1"))
mhvg <- FindVariableFeatures(mastro, nfeatures = 500) %>%
            VariableFeatures()


hastro <- readRDS(file = paste0("../recluster_v2/load_files/Cls_bysp_detail_Astro_Human.rds"))
DefaultAssay(hastro) <- "RNA"
dfchvg <- SplitObject(hastro, split.by = "samplename") %>%
            lapply(., function(x) FindVariableFeatures(x, nfeatures = 1200)) %>%
            SelectIntegrationFeatures(., nfeatures = 1000)


Idents(mastro) <- "cluster_label"
mavgs <- log(AverageExpression(mastro)$RNA + 1)


Idents(hastro) <- "fig1cluster"
havgs <- log(AverageExpression(hastro)$RNA + 1)



sh_genes <- intersect(rownames(mavgs), rownames(havgs))
hvg <- union(mhvg, dfchvg) %>%
            intersect(., sh_genes)


mat <- cor(mavgs[hvg, ], havgs[hvg, ], method = "p")


mord <- c("Astro L1 FGFR3 SERPINI2", "Astro L1-6 FGFR3 PLCG1", "Astro L1-6 FGFR3 AQP1")
hord <- c("iAstro GFAP FABP7", "pAstro AQP4 SLC1A2", "fAstro GFAP AQP1", "rAstro AQP4 OSMR")
pdf(paste0(outputdir, "Motor-dlPFC.human.match.pdf"), width = 7, height = 5)
pheatmap::pheatmap(mat[mord, hord], cluster_rows =FALSE, cluster_cols = FALSE, color = colorRampPalette(viridis::viridis(3))(30), border_color = NA, show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 5, fontsize_row = 8)
dev.off()




















