## Generate dendrograms across all cell types
source("../scripts/pfc.fun.R")
library(ggdendro)
library(dendextend)




updata_name <- function(x) {
    xx <- gsub("L5 BCL11B SYT2 POU3F1", "L5 FEZF2 BCL11B POU3F1", x) %>%
            gsub("InN SST TH GABRQ", "InN SST HGF GABRQ", .) %>%
            gsub("InN LHX6 TH STON2", "InN LHX6 HGF STON2", .) %>%
            gsub("iAstro", "Astro", .) %>%
            gsub("fAstro", "Astro", .) %>%
            gsub("pAstro", "Astro", .) %>%
            gsub("rAstro", "Astro", .) %>%
            gsub("cOPC PDGFRA TOP2A", "OPC PDGFRA PCDH15", .) %>%
            gsub("nOligo", "Oligo", .) %>%
            gsub("eOligo", "Oligo", .) %>%
            gsub("mOligo", "Oligo", .) %>%
            gsub("lOligo", "Oligo", .) %>%
            gsub("aEndo", "Endo", .) %>%
            gsub("vEndo", "Endo", .) %>%
            gsub("cEndo", "Endo", .) %>%
            gsub("aaSMC", "SMC", .) %>%
            gsub("aSMC", "SMC", .) %>%
            gsub("vSMC", "SMC", .) %>%
            gsub("ABC", "VLMC",.) %>%
            gsub("VLMCA8", "ABCA8",.)
    return(xx)
}



## Set cluster order and colors
cls_cols <- readRDS(file = "../data/cluster.color.rds")
cls_cols <- cls_cols[setdiff(rownames(cls_cols), "cOPC PDGFRA TOP2A"), ] ## cOPC was merged with OPC
rownames(cls_cols) <- updata_name(rownames(cls_cols))
cls_cols <- setNames(cls_cols$color, rownames(cls_cols))



cls_ord <- readRDS(file = paste0(dataDir, "PFC.cluster.order.rds"))
cls_ord <- cls_ord[c(1:6, 23:26, 15:17, 12:14, 7:11, 18:21, 22, 27:34, 36:37, 35, 38:40, 
	41:49, 51, 50, 52, 54, 53, 55, 59:60, 58, 62, 56:57, 61, 63, 64, 65, 70, 68:69, 66, 72, 67, 73:74, 71, 75:79, 82, 80:81, 83:86,
					87:90, 92:115)] ## removed idx == 93, cOPC
cls_ord <- updata_name(cls_ord)



avgs <- readRDS(file = paste0("../data/Expr_avg_ratio_by_all.rds"))[["fig1cluster"]]$avg
avgs <- avgs[, !grepl("cOPC", colnames(avgs))]
colnames(avgs) <- updata_name(colnames(avgs))



all.markers <- readRDS(file = paste0(dataDir, "Markers.hprc.", "fig1cluster", ".rds"))


##------------------------------------------------------------------------
## ExN tree
exn_ord <- cls_ord[1:40]
exavgs <- avgs[, grep("^L[0-9]", colnames(avgs), value = TRUE)]
mres.markers <- readRDS(file = paste0(dataDir, "Markers.hprc.", "mres", ".rds"))
exn_mars  <- mres.markers %>%
                    subset(group == "ExN") %>%
                    subset(pct.1 >= 0.25 & ratio_fc >= 1.5 & avg_logFC >= 0.35 & p_val_adj <= 0.01) %>%
                    group_by(cluster, gene) %>% 
                    summarize(nhits = n(), mfc = mean(ratio_fc)) %>%
                    filter(nhits >= 2) %>%
                    ungroup() %>%
                    group_by(cluster) %>%
                    top_n(100, wt = mfc) %>%
                    .$gene %>% unique()
exsimi <- cor(exavgs[exn_mars, ], method = "p")
exn_dend <- hclust(as.dist(1-exsimi), method = "ward.D2") %>%
            as.dendrogram()
exn_dend <- reorder(exn_dend, match(colnames(exavgs), exn_ord), agglo.FUN=mean)



##------------------------------------------------------------------------
## InN tree
inn_ord <- cls_ord[41:86]
inavgs <- avgs[, grep("^InN", colnames(avgs), value = TRUE)]
inn_mars  <- all.markers %>%
                    subset(group == "InN") %>%
                    subset(pct.1 >= 0.25 & ratio_fc >= 1.25 & avg_logFC >= 0.35 & p_val_adj <= 0.01) %>%
                    group_by(cluster, gene) %>% 
                    summarize(nhits = n(), mfc = mean(ratio_fc)) %>%
                    filter(nhits >= 2) %>%
                    ungroup() %>%
                    group_by(cluster) %>%
                    top_n(200, wt = mfc) %>%
                    .$gene %>% unique()
insimi <- cor(inavgs[inn_mars, ], method = "p")
inn_dend <- hclust(as.dist(1-insimi), method = "ward.D2") %>%
            as.dendrogram()
inn_dend <- reorder(inn_dend, match(colnames(inavgs), inn_ord), agglo.FUN=mean)



##------------------------------------------------------------------------
## NNC tree
nnc_ord <- cls_ord[87:114]
nncavgs <- avgs[, grep("^L[0-9]|InN", colnames(avgs), value = TRUE, invert = TRUE)]
nnc_mars  <- all.markers %>%
                    subset(group == "NNC") %>%
                    subset(pct.1 >= 0.25 & ratio_fc >= 1.5 & avg_logFC >= 0.35 & p_val_adj <= 0.01) %>%
                    group_by(cluster, gene) %>% 
                    summarize(nhits = n(), mfc = mean(ratio_fc)) %>%
                    filter(nhits >= 2) %>%
                    ungroup() %>%
                    group_by(cluster) %>%
                    top_n(150, wt = mfc) %>%
                    .$gene %>% unique()
nnsimi <- cor(nncavgs[nnc_mars, ], method = "p")
nnc_dend <- hclust(as.dist(1-nnsimi), method = "ward.D2") %>%
            as.dendrogram()
nnc_dend <- reorder(nnc_dend, match(colnames(nncavgs), nnc_ord), agglo.FUN=mean)





##------------------------------------------------------------------------
## NNC tree
gunion <- unique(c(exn_mars, inn_mars, nnc_mars))


major_dist <- dist(rbind(rowMeans(exavgs[gunion, ]), rowMeans(inavgs[gunion, ])))
major_dist[1] <- 0.9
dend_p = as.dendrogram(hclust(major_dist))
dend_neuron = merge_dendrogram(dend_p, list(exn_dend, inn_dend), only_parent = FALSE)



major_dist2 <- dist(rbind((rowMeans(exavgs[gunion, ]) + rowMeans(inavgs[gunion, ]))/2, rowMeans(nncavgs[gunion, ])))
major_dist2[1] <- 0.9
dend_p2 = as.dendrogram(hclust(major_dist2))
dend_final = merge_dendrogram(dend_p2, list(dend_neuron, nnc_dend), only_parent = FALSE) %>%
				rev()


pdf(paste0(outputdir, "MF1._dendrogram.demo.pdf"), width = 10, height = 15)
par(mar = c(4, 4, 4, 4))
##grid.dendrogram(dend_m, test = TRUE)
dend_final %>% set("labels_col", value = cls_cols[labels(dend_final)[order.dendrogram(dend_final)]]) %>% # change color
            set("labels_cex", 0.5) %>% 
            ##rev() %>% 
            plot(main = "MF1", horiz = TRUE) # plot
dev.off()
saveRDS(dend_final, file = paste0("./load_files/MF1_plot_dend_v3.rds"))
















