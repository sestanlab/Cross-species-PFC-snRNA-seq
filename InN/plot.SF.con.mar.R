## This script will Plot the species-conserved markers across clusters.
source("../scripts/pfc.fun.R")
source("../scripts/exp.fun.R")


mar_res <- readRDS(file = paste0(dataDir, "Markers.hprc.fig7cluster.rds"))
cls_ord <- readRDS(file = paste0(dataDir, "PFC.cluster.order.rds"))[c(41:42, 44:86)]

mres_ord <- c("LAMP5_LHX6", "LAMP5_RELN", "VIP", "ADARB2", "SST_LPN", "SST_RELN", "SST_THSD7B", "SST_CALB1", "TH", "PVALB", "PVALB_CHC")


sub_mar <- mar_res %>%
			filter(group == "InN") %>%
			mutate(ratio_fc >= 1.25) %>%
			mutate(pct.1 >= 0.35) %>%
			group_by(species, cluster) %>%
			mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
			filter(p_val_adj <= 0.01) %>%
			ungroup()
new_mar <- sub_mar %>%
			group_by(cluster, gene) %>%
			summarize(nhits = n(), mavg = mean(ratio_fc)) %>%
			ungroup() %>%
			filter(nhits == 4) %>%
			group_by(cluster) %>%
			arrange(desc(mavg), .by_group = TRUE) %>%
			top_n(n = 5, wt = mavg)
markers <- split(new_mar$gene, new_mar$cluster) %>%
			.[mres_ord] %>%
			unlist() %>% 
			unique()

genes <- c("ADARB2", "PROX1", "NR2F2", "LAMP5", "CALB2", "VIP", "PAX6", "SP8", "RELN", "LHX6", "SOX6", "SST", "TH", "PVALB", "UNC5B") %>%
            c(., markers) %>% unique()




object <- readRDS(file = paste0(dataDir, "DotHeatPlot.pre.rds"))
CirclePlot.horizontal(avg = object$hres$avg, ratio = object$hres$ratio, features = genes, dot.min = 0.05, file_name = paste0("InN_conserved_mars", "_radius"), dot.scale = 4, scale.by = "radius", cluster.order = cls_ord, stroke.size = 0, stroke.matrix = NULL, width.scale = 1.2, height.base = 2.5, height.unit = 0.075, col.cex.sf = 0.8, row.cex.sf = 0.8)




















