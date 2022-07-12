## Visualize species-exclusive markers on Heatmap
source("../scripts/pfc.fun.R")



resolution <- "fig1cluster"
all.markers <- readRDS(file = paste0(dataDir, "Markers.hprc.", resolution, ".rds"))
all.dex <- readRDS(file = paste0(dataDir, "DEX.hprc.", resolution, ".rds")) 



cls_order <- readRDS(file = paste0(dataDir, "PFC.cluster.order.rds"))
homocls <- readRDS(file = paste0(dataDir, "homolog.cluster.rds"))
cls_order <- cls_order[match(homocls, cls_order)]
cls_order <- cls_order[c(1:14, 16:20, 15, 21:length(cls_order))]

colss <- readRDS(file = paste0(dataDir, "cluster.color.rds"))
cls_col <- setNames(colss$color, rownames(colss))

ctp <- c("ExN", "InN", "NNC")[3]
homo.order <- switch(ctp, ExN = cls_order[1:36],
						InN = cls_order[37:81],
						NNC = cls_order[82:103])
all_species <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")



slim.dex <- all.dex %>%
					filter(pct.1 >= 0.3 & pct.2 <= 0.3 & ratio_fc >= 1.5 & avg_logFC >= log(2.5)) %>%
					group_by(cluster, speciespair) %>% 
					mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
					filter(p_val_adj <= 0.01)
saveRDS(slim.dex, file = paste0(dataDir, "HPRC_DEMarkers.", ctp, ".rds"))


dex.up <- lapply(all_species, function(sp) {
	sub.dex <- filter(slim.dex, species1 == sp) %>%
				group_by(cluster, gene) %>%
				summarize(nhits = n()) %>%
				arrange(desc(nhits)) %>%
				filter(nhits == 3) %>%
				mutate(color = cls_col[cluster])
	dex <- split(sub.dex, sub.dex$cluster) %>%
			.[homo.order] %>%
			do.call(rbind, .)
	return(dex)
	}) %>%
		do.call(rbind, .)
up_idx <- which(!duplicated(dex.up$gene))
up_genes <- dex.up$gene[up_idx]
up_cols <- dex.up$color[up_idx]



dex.down <- lapply(all_species, function(sp) {
	sub.dex <- filter(slim.dex, species2 == sp) %>%
				group_by(cluster, gene) %>%
				summarize(nhits = n()) %>%
				arrange(desc(nhits)) %>%
				filter(nhits == 3) %>%
				mutate(color = cls_col[cluster])
	dex <- split(sub.dex, sub.dex$cluster) %>%
			.[homo.order] %>%
			do.call(rbind, .)
	return(dex)
	}) %>%
		do.call(rbind, .)
down_idx <- which(!duplicated(dex.down$gene))
down_genes <- dex.down$gene[down_idx]
down_cols <- dex.down$color[down_idx]


sp.cls <- paste0(rep(all_species, each = length(homo.order)), "|",
			rep(homo.order, times = length(all_species)))
avg_data <- readRDS(file = paste0(dataDir, "Expr_avg_ratio_by_species.rds"))[["fig1cluster"]]$avg

avg_data <- avg_data[, sp.cls] %>% 
				as.matrix() %>% 
				t() %>% scale() %>% t() %>% 
				MinMax(., min = -2.5, max = 2.5)


up_data <- avg_data[up_genes, ]; rownames(up_data) <- paste0("up_", rownames(up_data))
down_data <- avg_data[down_genes, ]; rownames(down_data) <- paste0("down_", rownames(down_data))
data_use <- rbind(up_data, down_data)


gene_split <- setNames(c(rep("Up", length(up_genes)), rep("Down", length(down_genes))) %>% factor(., levels = c("Up", "Down")), rownames(data_use))


set.seed(4)
## Upregulated labels
lbs1 <- list(c("HEPN1", "VWCE", "HPSE2", "AQP1", "FREM2", "GALR1", "KLHL1", "PLPPR1", "FGF1", "PREX2", "DNAH17", "ACSBG1", "FOXP2", "CACNA1D", "DLEU7", "HERPUD1", "NR4A2", "C22orf34", "NUPR1", "ACKR3", "KCNMB1", "MYO16",
	"TNFRSF19", "CROCC", "DCX", "SEMA3C", "BFSP2", "SAMD5", "TNS3", "PLA2R1", "FGFR1", "JAM3", "SDC2", "EMB", "ALDH1A3", "SOX6", "DTX4", "DCHS2", "TTN"),
		c("MYO3A", "SNTG2", "NTNG1", "ARHGAP15", "SLC6A20", "EGFR", "SLC25A48", "ABCA10", "BTBD11", "CPM", "IL12RB2", "SHROOM3", "SLCO4A1", "DNAH11", "ABCC6",
			"MATN2", "MAP3K19", "LIG1", "PRIMA1", "PLD1", "FABP7", "DCDC1", "CDH23", "APBB1IP", "MLANA", "EYA2", "OCA2", "IL1RAP", "SDK2", "BOC", "PDE3A", "PCA3", "ADAM12","NHS", "WDR49", "CLEC10A", "SLC5A6", "GREM2", "F2R", "COLEC11", "FRMD3", "CCDC192", "CCBE1", "DSP", "GPC3")) %>%
			lapply(., function(x) intersect(x, up_genes)) %>%
			lapply(., function(x) sample(x, 20))

## Downregulated labels
lbs2 <- list(c("NLGN1", "GABRB2", "DNM1", "MDGA2", "SEMA3E", "GUCY1A3", "PLAT", "TMTC2", "MYO1B",
			"ALDH1L1", "SEMA5B", "FAM107B", "P2RX7", "ID3", "CEBPD", "CA4", "SMAD9", "S100A13","CLSTN2"),
			c("VAV3", "ETV6", "GAD1", "PTMA", "CFH", "FREM1", "FBLN5", "HS3ST5",
				"FAM184A", "AGT", "FOXP1", "ID4", "POU6F2", "SEMA3E", "COL4A5", "LAMP2", "TMTC4", "PPP2R3A", "ATM", "C3", "KCNQ3", "PTPRJ", "GALNT15", "GALNT18"))

label_genes <- lapply(1:2, function(x) setNames(c(lbs1[[x]], lbs2[[x]]), c(paste0("up_", lbs1[[x]]), paste0("down_", lbs2[[x]]))))

label_cols <- setNames(c(up_cols, down_cols), rownames(data_use)) %>%
				.[unlist(lapply(label_genes, names))]
source("../scripts/pfc.fun.R")
plot_dex_heatmap(data = data_use, split.by = "species", group.by = "cluster", label_genes = label_genes, label_col = label_cols, font_scale = 4, height_unit = 0.01, file_name = paste0("HPRC_DEX_", ctp), output_dir = outputdir, show_rownames = FALSE, color_limits = seq(-1, 2.5, 0.5), pdf_width = 18, pdf_height = 12, row_split = gene_split, heat_width = 7)






