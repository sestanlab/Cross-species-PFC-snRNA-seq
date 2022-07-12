source("../scripts/pfc.fun.R")
library(tidyr) ## gather
library(ggsankey) 


## Filter some cells
FilterTransfer <- function(meta, ref_col, query_col, thre = 0.02) {
	meta$pair <- paste0(meta[, ref_col], "|", meta[, query_col])
	kp_pairs <- meta %>% 
		group_by(!!sym(ref_col), !!sym(query_col)) %>%
		summarize(nhits = n()) %>%
		ungroup() %>%
		group_by(!!sym(ref_col)) %>%
		mutate(ratio = nhits/sum(nhits)) %>%
		filter(ratio >= thre) %>%
		mutate(pair = paste0(!!sym(ref_col), "|", !!sym(query_col))) %>%
		.$pair

 	meta <- meta %>%
		filter(pair %in% kp_pairs) %>%
		select(-pair)
  return(meta)
}




pfc <- readRDS(file = paste0(dataDir, "Final_PFC_HPRC.InN.rm1680.07102021.rds"))
layer_infor <- readRDS(file = paste0(inputdir, "PFC_transed_layers_by_HMTG.meta.rds"))
main_data <- pfc@meta.data[, c("fig1cluster", "mres", "species", "lres")]
main_data$layer <- layer_infor[rownames(main_data), "predicted.id"]



##---------------------------------------------------------------
## Sankey of SST clusters & laminar distributions
plot_data <- main_data %>%
				subset(species == "Human") %>%
				subset(mres %in% c("SST", "SST_NPY", "TH")) %>%
				subset(fig1cluster != "InN LHX6 TH STON2")

plot_data <- FilterTransfer(meta = plot_data, ref_col = "fig1cluster", query_col = "layer", thre = 0.025)

cls_ord <- c("InN SST HTR2C", "InN SST ADAMTSL1", "InN SST HS3ST5", "InN SST PLPP4", "InN SST TRPC7", "InN SST STK32A", "InN SST FREM1", "InN SST KLHL14", "InN SST THSD7B", "InN SST TH GABRQ", "InN SST NPY") %>% rev()
layer_ord <- paste0("L", 2:6) %>% rev()


##cls_anno <- readRDS(paste0(dataDir, "cluster.color.rds"))
cls_cols <- setNames(rep(NA, length(cls_ord)), cls_ord)
layer_cols <- c("#bf5b17", "#beaed4", "#386cb0", "#7fc97f", "#f0027f", "#fdc086") %>% setNames(., paste0("L", 1:6))



source("./inn.fun.R")
p <- plot_ggsankey(meta = plot_data, x_col = "layer", next_col = "fig1cluster", x_ord = layer_ord, next_ord = cls_ord, type = "sankey", space= 0.001, width= 0.1, x_cols = layer_cols[rev(layer_ord)], next_cols = cls_cols)

pdf(paste0(outputdir, "SST_cluster_layer_sankey_human.pdf"), width = 10, height = 10)
print(p)
dev.off()



##---------------------------------------------------------------
## Sankey of Layer 1 cells & cluster distributions
plot_data <- main_data %>%
				subset(species == "Human")
plot_data <- FilterTransfer(meta = plot_data, ref_col = "fig1cluster", query_col = "layer", thre = 0.025)


plot_data <- plot_data %>%
				subset(layer == "L1") %>%
				mutate(gp = gsub("^SST_NPY$", "SST", mres)) %>%
				mutate(gp = gsub("^TH$", "SST", mres)) %>%
				mutate(gp = gsub("^PVALB_CHC$", "PVALB", mres)) %>%
				mutate(gp = gsub("^PVALB_CHC$", "PVALB", mres))

cls_ord <- c("ADARB2","LAMP5_RELN", "VIP") %>% rev()


##cls_anno <- readRDS(paste0(dataDir, "cluster.color.rds"))
cls_cols <- setNames(rep(NA, length(cls_ord)), cls_ord)
layer_cols <- c("#bf5b17", "#beaed4", "#386cb0", "#7fc97f", "#f0027f", "#fdc086") %>% setNames(., paste0("L", 1:6))



source("./inn.fun.R")
p <- plot_ggsankey(meta = plot_data, x_col = "layer", next_col = "mres", x_ord = "L1", next_ord = cls_ord, type = "sankey", space= 0.001, width= 0.1, x_cols = "#bf5b17", next_cols = cls_cols)

pdf(paste0(outputdir, "SST_L1cells_cluster_sankey_human.pdf"), width = 10, height = 10)
print(p)
dev.off()
















