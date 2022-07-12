## Visualize species-exclusive markers on Heatmap
source("../scripts/pfc.fun.R")



##------------------------------------------------------------------
## data preparation/parameter setting

# load data
resolution <- "fig1cluster"
all.markers <- readRDS(file = paste0(dataDir, "Markers.hprc.", resolution, ".rds"))
all.dex <- readRDS(file = paste0(dataDir, "DEX.hprc.", resolution, ".rds")) 



# cluster order
cls_order <- readRDS(file = paste0(dataDir, "PFC.cluster.order.rds"))
homocls <- readRDS(file = paste0(dataDir, "homolog.cluster.rds"))
cls_order <- cls_order[match(homocls, cls_order)]
cls_order <- cls_order[c(1:14, 16:20, 15, 21:length(cls_order))]


ctp <- c("ExN", "InN", "NNC")[2]
homo.order <- switch(ctp, ExN = cls_order[1:36],
						InN = cls_order[37:81],
						NNC = cls_order[82:103])

all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")


###------------------------------------------------------------------------------------
## Get the species-exclusive markers [Both positive & negative]
demarkers <- readRDS(file = paste0(dataDir, "HPRC_DEMarkers.", ctp, ".rds"))

dex.up.base <- lapply(all_sps, function(sp) {
	sub.dem <- filter(demarkers, species1 == sp) %>%
				group_by(cluster, gene) %>%
				summarize(nhits = n(), mfc = min(ratio_fc)) %>%
				filter(nhits == 3) %>%
				ungroup() %>%
				group_by(cluster) %>%
				arrange(desc(mfc)) %>%
				ungroup() %>%
				mutate(species = sp)
	dem <- split(sub.dem, sub.dem$cluster) %>%
			.[homo.order] %>%
			do.call(rbind, .)
	dem
	}) %>% 
		do.call(rbind, .)
up.idx <- which(!duplicated(dex.up.base$gene))
up.genes <- dex.up.base$gene[up.idx]



dex.down.base <- lapply(all_sps, function(sp) {
	sub.dem <- filter(demarkers, species2 == sp) %>%
				group_by(cluster, gene) %>%
				summarize(nhits = n(), mfc = min(ratio_fc)) %>%
				filter(nhits == 3) %>%
				ungroup() %>%
				group_by(cluster) %>%
				arrange(desc(mfc)) %>%
				ungroup() %>%
				mutate(species = sp)
	dem <- split(sub.dem, sub.dem$cluster) %>%
			.[homo.order] %>%
			do.call(rbind, .)
	dem
	}) %>% 
		do.call(rbind, .)
down.idx <- which(!duplicated(dex.down.base$gene))
down.genes <- dex.down.base$gene[down.idx]





##------------------------------------------------------------------
## Plot
## Rename cluster names
updata_name <- function(x) {
    xx <- gsub("L5 BCL11B SYT2 POU3F1", "L5 FEZF2 BCL11B POU3F1", x) %>%
            gsub("InN SST TH GABRQ", "InN SST HGF GABRQ", .) %>%
            gsub("InN LHX6 TH STON2", "InN LHX6 HGF STON2", .) %>%
            gsub("iAstro", "Astro", .) %>%
            gsub("fAstro", "Astro", .) %>%
            gsub("pAstro", "Astro", .) %>%
            gsub("rAstro", "Astro", .) %>%
            gsub("cOPC", "OPC", .) %>%
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
            gsub("L6b", "L6B",.)
    return(xx)
}

homo.update <- updata_name(homo.order)
sp.cls <- paste0(rep(all_sps, each = length(homo.update)), "|", rep(homo.update, times = length(all_sps)))


avg_list <- lapply(all_sps, function(sp) {
	avg <- readRDS(file = paste0(inputdir, "Expr_avg_ratio_non-over_", sp, ".rds"))$avg
	colnames(avg) <- paste0(sp, "|", colnames(avg))

	allgenes <- extract_field(rownames(avg), "rm_start", "-") %>% setNames(., NULL)
	avg <- avg[remove_duplicates(allgenes, return_index = TRUE), ]
	rownames(avg) <- extract_field(rownames(avg), "rm_start", "-") %>% setNames(., NULL)
	return(avg)
	}) %>%
		setNames(., all_sps)

up.update <- lapply(avg_list, rownames) %>%
				Reduce("intersect", .) %>%
				intersect(up.genes, .)
down.update <- lapply(avg_list, rownames) %>%
				Reduce("intersect", .) %>%
				intersect(down.genes, .)
up_data <- lapply(avg_list, function(x) as.matrix(x[up.update, ,drop = FALSE])) %>%
				do.call(cbind, .) %>%
				.[, sp.cls] %>%
				as.matrix() %>% 
				t() %>% scale() %>% t() %>% 
				MinMax(., min = -2.5, max = 2.5)
rownames(up_data) <- paste0("up_", rownames(up_data))
down_data <- lapply(avg_list, function(x) as.matrix(x[down.update, ,drop = FALSE])) %>%
				do.call(cbind, .) %>%
				.[, sp.cls] %>%
				as.matrix() %>% 
				t() %>% scale() %>% t() %>% 
				MinMax(., min = -2.5, max = 2.5)
rownames(down_data) <- paste0("up_", rownames(down_data))

data_use <- rbind(up_data, down_data)
gene_split <- setNames(c(rep("Up", length(up.update)), rep("Down", length(down.update))) %>% factor(., levels = c("Up", "Down")), rownames(data_use))

source("./plot.dem.fun.R")
plot_dex_heatmap_nolabel(data = data_use, split.by = "species", group.by = "cluster", font_scale = 4, height_unit = 0.01, file_name = paste0("HPRC_DEX_non-over_", ctp), output_dir = outputdir, show_rownames = FALSE, color_limits = seq(-1, 2.5, 0.5), pdf_width = 18, pdf_height = 12, row_split = gene_split, heat_width = 6)




avg_over <- readRDS(file = paste0(dataDir, "Expr_avg_ratio_by_species.rds"))[["fig1cluster"]]$avg
colnames(avg_over) <- updata_name(colnames(avg_over))
avg_over <- avg_over[, sp.cls] %>% 
				as.matrix() %>% 
				t() %>% scale() %>% t() %>% 
				MinMax(., min = -2.5, max = 2.5)
up_over <- avg_over[up.update, ]; rownames(up_over) <- paste0("up_", rownames(up_over))
down_over <- avg_over[down.update, ]; rownames(down_over) <- paste0("down_", rownames(down_over))
over_use <- rbind(up_over, down_over)
plot_dex_heatmap_nolabel(data = over_use, split.by = "species", group.by = "cluster", font_scale = 4, height_unit = 0.01, file_name = paste0("HPRC_DEX_over_", ctp), output_dir = outputdir, show_rownames = FALSE, color_limits = seq(-1, 2.5, 0.5), pdf_width = 18, pdf_height = 12, row_split = gene_split, heat_width = 6)


## length
## up.genes: 2059
## down.genes: 668
## up.update: 1394
## down.update: 548



