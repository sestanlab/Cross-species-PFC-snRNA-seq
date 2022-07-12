source("../scripts/pfc.fun.R")
source("./grn.fun.R")




##----------------------------------------------------------------------------------------------
## Venn diagram illustrate the overlap of targets
library(ggVennDiagram)
links <- readRDS(file = paste0(inputdir, "Expr_fil_GRNs_combined_annot.rds"))
all_pts <- levels(as.factor(links$pattern))
plist <- lapply(all_pts, function(pt) {
	data <- links %>%
			filter(pattern == pt)
	if (pt %in% c("3_out_pos", "4_out_neg")){ 
		xx <- split(data$target, data$cluster)
	} else if (pt %in% c("1_in_pos", "2_in_neg")){
		xx <- split(data$source, data$cluster)
	}

	yy <- list(L5_6 = xx[c("L6b", "L6_CT", "L5_6_NP")] %>% unlist() %>% unique(),
				Micro = xx[["Micro"]],
				L3_5_IT = xx[["L3_5_IT"]])

	p <- ggVennDiagram(yy, label_alpha = 0) +
			scale_fill_gradientn(colors = c("white", "red"))
	return(p)
	})


pdf(paste0(outputdir, "Targets_overlap.pdf"), width = 10, height = 10)
plot <- patchwork::wrap_plots(plist, nrow = 2, ncol = 2)
print(plot)
dev.off()



##----------------------------------------------------------------------------------------------
## Networks showing the targets from all subtypes
links <- readRDS(file = paste0(inputdir, "Expr_fil_GRNs_combined_annot.rds")) %>%
			group_by(cluster) %>%
			mutate(weight = as.numeric(cut(coef_abs, 10))) %>%
			ungroup() %>%
			mutate(weight = sign(coef_mean) * weight)
links$group <- links$cluster %>%
			gsub("L6b|L6_CT|L5_6_NP", "L5_6", .)


all_pts <- levels(as.factor(links$pattern))
flinks <- lapply(all_pts, function(pt) {
	if (pt %in% c("3_out_pos", "1_in_pos")){
		sublink <- links %>%
				filter(pattern == pt) %>%
				group_by(cluster) %>%
				slice_max(order_by = spec, n = 25) %>%
				ungroup() %>%
				group_by(source, target, group) %>%
				summarize(weight = mean(weight)) %>%
				ungroup() %>%
				mutate(index = paste0(group, "|", source, "|", target))
	} else {
		sublink <- links %>%
				filter(pattern == pt) %>%
				group_by(cluster) %>%
				slice_min(order_by = spec, n = 25) %>%
				ungroup() %>%
				group_by(source, target, group) %>%
				summarize(weight = mean(weight)) %>%
				ungroup() %>%
				mutate(index = paste0(group, "|", source, "|", target))
	}

	interest_links <- list(`3_out_pos` = c("L3_5_IT|FOXP2|CNR1", "L3_5_IT|FOXP2|PHLDB2", "L5_6|FOXP2|PHLDB2", "L5_6|FOXP2|CALCRL", "L5_6|FOXP2|FGF1", "L5_6|FOXP2|CALD1", "L5_6|FOXP2|PCSK5"),
				`4_out_neg` = c("L5_6|FOXP2|SHISA6"))
	if (pt %in% names(interest_links)){
		mis_idx <- setdiff(interest_links[[pt]], sublink$index)
		mis_links <- links %>%
				filter(pattern == pt) %>%
				group_by(source, target, group) %>%
				summarize(weight = mean(weight)) %>%
				ungroup() %>%
				mutate(index = paste0(group, "|", source, "|", target)) %>%
				filter(index %in% mis_idx)
		sublink <- rbind(sublink, mis_links)
	}


	sublink
	}) %>%
	setNames(., all_pts)
sapply(flinks, function(x) length(unique(c(x$target, x$source))))


flinks2 <- lapply(all_pts, function(pt){
	sub1 <- flinks[[pt]] %>%
			mutate(index2 = paste0(source, "|", target)) %>%
			.$index2 %>% unique
	sublink <- links %>%
				filter(pattern == pt) %>%
				group_by(source, target, group) %>%
				summarize(weight = mean(weight)) %>%
				ungroup() %>%
				mutate(index2 = paste0(source, "|", target)) %>%
				filter(index2 %in% sub1)
	
	sublink
	}) %>%
	setNames(., all_pts) 


plink <- flinks2[[1]][, c("source", "target", "weight", "group")]
plink$source[plink$source == "FOXP2"] <- paste0(plink$group, "|", plink$source)[plink$source == "FOXP2"]
plink$target[plink$target == "FOXP2"] <- paste0(plink$group, "|", plink$target)[plink$target == "FOXP2"]

gset_list <- lapply(all_pts, function(pt) {
	sublink <- links %>%
				filter(pattern == pt)
	glist <- split(sublink, sublink$group)%>%
				lapply(., function(x) union(x$source, x$target))
	glist
	}) %>%
	setNames(., all_pts)

pmeta <- data.frame(gene = union(plink$source, plink$target),
				stringsAsFactors = FALSE) %>%
			mutate(gtype = ifelse(gene %in% links$source, "TF", "target"))
pmeta$gtype[grep("\\|FOXP2", pmeta$gene)] <- "FOXP2"



source("./grn.fun.R")
plot_tfnetwork.showctps(linkdata = plink, meta = pmeta, file_name = "TF_network_overlaps_pattern-1", width_scale = 0.3, output_dir = outputdir, plot.scale = 1, gset.path = gset_list[[1]], group_cols = NULL) 










