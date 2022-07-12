source("../scripts/pfc.fun.R")
source("./grn.fun.R")


##--------------------------------------------------------------------------
## Extract the FOXP2-related links that are significant
## in each cell type
load(file = paste0("./load_files/", "Expr_avg_ratio_multiome_4sps.rds"))
## logavgs, ratios, 


allctps <- c("L3_5_IT_1", "L3_5_IT_2", "L3_5_IT_3", "L5_6_NP", "L6_CT", "L6b", "Micro")
for (ctp in allctps){
	message(paste0("Working on cell type: ", ctp))
	links <- read.csv(file = paste0(inputdir, "mres=", ctp, ".raw_GRN.csv"), stringsAsFactors = FALSE) %>%
				filter(source == "FOXP2" | target == "FOXP2") %>%
				filter(p <= 0.01) %>%
				mutate(cluster = ctp)
	flinks <- PruneEdges(links, ctp = ctp, ratios = ratios, logavgs = logavgs) %>%
				filter(qc == 1)
	saveRDS(flinks, file = paste0(inputdir, "Expr_fil_GRNs_", ctp, ".rds"))
}


##--------------------------------------------------------------------------
## Merge the links from multiple cell types
## First merge the L3-5 IT links
links1 <- lapply(c("L3_5_IT_1", "L3_5_IT_2", "L3_5_IT_3"), function(ctp)
	readRDS(file = paste0(inputdir, "Expr_fil_GRNs_", ctp, ".rds"))
	) %>%
	do.call(rbind, .) %>%
	mutate(pair = paste0(source, "|", target)) %>%
	group_by(pair) %>%
	filter(coef_abs == max(coef_abs)) %>%
	as.data.frame()

## Then, other cell types
links2 <- lapply(c("L5_6_NP", "L6_CT", "L6b", "Micro"), function(ctp)
	readRDS(file = paste0(inputdir, "Expr_fil_GRNs_", ctp, ".rds"))
	) %>%
	do.call(rbind, .) %>%
	mutate(pair = paste0(source, "|", target)) %>%
	as.data.frame()

links <- rbind(links1, links2) %>%
		mutate(cluster = gsub("L3_5_IT_[0-9]", "L3_5_IT", cluster))
links$pattern <- case_when(
	(links$source == "FOXP2" & links$coef_mean > 0) ~ "3_out_pos", 
	(links$source == "FOXP2" & links$coef_mean < 0) ~ "4_out_neg", 
	(links$target == "FOXP2" & links$coef_mean > 0) ~ "1_in_pos", 
	(links$target == "FOXP2" & links$coef_mean < 0) ~ "2_in_neg")
saveRDS(links, file = paste0(inputdir, "Expr_fil_GRNs_combined.rds"))



##--------------------------------------------------------------------------
## Annotate human-specific microglia DEGs and mouse IUE genes
links <- readRDS(file = paste0(inputdir, "Expr_fil_GRNs_combined.rds"))


## Intersect with human-specific genes/Mouse IUE results
degs_micro <- readRDS(file = "~/myshare/PFC/FOXP2/Human_microglia_DEGs.rds") %>%
			lapply(., function(x) setdiff(x, "FOXP2")) %>%
			setNames(., c("Micro.up", "Micro.down"))
degs_iue <- readRDS(paste0("~/myshare/PFC/FOXP2/Mouse_DE_v2.rds")) %>%
				do.call(c, .) %>%
				lapply(., function(x) setdiff(x, "FOXP2"))
degs <- c(degs_micro, degs_iue)
alldegs <- unlist(degs, recursive = TRUE) %>% unique()
degtb <- lapply(degs, function(x) setNames(as.numeric(alldegs %in% x), alldegs)) %>%
				as.data.frame(., check.names = FALSE) %>%
				as.matrix()

for (ii in colnames(degtb)){
	ggs <- rownames(degtb)[degtb[, ii] == 1]
	links[, ii] <- sapply(1:nrow(links), function(x) {
		ctp <- links$cluster[x]
		coef <- links$coef_mean[x]
		gene <- c(links$source[x], links$target[x]) %>% setdiff(., c('FOXP2'))
		type <- strsplit(ii, ".", fixed = TRUE)[[1]][2]
		model_match <- (ii %in% c("Micro.up", "Micro.down") & ctp %in% "Micro") | (ii %in% c("L2-3.up", "L2-3.down", "L4.up", "L4.down", "L5-6.up", "L5-6.down") & ctp %in% c("L3_5_IT", "L5_6_NP", "L6_CT", "L6b"))

		if (model_match){
			if (type == "up"){
				yy <- ifelse(coef > 0 & gene %in% ggs, 1, 0)
			} else {
				yy <- ifelse(coef < 0 & gene %in% ggs, 1, 0)
			}
		} else {
			yy <- 0
		}
		return(yy)
		})
}

## Remove pattern 1_in_neg
links <- links %>%
			filter(pattern != "2_in_neg")
saveRDS(links, file = paste0(inputdir, "Expr_fil_GRNs_combined_annot.rds"))


















