source("../scripts/pfc.fun.R")
library(tibble)


spdf2arry <- function(df, all_pairs) {
	## Add empty cols for missing pairs (diagonal Human-Human, Chimp-Chimp, etc)
	mis_pairs <- setdiff(all_pairs, colnames(df))
	if (length(mis_pairs) > 0){
		mis_df <- matrix(0, nrow = nrow(df), ncol = length(mis_pairs), dimnames = list(rownames(df), mis_pairs))
		df <- cbind(as.matrix(df), mis_df)
	}
	df <- df[, all_pairs] %>%
			as.matrix()

	dexinfo <- lapply(1:nrow(df), function(i) {matrix(df[i, ], nrow = 4, ncol = 4, byrow = TRUE)}) %>%
			do.call(c, .) %>%
			array(., dim = c(4, 4, nrow(df)), dimnames = list(all_sps, all_sps, rownames(df)))
	return(dexinfo)
}


SummariseArray <- function(ary, all_sps) {
	genes <- dimnames(ary)[[3]]
	rSUM <- lapply(genes, function(x) rowSums(ary[, , x])) %>%
			setNames(., genes) %>%
			as.data.frame(., check.names = FALSE) %>%
			as.matrix() %>% t()
	cSUM <- lapply(genes, function(x) colSums(ary[, , x])) %>%
			setNames(., genes) %>%
			as.data.frame(., check.names = FALSE) %>%
			as.matrix() %>% t()


	ept_mat <- matrix(0, nrow = length(genes), ncol = length(all_sps), dimnames = list(genes, all_sps))

	##------------------------------------------------
	mod3_idx <- which(apply(cSUM, 1, function(x) sum(x == 3) == 1))
	tem_mat <- cSUM[mod3_idx, ]
	tem_mat[tem_mat != 3] <- 1
	tem_mat[tem_mat == 3] <- 0
	ept_mat[mod3_idx, ] <- tem_mat


	mod1_idx <- which(apply(rSUM, 1, function(x) sum(x == 3) == 1))
	tem_mat <- rSUM[mod1_idx, ]
	tem_mat[tem_mat != 3] <- 0
	tem_mat[tem_mat == 3] <- 1
	ept_mat[mod1_idx, ] <- tem_mat


	mod2_idx <- which(apply(cSUM, 1, function(x) sum(x == 2) == 2) & 
						apply(rSUM, 1, function(x) sum(x == 2) == 2))
	tem_mat <- rSUM[mod2_idx, ]
	tem_mat[tem_mat != 2] <- 0
	tem_mat[tem_mat == 2] <- 1
	ept_mat[mod2_idx, ] <- tem_mat
	return(ept_mat)
}


pdata <- readRDS(file = paste0("Prop.difs.hres.custom_ref.FDR-02.rds"))


pair_ord <- c("Human_Chimpanzee", "Human_Rhesus", "Human_Marmoset",
			"Chimpanzee_Human", "Chimpanzee_Rhesus", "Chimpanzee_Marmoset", 
			"Rhesus_Human", "Rhesus_Chimpanzee", "Rhesus_Marmoset",
			"Marmoset_Human", "Marmoset_Chimpanzee", "Marmoset_Rhesus")
all_pairs <- c("Human_Human", pair_ord[1:4], "Chimpanzee_Chimpanzee", pair_ord[5:8],
				"Rhesus_Rhesus", pair_ord[9:12], "Marmoset_Marmoset")
dend <- readRDS(file = paste0("../MF1_overview/load_files/MF1_plot_dend_v3.rds"))
cls_ord <- rev(labels(dend))
all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")



## First perform a pair-wise test to get a p value data frame (row: cluster; col: species-pairs)
comp_df <- pdata %>%
		select(Cell.Type, Final.Parameter, pair) %>%
		reshape2::dcast(., Cell.Type ~ pair, value.var = "Final.Parameter") %>%
		column_to_rownames("Cell.Type")
comp_df[is.na(comp_df)] <- 0
comp_df[comp_df < 0] <- 0
comp_df[comp_df > 0] <- 1

comp_array <- spdf2arry(df = comp_df, all_pairs = all_pairs)
ept_mat <- SummariseArray(ary = comp_array, all_sps = all_sps)
isshare <- rowSums(ept_mat) == 0
pre_df <- ept_mat %>%
				as.data.frame(., check.names = FALSE) %>%
				rownames_to_column("cluster") %>%
				mutate(Shared = ifelse(isshare, 1, 0)) %>%
				mutate(radius = 0.5) %>%
				mutate(gene = "scCODA")
source("../MF7_contact/pie.fun.R")

p <- PlotScatterPie2(pie.data = pre_df, group.col = "cluster", feature.col = "gene", r.col = "radius", split.order = c("Human", "Chimpanzee", "Rhesus", "Marmoset", "Shared"), group.order = cls_ord, feature.order = "scCODA", scale.expression = FALSE, rsf = 1/0.5)
pdf(paste0(outputdir, "Prop.difs.hres.custom_ref.FDR-0.2.pie.summary.pdf"), width = 12, height = 4)
print(p);
dev.off() 







