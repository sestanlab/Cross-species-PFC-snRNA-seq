SinglePWTest <- function(response, group, p.adjust.method) {
	p.val <- pairwise.wilcox.test(response, group, p.adjust.method = p.adjust.method)$p.val
	df <- p.val %>%
			as.matrix() %>%
			reshape2::melt() %>%
			setNames(., c("species1", "species2", "pval")) %>%
			mutate(species1 = as.character(species1), species2 = as.character(species2)) %>%
			mutate(pair = paste0(species1, "_", species2)) %>%
			filter(!is.na(pval))
	df
}


PairWiseTest <- function(data, group.by, split.by, value.col, p.adjust.method = "fdr", fc.thre = 0.15, pseudo.count = 1){
	data <- as.data.frame(data)
	data.list <- split(data, data[, group.by])
	pval_data <- lapply(names(data.list), function(cls) {
		dd <- as.data.frame(data.list[[cls]])
		df <- SinglePWTest(response = dd[, value.col], group = dd[, split.by], p.adjust.method = p.adjust.method) %>%
				mutate(cluster = cls)
		df
		}) %>%
			do.call(rbind, .)

	## Duplicate the data with the reverse comparison direction 
	pval_data2 <- data.frame(species1 = pval_data$species2,
					species2 = pval_data$species1,
					pval = pval_data$pval,
					pair = paste0(pval_data$species2, "_", pval_data$species1),
					cluster = pval_data$cluster,
					stringsAsFactors = FALSE)
	pval_data <- rbind(pval_data, pval_data2)


	## Further correct the p values based on the fold changes [If fold changes does not exceed threshoold, p values will be set to 1]
	pval_data$fc <- sapply(1:nrow(pval_data), function(idx) {
		sp1 <- pval_data$species1[idx]
		sp2 <- pval_data$species2[idx]
		cls <- pval_data$cluster[idx]
		enp1 <- data[data[, split.by] == sp1 & data[, group.by] == cls, value.col]
		enp2 <- data[data[, split.by] == sp2 & data[, group.by] == cls, value.col]

		fc <- log2((mean(enp1) + pseudo.count)/(mean(enp2) + pseudo.count))
		return(fc)
		})
	pval_data$pval_final <- ifelse(pval_data$fc < fc.thre, 1, pval_data$pval)

	pval_mat <- reshape2::dcast(pval_data, as.formula(paste0("cluster", " ~ pair")), value.var = "pval_final") %>%
				column_to_rownames("cluster")
	return(list(df = pval_data, matrix = pval_mat))
}




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
	features <- dimnames(ary)[[3]]
	##------------------------------------------------
	## Transfer rSUM & cSUM to features * 4sp
	## 3 models: model 1:species exclusively enriched; mod 2: enriched in 2 species; mod3 depleted in one species

	## For each feature, evaluate the species-erichment pattern
	test_mod3 <- function(mat) {
		rsum <- rowSums(mat)
		csum <- colSums(mat)
		if (sum(csum == 3) == 1){
			pat_vec <- setNames(ifelse(csum == 3, 0, 1), colnames(mat))
			return(pat_vec)
		} else {
			return(setNames(rep(0, length(all_sps)), all_sps))
		}
	}


	test_mod1 <- function(mat) {
		rsum <- rowSums(mat)
		csum <- colSums(mat)
		if (sum(rsum == 3) == 1){
			pat_vec <- setNames(ifelse(rsum == 3, 1, 0), colnames(mat))
			return(pat_vec)
		} else {
			return(setNames(rep(0, length(all_sps)), all_sps))
		}
	}


	test_mod2 <- function(mat) {
		rsum <- rowSums(mat)
		ridx <- which(rsum == 2)
		bgidx <- setdiff(1:length(rsum), ridx)
		csum <- colSums(mat)
		if (sum(rsum == 2) == 2 & sum(csum[bgidx] >= 2) == 2){
			pat_vec <- setNames(ifelse(rsum == 2, 1, 0), colnames(mat))
			return(pat_vec)
		} else {
			return(setNames(rep(0, length(all_sps)), all_sps))
		}
	}

	all_df <- lapply(features, function(x) {
		mat <- ary[, , x]
		vec_list <- list(test_mod1(mat),
						test_mod2(mat),
						test_mod3(mat))
		pat_df <- vec_list %>%
					do.call(rbind, .) %>%
					as.data.frame() %>%
					mutate(model = c("m1", "m2", "m3"), feature = x)
		return(pat_df)
		}) %>%
			do.call(rbind, .)
	return(all_df)
}




