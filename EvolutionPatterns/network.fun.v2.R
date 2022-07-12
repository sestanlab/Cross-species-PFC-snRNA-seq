library(igraph)

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





PrepNet <- function(data, qt = 0.9) {
	simi <- as.matrix(cor(data, method = "pearson"))
	simi[upper.tri(simi, diag = TRUE)] <- NA
	df <- simi %>%
			reshape2::melt() %>%
			setNames(., c("from", "to", "weight")) %>%
			mutate(from = as.character(from), to = as.character(to)) %>%
			filter(!is.na(weight))

	df <- subset(df, weight >= quantile(df$weight, qt)) %>%
			mutate(weight = MinMax(weight, min = 0, max = 1))

	node_meta <- data.frame(cluster = rownames(simi), grou = "ExN", stringsAsFactors = FALSE)
	return(list(links = df, meta = node_meta))
}




CalcProps <- function(meta, ctp) {
	## Get the prop per sample in each major group
	gp_data <- meta %>%
				filter(group == ctp) %>%
				group_by(samplename, fig1cluster) %>%
                summarize(ncells = n(), species = unique(species)) %>%
                ungroup()


    if (ctp %in% c('ExN', 'InN')){
    	gp_data <- gp_data %>%
    			group_by(samplename) %>%
                mutate(percentage = ncells * 100/sum(ncells)) 
    } else if (ctp == 'NNC'){
    	gp_data$group <- "Vas"
    	gp_data$group[grep("Astro", gp_data$fig1cluster)] <- "Astro"
    	gp_data$group[grep("OPC|COP|Oligo", gp_data$fig1cluster)] <- "Oligo"
    	gp_data$group[grep("EBF1|Micro|SKAP1|Myeloid|Macro", gp_data$fig1cluster)] <- "Immune"
    	gp_data <- gp_data %>%
    			group_by(samplename, group) %>%
                mutate(percentage = ncells * 100/sum(ncells))
    }



	prop_mat <- gp_data %>%
				ungroup() %>%
                group_by(species, fig1cluster) %>%
                summarize(msize = mean(percentage)) %>%
                reshape2::dcast(., fig1cluster ~ species, value.var = "msize") %>%
                column_to_rownames("fig1cluster") %>%
                as.matrix()
    prop_mat[is.na(prop_mat)] <- 0
	all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")


	cls_use <- meta %>%
			filter(group == ctp) %>%
			.$fig1cluster %>% as.factor() %>% levels()
	prop_list <- lapply(cls_use, function(x) prop_mat[x, all_sps]) %>%
				setNames(., cls_use)

    return(prop_list)
}




## Plot the network
PlotNet <- function(res, file_name, dot_scale = 8){
	net <- graph_from_data_frame(d=res$links, vertices=res$meta, directed=FALSE)


	## Vertices parameters
	V(net)$size <- res$size[res$meta$cluster] * dot_scale
	sp_cols <- c(Human = "#FF420E", Chimpanzee = "#4CB5F5", Rhesus = "#89DA59", Marmoset = "#FFBB00")
	V(net)$pie.color <- list(sp_cols)


	## Edge parameters
	width_scale <- 1
	E(net)$width <- E(net)$weight * width_scale


	## Remove certain edges
	set.seed(0)
	lo <- do.call("layout_with_fr", list(net)) 


	## Determine the shape of the vertex (pie or circle)
	shape_vec <- sapply(res$props, function(x) ifelse(sum(x > 0) > 1, "pie", "circle")) %>% 
					setNames(., names(res$props))
	shape_col <- sapply(res$props, function(x) sp_cols[names(x)[which(x > 0)[1]]]) %>% 
					setNames(., names(res$props))


	pdf(paste0(outputdir, file_name, ".pdf"), width = 10, height = 10, useDingbats = FALSE)
	plot(net, edge.arrow.size=0.5, edge.curved=0, layout=lo, main="layout", vertex.label.cex = 0.5, vertex.frame.color= NA, asp = 1, add=FALSE, 
			vertex.shape = shape_vec, 
			vertex.pie = res$props,
			vertex.label = NA,
			vertex.color = shape_col,
			vertex.label.dist = 0, 
			edge.color = "#000000") 
	plot(net, edge.arrow.size=0.5, edge.curved=0, layout=lo, main="layout", vertex.label.cex = 0.5, vertex.frame.color= NA, asp = 1, add=FALSE, 
			vertex.shape = shape_vec, 
			vertex.pie = res$props,
			vertex.label = res$meta$cluster,
			vertex.color = shape_col,
			vertex.label.dist = 0, 
			edge.color = "#000000")
	dev.off()
}




