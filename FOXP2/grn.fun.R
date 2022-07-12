GeneEnrich <- function(ratios, gene, ctp) {
	x <- which(rownames(ratios) == gene)
	if (length(x) == 0){
		stop("The input gene: ", gene, "is not contained in the ratios matrix")
	}

	expr <- ratios[x, ]
	mu <- mean(expr)
	thre <- mu * 2
	thre_floor <- 0.1
	thre_ceil <- quantile(expr, 0.9)
	if (thre < thre_floor) {
		thre <- thre_floor
	} 
	if (thre > thre_ceil){
		thre <- thre_ceil
	}

	## Highly expressed genes are define by expression ratio >= 0.2
	## lowly expressed genes: enriched i
	yy <- as.numeric(expr[ctp] >= thre | expr[ctp] >= 0.25)
	return(yy)
}


GeneDeplete <- function(ratios, gene, ctp) {
	x <- which(rownames(ratios) == gene)
	if (length(x) == 0){
		stop("The input gene: ", gene, "is not contained in the ratios matrix")
	}

	expr <- ratios[x, ]
	mu <- mean(expr)
	mu_bot <- mean(sort(expr, decreasing = FALSE)[1:6])
	thre <- (mu - mu_bot) * 0.6 + mu_bot


	## lowly expressed genes are define by expression ratio <= 0.1
	## highly expressed genes: using background expression
	yy <- as.numeric(expr[ctp] < thre | expr[ctp] < 0.1)
	return(yy)
}

EvaluateEdge <- function(sourceg, targetg, coef, ctp, ratios) {
	if (sourceg == "FOXP2" & coef > 0){
		yy <- GeneEnrich(ratios = ratios, gene = targetg, ctp = ctp)
	} else if (sourceg == "FOXP2" & coef < 0){
		yy <- GeneDeplete(ratios = ratios, gene = targetg, ctp = ctp)
	} else if (targetg == "FOXP2" & coef > 0){
		yy <- GeneEnrich(ratios = ratios, gene = sourceg, ctp = ctp)
	} else if (targetg == "FOXP2" & coef < 0){
		yy <- GeneDeplete(ratios = ratios, gene = sourceg, ctp = ctp)
	} else {
		yy <- 0
	}

	return(yy)
}


CalcSpecificity <- function(sourceg, targetg, ctp, logavgs) {
	gene <- c(sourceg, targetg) %>%
				setdiff(., "FOXP2")
	x <- which(rownames(logavgs) == gene)
	if (length(x) == 0){
		stop("The input gene: ", gene, "is not contained in the ratios matrix")
	}

	expr <- logavgs[x, ]
	mu <- mean(expr)
	spec <- (expr[ctp] + 0.1)/(mu + 0.1)
	return(spec)
}



PruneEdges <- function(df, ctp, ratios, logavgs) {
	df$qc <- sapply(1:nrow(df), function(x) EvaluateEdge(sourceg = df$source[x], targetg = df$target[x], coef = df$coef_mean[x], ctp = ctp, ratios = ratios))
	df$spec <- sapply(1:nrow(df), function(x) CalcSpecificity(sourceg = df$source[x], targetg = df$target[x], ctp = ctp, logavgs = logavgs))
	return(df)
}



library(igraph)
plot_tfnetwork.showctps <- function(linkdata, meta, file_name, color_use = c("lightgrey", "orange"), width_scale = 1, output_dir = new_outputdir, weight_thre = 1.30, plot.scale = 1, gset.path, group_cols = NULL) {
    label_genes <- meta$gene
    linkdata$coluse <- ifelse(linkdata$weight > 0, "red", "blue")
    linkdata$weight <- abs(linkdata$weight)

    set.seed(42)
    ## Build the net object
    net <- graph_from_data_frame(d=linkdata, vertices=meta, directed=TRUE)

    ##1. Colors
    V(net)$size <- case_when(
    	V(net)$gtype == "target" ~ 12,
    	V(net)$gtype == "TF" ~ 18, 
    	V(net)$gtype == "FOXP2" ~ 18)
    V(net)$label.color <- "black"
    V(net)$label.cex <- 0.5


    E(net)$color <- "#b0acac"
    E(net)$width <- 1


    ##--------------------------------------------------------------------------
    ## Use pie charts to visualize the vertices/nodes
    gset <- ifelse_check(is.character(gset.path), readRDS(gset.path), gset.path)
    gtable <- lapply(names(gset), function(x) data.frame(gene = gset[[x]], xx = 1, stringsAsFactors = FALSE) %>% setNames(., c("gene", x))) %>%
            purrr::reduce(full_join, by = "gene") %>%
            replace(is.na(.), 0) %>%
            column_to_rownames("gene")

    ## Add the FOXP2 part
    fp2 <- grep("FOXP2$", names(V(net)), value = TRUE)
    mat <- matrix(0, nrow = length(fp2), ncol = ncol(gtable), dimnames = list(fp2, colnames(gtable)))
    gtable <- as.data.frame(rbind(gtable, mat))
    gtable[, "Center"] <- ifelse(rownames(gtable) %in% fp2, 1, 0)

    
    shape.vec <- ifelse(rowSums(gtable[names(V(net)), ]) >= 2, "pie", "circle")
    shape.vec[V(net)$gtype == "FOXP2"] <- "csquare"
    new_gt <- gtable[names(V(net)), ] %>%
                    t() %>% 
                    as.data.frame(check.names = FALSE) %>% 
                    as.list()


    if (is.null(group_cols)){ 
        cls_cols <- gg_color_hue(ncol(gtable)) %>% 
                    setNames(., colnames(gtable))
    } else {
        cls_cols <- group_cols[colnames(gtable)]
    }
    V(net)$pie.color=list(c(cls_cols))



    ## Delete some unnecessary edges
    pdf(paste0(output_dir, file_name, ".pdf"), 1 * plot.scale * 6, 6 * plot.scale, useDingbats = FALSE)
    set.seed(42)
    l <- do.call("layout_with_fr", list(net)) 
    plot(net, edge.arrow.size=0.1, edge.curved=0.2,
    		arrow.width = 0.075, layout=l, main="layout", 
    		vertex.frame.color=NA,
            vertex.shape = shape.vec, 
            vertex.label = ifelse(names(V(net)) %in% label_genes, names(V(net)), ""), 
            vertex.pie = new_gt,
            vertex.label.dist = ifelse(V(net)$gtype == "target", 0.5, 0), 
            vertex.color = cls_cols[apply(gtable[names(V(net)), ], 1, function(x) which(x == 1)[1])],
            vertex.label.degree = pi/2) 
    legend("bottomleft", inset=.02, title="Pathways", colnames(gtable), fill=cls_cols, cex=0.3, ncol = 4)
    dev.off()
} 




