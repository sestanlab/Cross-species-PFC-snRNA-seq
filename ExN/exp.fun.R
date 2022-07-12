library(ComplexHeatmap)


###############################################################################################

## Dot Plot (modified version)

###############################################################################################
## Validate if the genes are contained in the dataset
ValidateFeatures <- function(features, all.genes) {
    ## Filter genes
    rm_genes <- setdiff(features, all.genes)

    ## Check features
    if (length(rm_genes) >= 1){
        warning(paste0("The following features are missing: ", paste(rm_genes, collapse = ", ")))
        features <- features[features %in% all.genes]
    }
    if (length(features) == 0){
        stop("No genes were found in the dataset")
    }
    return(features)
}



## Because some clusters are not conserved across species, need to add some NA cols to make them pseudo-conserved
FillNACol <- function(avg, ratio = NULL, subset.ident = NULL) {
    ## Metadata 
    plot_meta <- data.frame(spcls = colnames(avg), stringsAsFactors = FALSE) %>%
                        mutate(cluster = extract_field(spcls, 2, "|"), species = extract_field(spcls, 1, "|"))

    ## Subset the data if necessry
    if (!is.null(subset.ident)){
        plot_meta <- plot_meta %>%
                        filter(cluster %in% subset.ident) %>%
                        column_to_rownames("spcls")

        avg <- avg[, rownames(plot_meta)]
        if (!is.null(ratio)){
            ratio <- ratio[, rownames(plot_meta)]
        }

        dif_cls <- table(plot_meta$cluster) %>% .[. < 4] %>% names() %>% 
                    intersect(., subset.ident)
    } else {
        plot_meta <- plot_meta %>%
                        column_to_rownames("spcls")

        dif_cls <- table(plot_meta$cluster) %>% .[. < 4] %>% names()
    }



    ## Make sure we have NA columns for non-homologous clusters
    sp_use <- unique(plot_meta$species)


    if (length(dif_cls) > 0) {
        dif_cols <- paste0(rep(sp_use, times = length(dif_cls)), "|", 
                                    rep(dif_cls, each = 4)) %>%
                            setdiff(., colnames(avg))
        dif_mat <- matrix(NA, nrow = nrow(avg), ncol = length(dif_cols), dimnames = list(rownames(avg), dif_cols))
        avg <- cbind(avg, dif_mat)
        if (!is.null(ratio)){
            ratio <- cbind(ratio, dif_mat)
        }

        plot_meta <- data.frame(spcls = dif_cols, stringsAsFactors = FALSE) %>%
                        mutate(cluster = extract_field(spcls, 2, "|"), species = extract_field(spcls, 1, "|")) %>%
                        column_to_rownames("spcls") %>%
                        rbind(., plot_meta)
    }
    return(list(avg = avg, ratio = ratio, meta = plot_meta))
}


## Combinatory melt of avg and ratio matrices
CombnMelt <- function(avg, ratio, stroke = NULL) {
    Mat2Df <- function(mat, name) {
        df <- mat %>% 
                    as.matrix() %>%
                    reshape2::melt() %>%
                    setNames(., c("features.plot", "newcls", name)) %>%
                    mutate(features.plot = as.character(features.plot), newcls = as.character(newcls))
        return(df)
    }


    ## Combine average and expression ratio
    avg.exp <- Mat2Df(mat = avg, name = "avg.exp.scaled")
    pct_exp <- Mat2Df(mat = ratio, name = "pct.exp")
    pre.plot <- dplyr::full_join(avg.exp, pct_exp, by = c("features.plot", "newcls"))


    ## Add stroke data if present
    if (!is.null(stroke)){
        pre.plot <- Mat2Df(mat = ratio, name = "stroke") %>%
                        dplyr::full_join(pre.plot, ., by = c("features.plot", "newcls")) 
    }


    ## Further format the columns
    data.plot <- pre.plot  %>%
                    mutate(species = extract_field(newcls, 1, "|"),  cluster = extract_field(newcls, 2, "|"))%>%
                    select(-newcls)

    return(data.plot)
}



CirclePlot.horizontal <- function(avg, ratio, features, file_name, dot.min = 0, dot.scale = 4, scale.by = "radius", shape = 16, cluster.order = NULL, stroke.col = "black", stroke.size = 0.5, col.min = -2.5, col.max = 2.5, scale = TRUE, scale.min = NA, scale.max = NA, stroke.matrix = NULL, return.plot = FALSE, width.scale = 1, height.base = 1.5, font.scale = 1, height.unit = 0.3){
    split.order <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
    cols <- setNames(c("#FF420E","#4CB5F5","#89DA59","#FFBB00"), c("Human", "Chimpanzee", "Rhesus", "Marmoset"))
    ngradient <- 20

    ## Validate features
    features <- ValidateFeatures(features = features, all.genes = rownames(avg))
    

    ##-------------------------------------------------------------------------
    ## Scale data
    if (scale) {
        data.use <- avg[features, , drop = FALSE] %>%
                        t() %>% scale() %>% t() %>%
                        MinMax(data = ., min = col.min, max = col.max)
    } else {
        data.use <- log(x = avg[features, , drop = FALSE])
    }


    ## Fill NA cols (for non-conserved clusters)
    res <- FillNACol(avg = data.use, ratio = ratio[features, , drop = FALSE], subset.ident = cluster.order)


    ## Matrix to DF
    data.plot <- CombnMelt(avg = res$avg, ratio = res$ratio, stroke = stroke.matrix)


    data.plot$features.plot <- factor(x = data.plot$features.plot, levels = rev(x = features))
    data.plot$avg.exp.scaled <- as.numeric(x = cut(x = data.plot$avg.exp.scaled, breaks = ngradient))


    ## Set colors
    if (is.null(names(cols))){
        warning("The cols should be a named vector")
        cols <- setNames(cols, split.order)
    }
    data.plot$colors <- mapply(FUN = function(color, value) {
                    return(colorRampPalette(colors = c("grey", color))(ngradient)[value])
                    }, color = cols[data.plot$species], value = data.plot$avg.exp.scaled)


    ## Set pct.exp range
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    data.plot$pct.exp[data.plot$pct.exp < (dot.min*100)] <- NA



    ## Define the order the species & cluster
    if (!is.null(split.order)){
        data.plot$species <- factor(as.character(data.plot$species), levels = split.order)
    } 
    if (!is.null(cluster.order)){
        data.plot$cluster <- factor(as.character(data.plot$cluster), levels = cluster.order)
    }


    ##------------------------------------------------------------------
    ## Plot 
    ## Set dot scale functions
    scale.func <- switch(EXPR = scale.by, size = scale_size,
                        radius = scale_radius, 
                        stop("'scale.by' must be either 'size' or 'radius'"))

    if (shape == 21){ 
        if (!is.null(stroke.matrix)){
            data.plot$fills <- data.plot$colors
            data.plot$colors[data.plot$stroke != 0] <- "#000000"
            plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "cluster")) + 
            geom_point(mapping = aes_string(size = "pct.exp", fill = "fills",color = "colors", stroke = "stroke", shape = "shape")) + 
            scale_fill_identity()+
            scale_color_identity()+
            scale_shape_identity()
        } else {
            plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "cluster")) + 
            geom_point(mapping = aes_string(size = "pct.exp", fill = "colors"), color = stroke.col, stroke = stroke.size, shape = 21) + 
            scale_fill_identity()
        }
    } else {
        plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "cluster")) + 
            geom_point(mapping = aes_string(size = "pct.exp", color = "colors"), shape = shape) + 
            scale_color_identity()
    }
    plot <- plot +
            scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
            coord_flip() +
            theme_cowplot() +
            RotatedAxis() +
            facet_wrap(vars(species), nrow = 1, ncol = 4) +
            theme(panel.grid.minor = element_line(size = 0.1, color = "grey"),
                panel.grid.major = element_line(size = 0.1, color = "grey"),
                axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y=element_text(size = rel(0.75)), axis.text.x=element_text(size = rel(0.6 * font.scale)), strip.background = element_blank(), strip.text = element_blank(), panel.spacing = unit(0.1, "in"), legend.position = "bottom", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2)) + 
            guides(size = guide_legend(title = "Percent Expressed"))


    if (return.plot){
        return(plot)
    } else {
        pdf(paste0(outputdir, file_name, "_DotHori.pdf"), width = 15 * width.scale, height = ceiling(height.unit * length(features)) + height.base, useDingbats = FALSE)
        print(plot)
        dev.off()
    }
    
}
 

CirclePlot.vertical <- function(avg, ratio, features, file_name, dot.min = 0, dot.scale = 4, scale.by = "radius", shape = 16, cluster.order = NULL, stroke.col = "black", col.min = -2.5, col.max = 2.5, scale = TRUE, stroke.size = 0.5, stroke.matrix = NULL, scale.min = NA, scale.max = NA, return.plot = FALSE, width.scale = 1, height.base = 1.5, font.scale = font.scale, height.unit = 0.6){
    split.order <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
    cols <- c("#FF420E", "#4CB5F5", "#89DA59", "#FFBB00", "#AE017E") %>% setNames(., c("Human", "Chimpanzee", "Rhesus", "Marmoset", "Mouse"))
    ngradient <- 20

    ## Check features
    if (length(setdiff(features, rownames(avg))) >= 1){
        warning(paste0("The following features are missing: ", paste(setdiff(features, rownames(avg)), collapse = ", ")))
    }
    if (sum(features %in% rownames(avg)) == 0){
        stop("No genes are found in the avg matrix")
    }
    features <- features[features %in% rownames(avg)]


    ##-------------------------------------------------------------------------
    ## Set data
    if (scale) {
        data.use <- avg[features, , drop = FALSE] %>%
                        t() %>% scale() %>% t() %>%
                        MinMax(data = ., min = col.min, max = col.max)
    } else {
        data.use <- log(x = avg[features, , drop = FALSE])
    }


    ## Fill NA cols (for non-conserved clusters)
    res <- FillNACol(avg = data.use, ratio = ratio[features, , drop = FALSE], subset.ident = cluster.order)


    ## Matrix to DF
    data.plot <- CombnMelt(avg = res$avg, ratio = res$ratio, stroke = stroke.matrix)


    ## Set new feature column
    data.plot$features.plot <- factor(x = data.plot$features.plot, levels = features)
    data.plot$avg.exp.scaled <- as.numeric(x = cut(x = data.plot$avg.exp.scaled, breaks = ngradient))



    ## Set colors
    if (is.null(names(cols))){
        warning("The cols should be a named vector")
        cols <- setNames(cols, split.order)
    }
    data.plot$colors <- mapply(FUN = function(color, value) {
                    return(colorRampPalette(colors = c("lightgrey", color))(ngradient)[value])
                    }, color = cols[data.plot$species], value = data.plot$avg.exp.scaled)
    data.plot$colors[is.na(data.plot$colors)] <- "#FFFFFF"


    ## Set pct.exp range
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    data.plot$pct.exp[data.plot$pct.exp < (dot.min*100)] <- NA
    #data.plot$pct.exp[is.na(data.plot$pct.exp)] <- 0


    ## Define the order the species & cluster
    if (!is.null(split.order)){
        data.plot$species <- factor(as.character(data.plot$species), levels = rev(split.order))
    } 
    if (!is.null(cluster.order)){
        data.plot$cluster <- factor(as.character(data.plot$cluster), levels = cluster.order)
    }
    


    ##------------------------------------------------------------------
    ## Plot 
    ## Set dot scale functions
    scale.func <- switch(EXPR = scale.by, size = scale_size,
                        radius = scale_radius, 
                        stop("'scale.by' must be either 'size' or 'radius'"))

    if (shape == 21){ 
        if (!is.null(stroke.matrix)){
            data.plot$fills <- data.plot$colors
            data.plot$colors[data.plot$stroke != 0] <- "#000000"
            plot <- ggplot(data = data.plot, mapping = aes_string(x = "species", y = "cluster")) + 
            geom_point(mapping = aes_string(size = "pct.exp", fill = "fills",color = "colors", stroke = "stroke", shape = "shape")) + 
            scale_fill_identity()+
            scale_color_identity()+
            scale_shape_identity()
        } else {
            plot <- ggplot(data = data.plot, mapping = aes_string(x = "species", y = "cluster")) + 
            geom_point(mapping = aes_string(size = "pct.exp", fill = "colors"), color = stroke.col, stroke = stroke.size, shape = 21) + 
            scale_fill_identity()
        }
    } else {
        plot <- ggplot(data = data.plot, mapping = aes_string(x = "species", y = "cluster")) + 
            geom_point(mapping = aes_string(size = "pct.exp", color = "colors"), shape = shape) + 
            scale_color_identity()
    }
    plot <- plot +
            scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
            coord_flip() +
            theme_cowplot() +
            RotatedAxis() +
            facet_wrap(vars(features.plot), nrow = length(features), ncol = 1, strip.position = 'left') +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y=element_blank(), axis.text.x=element_text(size = rel(0.8 * font.scale)), strip.background = element_blank(), strip.text = element_text(size = rel(0.8), angle = 90), strip.text.y.left = element_text(size = rel(0.8), angle = 0), strip.placement = 'outside', panel.spacing = unit(0.05, "in"), legend.position = "bottom", axis.ticks.y = element_blank(), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2)) + 
            guides(size = guide_legend(title = "Percent Expressed"))

    if (return.plot){
        return(plot)
    } else {
        pdf(paste0(outputdir, file_name, "_DotVert.pdf"), width = 6 * width.scale, height = ceiling(height.unit * length(features)) + height.base)
        print(plot)
        dev.off()
    }
}




plot_sp_heatmap <- function(avg, features, scale.exp = TRUE, min_max = c(-2.5, 2.5), file_name, cluster.order = cls_order, font.scale = 1, height.scale = 1, width.scale = 1, return.plot = FALSE) {
    ## Set the parameters
    split.order <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
    sp.cols <- c("#FF420E","#4CB5F5","#89DA59","#89DA59","#FFBB00") %>% 
                setNames(., c("Human", "Chimpanzee", "Rhesus", "Macaque", "Marmoset"))


    ## Validate features
    features <- ValidateFeatures(features = features, all.genes = rownames(avg))
    data <- avg[features, ]


    
    ## Set the MinMax Values 
    if (scale.exp){
        data <- data %>% 
                as.matrix() %>%
                t() %>% scale() %>% t()
    }
    if (!is.null(min_max)){
        data <- data %>% MinMax(., min = min_max[1], max = min_max[2])
    }


    ## Fill NA cols (for non-conserved clusters)
    res <- FillNACol(avg = data, ratio = NULL)
    data <- res$avg
    plot_meta <- res$meta


    ## Reorganize the plot data
    plot_data <- lapply(split.order, function(x) {
            submeta <- plot_meta[plot_meta$species == x, ]
            td <- data[, rownames(submeta), drop = FALSE] %>% as.matrix()
            rownames(td) <- paste0(x, "|", rownames(td))
            colnames(td) <- as.character(submeta$cluster)
            td <- td[, cluster.order[cluster.order %in% colnames(td)], drop = FALSE]
            return(td)
            }) %>% 
                do.call(rbind, .)
    gene_orders <- lapply(rownames(data), function(x) grep(paste0("\\|", x, "$"), rownames(plot_data))) %>%
                    unlist(., use.names = FALSE)
    plot_data <- plot_data[gene_orders, ]



    ## Get the range of the heatmap legends
    legend_limits <- c(-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5)
    color_breaks <- legend_limits[legend_limits >= min(data[!is.na(data)]) & legend_limits <= max(data[!is.na(data)])]
    font_size <- font.scale * 8


    ##cls.cols <- gg_color_hue(length(group.order)) %>% setNames(., group.order)
    ##column_ha <- HeatmapAnnotation(cluster = factor(colnames(plot_data), group.order), col = list(cluster = cls.cols), annotation_height = unit(0.01, "in"))
    htlist <- Heatmap(plot_data, name = "scaled_expr", 
                    col = colorRamp2(color_breaks, viridis(length(color_breaks))), na_col = "white",
                    row_split = factor(extract_field(rownames(plot_data), 2, "|"), levels = rownames(data)), 
                    cluster_rows = FALSE, cluster_columns = FALSE, 
                    row_title_rot = 0, column_title = NULL, row_title_gp = gpar(fontsize = font_size), 
                    show_column_names = TRUE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 10),column_names_rot = 45, 
                    width = unit(5 * width.scale, "in"),
                    top_annotation = NULL, 
                    heatmap_legend_param = list(title = "Scaled expr", color_bar = "continuous", at = color_breaks))


    if (return.plot){
        return(htlist)
    } else {
        pdf_height <- ceiling(height.scale * 0.5 * length(features)) %>% MinMax(., min = 3, max = 30)
        pdf(paste0(outputdir, file_name, "_heatmap.pdf"), width = 5 * width.scale + 2, height = pdf_height)
        draw(htlist)
        dev.off()
    }
}


 

plot_sp_all <- function(object = object, resolution = c("mres", "hres")[1], features, file_name, ident.use = c("ExN", "InN", "NNC", "All", NULL)[4]) {
    if (is.null(ident.use)){
        ident.use <- "All"
    }


    ## Set object
    avgs <- object[[resolution]]$avg
    ratios <- object[[resolution]]$ratio
    dex <- object[[resolution]]$dex



    ## Set order of clusters
    if (resolution == "mres"){
        ## Set the identity to plot  
        cls_order <- switch(ident.use, 
                    All = c("L2_3_IT","L3_5_IT_1","L3_5_IT_2","L3_5_IT_3","L5_PT","L5_6_NP","L6_CT","L6_IT_1","L6_IT_2","L6b","LAMP5_LHX6","LAMP5_RELN","ADARB2","VIP","SST_NPY","SST","TH","PVALB","PVALB_CHC","Astro","OPC","Oligo","Micro","immune","blood","Endo","PC","SMC","VLMC"), 
                    ExN = c("L2_3_IT","L3_5_IT_1","L3_5_IT_2","L3_5_IT_3","L5_PT","L5_6_NP","L6_CT","L6_IT_1","L6_IT_2","L6b"), 
                    InN = c("LAMP5_LHX6","LAMP5_RELN","ADARB2","VIP","SST_NPY","SST","TH","PVALB","PVALB_CHC"), 
                    NNC = c("Astro","OPC","Oligo","Micro","immune","blood","Endo","PC","SMC","VLMC"))
        sf <- switch(ident.use, All = 1, ExN = 0.5, InN = 0.5, NNC = 0.5)
        font.sf <- 1
        height.base <- 1.5
    } else if (resolution == "hres"){
        all_cls <- readRDS(file = paste0(dataDir, "PFC.cluster.order.rds"))
        cls_order <- switch(ident.use, 
                    All = all_cls, 
                    ExN = all_cls[1:40],
                    InN = all_cls[41:86],
                    NNC = all_cls[87:115]) 
        sf <- switch(ident.use, All = 3, ExN = 1.6, InN = 1.6, NNC = 1)
        font.sf <- 0.8
        height.base <- 2.5
    } else {
        stop("Unsupported type of resolution")
    }
    

    ## Filter genes
    features <- ValidateFeatures(features = features, all.genes = rownames(avgs))



    plot_sp_heatmap(avg = avgs, features = features, scale.exp = TRUE, min_max = c(-2.5, 2.5), file_name = file_name, cluster.order = cls_order, font.scale = 1, height.scale = 1, width.scale = 1 * sf, return.plot = FALSE)


    ## scaled by size
    CirclePlot.horizontal(avg = avgs, ratio = ratios, features = features, dot.min = 0, file_name = paste0(file_name, "_bysize"), dot.scale = 4, scale.by = "size", cluster.order = cls_order, stroke.size = 0, stroke.matrix = NULL, width.scale = 1 * sf, height.base = height.base, font.scale = 1 * font.sf)
    CirclePlot.vertical(avg = avgs, ratio = ratios, features = features, dot.min = 0, file_name = paste0(file_name, "_bysize"), dot.scale = 4, scale.by = "size", cluster.order = cls_order, stroke.size = 0, stroke.matrix = NULL, width.scale = 1 * sf, height.base = height.base, font.scale = 1 * font.sf)


    ## scaled by radius
    CirclePlot.horizontal(avg = avgs, ratio = ratios, features = features, dot.min = 0, file_name = paste0(file_name, "_radius"), dot.scale = 4, scale.by = "radius", cluster.order = cls_order, stroke.size = 0, stroke.matrix = NULL, width.scale = 1 * sf, height.base = height.base, font.scale = 1 * font.sf)
    CirclePlot.vertical(avg = avgs, ratio = ratios, features = features, dot.min = 0, file_name = paste0(file_name, "_radius"), dot.scale = 4, scale.by = "radius", cluster.order = cls_order, stroke.size = 0, stroke.matrix = NULL, width.scale = 1 * sf, height.base = height.base, font.scale = 1 * font.sf)

}


