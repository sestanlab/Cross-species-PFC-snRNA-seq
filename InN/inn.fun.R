gset_enrich <- function(universe, interesting_genes, gset, min_nodeSize = 5){
    ## First filter the gset & interseting_genes based on the all_genes (gene universe)
    gfset <- lapply(gset, function(x) intersect(x, universe))
    gfset <- gfset[sapply(gfset, function(x) length(x) >= min_nodeSize)]
    int_genes <- intersect(interesting_genes, universe)



    ## Get the value of all_sig(a) & all_non_sig(b)
    a <- length(int_genes)
    b <- length(universe) - length(int_genes)


    ## For each gene set, get the scores (0, 1) for genes
    gf_res <- lapply(gfset, function(x) {
        ## contingency table
        sig <- sum(int_genes %in% x)
        non_sig <- length(x) - sig
        ctable <- c(sig, a - sig, non_sig, b - non_sig) %>% matrix(., nrow = 2, ncol = 2, byrow = TRUE)

        ## Fisher exact test
        pval <- fisher.test(ctable, alternative = "greater")$p.value
        pval
        })


    pval_res <- data.frame(GFname = names(gf_res), 
                        size = sapply(gfset, length), 
                        expec_hits = sapply(gfset, function(x) length(int_genes) * length(x) / length(universe)),
                        obser_hits = sapply(gfset, function(x) sum(int_genes %in% x)), 
                        pval = sapply(gf_res, function(x) x), 
                        stringsAsFactors = FALSE) %>% 
                            mutate(mlogp = -log10(pval)) %>%
                            mutate(FDR = p.adjust(pval, method = "fdr")) %>%
                            arrange(pval)


    return(pval_res)
}


## Plot Sankey plot
plot_ggsankey <- function(meta, x_col, next_col, x_ord = NULL, next_ord = NULL, type = "sankey", space= 0.01, width= 0.1, x_cols, next_cols) {
    let_ord <- paste0(rep(letters, each = 9), rep(1:9, times = length(letters)))
    if (is.null(x_ord)){
        x_ord <- pull(meta, x_col) %>% as.character() %>% as.factor() %>% levels()
    }
    x_trans <- paste0(let_ord[1:length(x_ord)], x_ord) %>%
                    setNames(., x_ord)
    meta$newx <- x_trans[as.character(meta[, x_col])]



    if (is.null(next_ord)){
        next_ord <- pull(meta, next_col) %>% as.character() %>% as.factor() %>% levels()
    }
    next_trans <- paste0(let_ord[1:length(next_ord)], next_ord) %>%
                    setNames(., next_ord)
    meta$newnext <- next_trans[as.character(meta[, next_col])]


    df <- meta %>%
                make_long(newx, newnext) %>%
                mutate(nlabel = substring(node, 3))


    cols <- c(x_cols, next_cols) %>%
                setNames(., c(setNames(x_trans, NULL), setNames(next_trans, NULL)))
    p <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = nlabel)) +
            geom_sankey(flow.alpha = .6, node.color = "gray30", type = type, width = width) +
            geom_sankey_text(size = 3, color = "black") +
            scale_fill_manual(values = cols) +
            theme_sankey(base_size = 18) +
            labs(x = NULL) +
            theme(legend.position = "none",  plot.title = element_text(hjust = .5), axis.text = element_blank())
    return(p)
}




PlotScatterPie2_custominterval <- function(pie.data, group.col = "cluster", feature.col = "gene", r.col = "radius", split.order = c("Human", "Chimpanzee", "Rhesus", "Marmoset"), group.order = NULL, feature.order = NULL, rsf = 2, scale.expression = TRUE, bg.col = "lightgrey", x_scale = 1, y_scale = 1){ ## rsf:radius-scale-factor
   ## set the colors 
   all.sp <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
   sp.colors <- setNames(c("#FF420E", "#4CB5F5", "#89DA59", "#FFBB00", "#D3D3D3", "#FFFFFF"), c(all.sp, "Shared", "empty")) %>%   
            .[split.order]


   ## Set the order of x & y Axis
   if (is.null(group.order)){
      group.order <- levels(as.factor(pie.data[, group.col]))
   }
   if (is.null(feature.order)){
      feature.order <- levels(as.factor(pie.data[, feature.col]))
   }

   pie.plot <- pie.data %>%
               mutate(!!group.col := as.character(!!sym(group.col))) %>%
               mutate(!!feature.col := as.character(!!sym(feature.col))) %>%
               mutate(!!group.col := as.numeric(factor(!!sym(group.col), levels = group.order))) %>%
               mutate(!!feature.col := as.numeric(factor(!!sym(feature.col), levels = rev(feature.order))))


    pie.plot[, group.col] <- pie.plot[, group.col] * x_scale
    pie.plot[, feature.col] <- pie.plot[, feature.col] * y_scale

   ## legend data
   lg.ra <- pie.plot[, r.col]
   lg.ra[lg.ra == 0] <- 0.00001




   ## Plots
   p <- ggplot(data = pie.plot) + 
      geom_scatterpie_new(mapping = aes_string(x = group.col, y = feature.col, r0 = "r0", r = r.col, amount = "value", fill = "fill_color", color = "border_color"), data = pie.plot, cols = sp.colors, scale.expression = scale.expression) + 
      theme_cowplot() + 
      coord_fixed(xlim = c(0, length(group.order)*x_scale + 4), ylim = c(0, length(feature.order)*y_scale + 1.5), expand = FALSE) + 
      scale_x_continuous(breaks = (1:length(group.order)) * x_scale, labels = group.order) +
      scale_y_continuous(breaks = (1:length(feature.order)) * y_scale, labels = rev(feature.order)) +
      ##scale_fill_manual(values = sp.colors)+
      scale_fill_identity() + 
      scale_color_identity() + 
      RotatedAxis() + 
      theme(panel.grid.major = element_line(colour=bg.col, size = 0.2), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.title = element_blank(), axis.text = element_text(size = 11)) +
      geom_scatterpie_legend(lg.ra, x=length(group.order)*x_scale+1, y=ceiling(length(feature.order)/4), n = 4, labeller = function(x) rsf * x)
   return(p)
}


