source("../scripts/pfc.fun.R")



avg <- readRDS(file = paste0(inputdir, "Subset.avg.rds"))
avg_sup <- readRDS(file = paste0(inputdir, "Subset.avg.sup.rds"))
avg <- cbind(avg, avg_sup[rownames(avg), , drop = FALSE])

hvg <- lapply(c("ExN", "InN", "NNC"), function(ctp) readRDS(file = paste0(inputdir, "HVG_CTP_Marker_", ctp, "_metrics.rds"))$gene) %>% 
            unlist() %>%
            unique()


## Pair-wise correlation of clusters
sp_pairs <- c("Human|Chimpanzee", "Human|Rhesus", "Human|Marmoset", "Chimpanzee|Rhesus", "Chimpanzee|Marmoset", "Rhesus|Marmoset")
cls_order <- readRDS(file = paste0(dataDir, "PFC.cluster.order.rds"))
cls_order <- setdiff(cls_order, c("L2-3 CUX2 ARHGAP18", "InN LAMP5 SYT10", "rAstro AQP4 OSMR","nOligo MOG CDH7", "Micro P2RY12 CCL3", "Micro P2RY12 GLDN", "B EBF1 IGKC", "cOPC PDGFRA TOP2A"))
cls_order <- cls_order[c(1:13, 15:21, 14, 22:length(cls_order))]




cor_res <- lapply(sp_pairs, function(pair) {
    sps <- strsplit(pair, "|", fixed = TRUE)[[1]]
    sp1 <- sps[1]; sp2 <- sps[2]

    cls1 <- paste0(sp1, "|", cls_order);
    cls2 <- paste0(sp2, "|", cls_order)

    aa <- cor(avg[hvg, cls1], avg[hvg, cls2], method = "pearson") 
    res <- diag(aa)%>%
            setNames(., cls_order)
    return(res)
    }) %>%
        setNames(., sp_pairs) %>%
        as.data.frame(., check.names = FALSE)
div_mat <- (1-cor_res) %>% as.matrix()




pair_order <- rev(c("Human|Chimpanzee", "Human|Rhesus", "Chimpanzee|Rhesus", "Rhesus|Marmoset", "Human|Marmoset", "Chimpanzee|Marmoset"))


col_use <- viridis(10)[c(1,2,5,6,7,8,10)]

gg_heat <- function(data, div_col = "div", cols) {
    p <- ggplot(data, aes_string(x = "cluster", y = "pair", fill = div_col)) +
                geom_tile(width = 1, height = 1, size = 0.1, color = "black") +
                scale_fill_gradientn(colors = cols) +
                theme_classic() +
                RotatedAxis() + 
                coord_fixed() +
                theme(legend.position = "right", axis.text.x = element_text(size = rel(0.5)), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.y = element_blank())
    return(p)
}





## Get the raw divergence scores
raw_data <- div_mat %>%
                reshape2::melt() %>%
                setNames(., c("cluster", "pair", "div")) %>%
                mutate(cluster = factor(as.character(cluster), levels = cls_order)) %>%
                mutate(pair = factor(as.character(pair), levels = pair_order))


raw_data$div <- MinMax(raw_data$div, min = 0, max = quantile(raw_data$div, 0.98))
p_raw <- gg_heat(data = raw_data, div_col = "div", cols = col_use) 




## Get the scaled divergence scores
scale_data <- div_mat %>%
        scale() %>%
        reshape2::melt() %>%
        setNames(., c("cluster", "pair", "div")) %>%
        mutate(cluster = factor(as.character(cluster), levels = cls_order)) %>%
        mutate(pair = factor(as.character(pair), levels = pair_order))

p_scale <- gg_heat(data = scale_data, div_col = "div", cols = col_use) 




p_raw <- p_raw + 
            theme(axis.text.x = element_blank())
p_scale <- p_scale + 
            theme(axis.text.x = element_blank())
pdf(paste0(outputdir, "DIV_1minusCor_Subset_combine_v2.pdf"), width = 10, height = 4)
plot_grid(p_raw, p_scale, nrow = 2, ncol = 1, align = "v") %>% print()
dev.off()


 




















