## Plot SF1 quality overview
source("../scripts/pfc.fun.R")
library(lemon)

###############################################################################################

##A. Plot the quality aross clusters

###############################################################################################
meta_use <- readRDS(file = paste0(dataDir, "PFC.filtered.meta.03022022.rds"))
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
            gsub("VLMCA8", "ABCA8",.) %>%
            gsub("L6b", "L6B",.)
    return(xx)
}



plot_cols <- c("mapped_reads", "nUMI", "nGene")
group.by <- "fig1cluster"
split.by <- "species"
rep.by <- "repname"
ind.by <- "samplename"
sp_colors <- c("#FF420E", "#4CB5F5", "#89DA59", "#FFBB00") %>% setNames(., c("Human", "Chimpanzee", "Rhesus", "Marmoset"))
cls_order <- readRDS(file = paste0(dataDir, "PFC.cluster.order.rds"))
cls_order <- updata_name(cls_order)
meta_use$fig1cluster <-  updata_name(as.character(meta_use$fig1cluster))
cluster_colors <- rep(c("#C51B7D", "#4D9221", "#984ea3", "#999999"), c(40, 46, 14, 15)) %>% setNames(., cls_order)



###----------------------------------------------------------------------------------------
## Plot nUMI, nGene, mapped reads by cluster
plot_data <- meta_use[, c(plot_cols, group.by, split.by)] %>%
                reshape2::melt(id.vars = c(group.by, split.by), measure.vars = plot_cols, variable.name = "quality") %>%
                mutate(!!group.by := as.character(!!sym(group.by))) %>%
                mutate(!!split.by := as.character(!!sym(split.by))) %>%
                mutate(!!group.by := factor(!!sym(group.by), levels = rev(cls_order))) %>%
                mutate(!!split.by := factor(!!sym(split.by), levels = rev(names(sp_colors))))

plot_data$value <- plot_data$value/1000
plot_data$quality <- plot_data$quality %>%
                        gsub("mapped_reads", "nReads ( * 1000)", .) %>%
                        gsub("nUMI", "nUMIs ( * 1000)", .) %>%
                        gsub("nGene", "nGenes ( * 1000)", .)
plot_data$quality <- factor(plot_data$quality, levels = c("nReads ( * 1000)", "nUMIs ( * 1000)", "nGenes ( * 1000)"))


p <- ggplot(plot_data, aes_string(y = group.by, x = "value")) +
        geom_violin(aes_string(fill = split.by), scale = "width", size = 0.05, adjust = 1.5, trim =TRUE, alpha = 0.8) + 
        geom_boxplot(width=0.3, outlier.shape = NA, lwd= 0.05, alpha = 0.3) + 
        coord_capped_cart(top='both', left='both') +
        scale_fill_manual(values = sp_colors) +
        theme_classic() + 
        scale_x_continuous(position = "top") + 
        theme(panel.border=element_blank(), axis.line=element_line(size = 0.25), axis.ticks = element_line(size = 0.25),axis.text.x = element_text(size = 10, hjust = 0, angle = 45), axis.text.y = element_text(size = 9, colour = rev(cluster_colors)), axis.title = element_blank()) + 
        facet_wrap(facets = vars(quality), nrow = 1, ncol = length(plot_cols), strip.position = "top", scales = "free_x") + 
        theme(strip.background = element_blank(), strip.text = element_text(size = 11, face = "bold"), panel.spacing = unit(0.05, "in"), legend.title = element_text(size = 10), legend.text = element_text(size = 12), strip.placement = "outside")
pdf(paste0(outputdir, "Cluster_quality.pdf"), width = 7.5, height = 12)##, units = "in",  res = 300)
print(p)
dev.off()


size_data <- meta_use %>%
        group_by(!!sym(split.by), !!sym(group.by)) %>%
        summarise(size = n()) %>%
        ungroup() %>%
        mutate(size = sqrt(size)) %>%
        mutate(!!group.by := factor(!!sym(group.by), levels = rev(cls_order))) %>%
        mutate(!!split.by := factor(!!sym(split.by), levels = rev(names(sp_colors))))
size_max <- size_data %>%
        group_by(!!sym(group.by)) %>%
        summarise(size = sum(size)) %>%
        .$size %>% max()
size_max <- ceiling(sqrt(size_max)/20)*20
size_bks <- seq(0, size_max, length.out = 5)



q <- ggplot(data = size_data, mapping = aes_string(x = "size", y = group.by, color = split.by, fill = split.by)) +
        geom_bar(stat = "identity", size = 0.05, position = "stack") + 
        scale_fill_manual(values = sp_colors) +
        scale_color_manual(values = sp_colors) +
        theme_classic() + 
        scale_x_continuous(name = expression("size"), breaks = size_bks, labels = as.character(round(size_bks^2, digits = 0)), limits = c(0, size_max))
pdf(paste0(outputdir, "Cluster_size_v3.pdf"), width = 8, height = 12)
print(q)
dev.off()



new_p <- p + theme(legend.position = "bottom", axis.text.y = element_blank())
new_q <- q + theme(legend.position = "bottom")
cluster_p <- plot_grid(new_q, new_p, nrow = 1, ncol = 2, align = "h", rel_widths = c(1, 1.25))

pdf(paste0(outputdir, "Cluster_plots.pdf"), width = 12, height = 14)##, units = "in",  res = 300)
cluster_p %>% print()
dev.off()




###----------------------------------------------------------------------------------------
## Plot number of cells by samples
plot_data <- meta_use %>%
                rownames_to_column("cell") %>%
                group_by(!!sym(rep.by), !!sym(split.by)) %>%
                summarize(size = n(), !!ind.by := unique(!!sym(ind.by))) %>%
                ungroup() %>% 
                mutate(!!split.by := factor(!!sym(split.by), levels = names(sp_colors))) %>%
                mutate(!!rep.by := extract_field(!!sym(rep.by), -1, "_")) %>%
                mutate(!!rep.by := as.numeric(as.factor(!!sym(rep.by)))) %>%
                mutate(!!rep.by := paste0("REP-", !!sym(rep.by))) 
sp_ncells <- plot_data %>% 
                group_by(!!sym(split.by)) %>%
                summarize(size = sum(size)) %>%
                mutate(y_loc = 65000) %>% 
                mutate(!!split.by := factor(!!sym(split.by), levels = names(sp_colors))) %>%
                mutate(cell_origin = c("HSB189", "PTB166", "RMB196", "CJB1540"))

all_samples <- plot_data[, ind.by][[1]] %>% unique()
sample_order <- all_samples[grep("HSB", all_samples) %>% c(., grep("PTB", all_samples)) %>% c(., grep("RMB", all_samples)) %>% c(., grep("CJB", all_samples))]


ind.cells <- plot_data %>% 
                group_by(!!sym(ind.by)) %>%
                summarize(size = sum(size), species = unique(species)) 


mm <- ggplot(plot_data, aes_string(x = ind.by, y = "size")) +
        geom_bar(aes_string(x = ind.by, y = "size", fill = split.by), color = "black", position = position_stack(reverse = FALSE), stat = "identity", lwd = 0.3) +
        geom_text(data = ind.cells, mapping = aes_string(x = ind.by, y = "size", label = "size"), nudge_y = 1500, angle = 45, vjust = 0.5, hjust = 0.1, size = 2.8) +
        coord_capped_cart(bottom='both', left='both') +
        scale_fill_manual(values = sp_colors) +
        scale_x_discrete(limits = sample_order) +  
        geom_label(data = sp_ncells, aes_string(x = "cell_origin", y = "y_loc", label = "size", fill = split.by), nudge_x = 0.5) + 
        
        theme_classic() + 
        RotatedAxis() + 
        labs(y = "Sample size", x = "Individual") +
        theme(axis.line=element_line(size = 0.25), axis.ticks=element_line(size = 0.25))


pdf(paste0(outputdir, "Sample_size.pdf"), width = 7, height = 5)
print(mm)
dev.off()





###----------------------------------------------------------------------------------------
## Plot quality by samples
plot_data <- meta_use[, c(plot_cols, split.by, ind.by)] %>%
                reshape2::melt(id.vars = c(split.by, ind.by), measure.vars = plot_cols, variable.name = "quality") %>%
                mutate(!!split.by := factor(!!sym(split.by), levels = names(sp_colors)))

plot_data$value <- plot_data$value/1000
plot_data$quality <- plot_data$quality %>%
                        gsub("mapped_reads", "nReads ( * 1000)", .) %>%
                        gsub("nUMI", "nUMIs ( * 1000)", .) %>%
                        gsub("nGene", "nGenes ( * 1000)", .)
plot_data$quality <- factor(plot_data$quality, levels = c("nReads ( * 1000)", "nUMIs ( * 1000)", "nGenes ( * 1000)"))



nn <- ggplot(plot_data, aes_string(x = ind.by, y = "value", fill = split.by)) +
        geom_violin(scale = "width", size = 0.1,adjust = 2,trim =TRUE) + 
        geom_boxplot(fill = "white", color = "black", width=0.15, outlier.shape = NA, lwd= 0.2, alpha = 0.5) + 
        coord_capped_cart(top='both', left='both') +
        scale_fill_manual(values = sp_colors) +
        theme_classic() + 
        scale_x_discrete(limits = sample_order) +
        ##RotatedAxis() + 
        theme(panel.border=element_blank(), axis.line=element_line(size = 0.25), axis.ticks = element_line(size = 0.25),axis.text.x = element_text(size = 10, hjust = 1, angle = 45), axis.text.y = element_text(size = 9), axis.title = element_blank()) + 
        facet_wrap(facets = vars(quality), nrow = length(plot_cols), ncol = 1, strip.position = "left", scales = "free_y") + 
        theme(strip.background = element_blank(), strip.text = element_text(size = 8, face = "bold"), panel.spacing = unit(0.05, "in"), legend.title = element_text(size = 8), legend.text = element_text(size = 12), strip.placement = "outside")
pdf(paste0(outputdir, "Sample_quality.pdf"), width = 6, height = 5)
print(nn)
dev.off()




sample_p <- plot_grid(mm, nn, nrow = 1, ncol = 2, rel_widths = c(1.1, 1), scale = c(1, 0.95))

pdf(paste0(outputdir, "Sample_Plots.pdf"), width = 13, height = 4)##, units = "in",  res = 300)
sample_p %>% print()
dev.off()



pdf(paste0(outputdir, "SF1_V1.pdf"), width = 13, height = 17)
plot_grid(sample_p, cluster_p, nrow = 2, ncol = 1, rel_heights = c(1, 3.5)) %>% print()
dev.off()





