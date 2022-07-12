source("../scripts/pfc.fun.R")


all_cls <- c("L2_3_IT", "L3_5_IT", "L56", "MGE", "CGE", "Astro", "Immune", "Oligo", "Vas")
all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")

## Rename cluster names
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


dend <- readRDS(file = paste0("../MF1_overview/load_files/MF1_plot_dend_v2.rds"))
cls_ord <- rev(labels(dend))



all_dfs <- lapply(all_cls, function(ctp) {
    sp_mats <- lapply(all_sps, function(sp) {
        aa <- readRDS(file = paste0(inputdir, "RobustAUC_", sp, "_", ctp, ".rds")) %>%
        		as.matrix() %>%
				reshape2::melt() %>%
				setNames(., c("cluster", "samplename", "AUC")) %>%
				mutate(cluster = as.character(cluster), samplename = as.character(samplename)) %>%
				mutate(species = sp) %>%
				mutate(AUC = MinMax(AUC, min = 0.5, max = 1))
        aa
        }) %>%
            do.call(rbind, .)
    sp_mats
    }) %>%
        do.call(rbind, .) %>%
        mutate(cluster = updata_name(cluster)) %>%
        mutate(cluster = factor(cluster, levels = cls_ord)) %>%
        mutate(species = factor(species, levels = all_sps))

all_dfs$group <- "NNC"
all_dfs$group[grep("^L[0-9]", all_dfs$cluster)] <- "ExN"
all_dfs$group[grep("^InN", all_dfs$cluster)] <- "InN"
all_dfs$group[grep("Astro|OPC|COP|Oligo|Micro", all_dfs$cluster)] <- "Glia"


gp_cols <- setNames(c("#C51B7D", "#4D9221", "#08519c", "#999999"), c("ExN", "InN", "Glia", "NNC"))




all_dfs <- all_dfs %>% mutate(cluster = factor(as.character(cluster), levels = rev(cls_ord)))
p <- ggplot(all_dfs, aes_string(y = "cluster", x = "AUC", fill = "group", color = "group")) +
            geom_boxplot(width = 0.8, size = 0.2, outlier.shape = NA, alpha = 0.7) +
            scale_fill_manual(values = gp_cols)+
            scale_color_manual(values = gp_cols)+
            theme_classic() +
            scale_x_continuous(limits = c(0.5, 1)) +
            facet_wrap(vars(species), nrow = 1, ncol = 4, strip.position = "top") +            
            theme(legend.position = "right", strip.placement = 'outside', axis.text.x = element_text(size = rel(0.5)), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.title = element_blank(), axis.text.y = element_text(size = rel(0.8)))

pdf(paste0(outputdir, "Cluster_robustness_AUC_box_v2.pdf"), width = 5, height = 10)
print(p)
dev.off()



