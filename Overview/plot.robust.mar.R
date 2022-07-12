source("../scripts/pfc.fun.R") 


all_cls <- c("L2_3_IT", "L6_IT", "L3_5_IT", "L5_6", "LAMP5", "SST", "PVALB", "ADARB2", "Astro", "immune", "Oligo", "Vas")
all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")

all_list <- lapply(all_cls, function(ctp) {
    sp_list <- lapply(all_sps, function(sp) {
        res_list <- readRDS(file = paste0("../MF1_overview/load_files/", "Cluster_Marres_", ctp, "_", sp, ".rds"))
        }) %>%
        setNames(., all_sps)
    sp_list
    }) %>%
        setNames(., all_cls)

new_list <- lapply(all_sps, function(sp) {
    mar_list <- lapply(all_list, function(x) x[[sp]]) %>%
            do.call(c, .)

    nmars <- sapply(mar_list, function(res){
        sum(res$ratio_fc >= 1.2 & res$p_val_adj <= 0.01 & res$avg_logFC >= 0.5)
        })
    nmars
    }) %>%
        setNames(., all_sps)




##  Get the plot data
cls_ord <- readRDS(file = paste0(dataDir, "PFC.cluster.order.rds"))


pdata <- lapply(new_list, function(x) {
    names(x) <- extract_field(names(x), 2, ".") %>% setNames(., NULL)
    emp <- setNames(rep(NA, length(cls_ord)), cls_ord)
    emp[names(x)] <- x
    emp
    }) %>%
        as.data.frame(., check.names = FALSE) %>%
        as.matrix() %>%
        reshape2::melt() %>%
        setNames(., c("cluster", "species", "nmar")) %>%
        mutate(species = factor(as.character(species), levels = rev(all_sps))) %>%
        mutate(nmar = MinMax(nmar, min = 0, max = 50))



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
            gsub("L6b", "L6B",.)
    return(xx)
}
dend <- readRDS(file = paste0("../MF1_overview/load_files/MF1_plot_dend_v2.rds"))
cls_ord <- rev(labels(dend))


pdata <- pdata %>%
			mutate(cluster = updata_name(cluster)) %>%
            mutate(cluster = factor(as.character(cluster), levels = cls_ord))



colors <- colorRampPalette(colors = c("white", "darkred"))(6)
colors[1] <- "#41ae76"
pdata$nmar_range <- cut(pdata$nmar, breaks = c(0, 0.9, 5.1, 10, 20, 50, Inf), include.lowest = TRUE)


p <- ggplot(pdata, aes_string(x = "cluster", y = "species", fill = "nmar_range")) +
                geom_tile(width = 1, height = 1, size = 0.1, color = "black") +
                scale_fill_manual(values = setNames(colors, levels(pdata$nmar_range)), na.value = "grey") +
                theme_classic() +
                RotatedAxis() + 
                coord_fixed() +
                theme(legend.position = "right", axis.text.x = element_text(size = rel(0.5)), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.y = element_blank())

pdf(paste0(outputdir, "Cluster_robustness_numbermarkers_v1.pdf"), width = 10, height = 5)
print(p)
dev.off()


