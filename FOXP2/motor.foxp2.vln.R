source("../scripts/pfc.fun.R")


## Load data
all_sps <- c("Human", "Marmoset", "Mouse")
seu_list <- lapply(all_sps, function(sp) {
    readRDS(file = paste0("./motor_cortex/", "Motor.", sp, ".full.rds"))
    }) %>%
        setNames(., all_sps)


gene <- 'FOXP2'

## Extract FOXP2 expression & cell identity
pdata <- lapply(names(seu_list), function(sp) {
    gene <- grep(paste0("^", gene, "$"), rownames(seu_list[[sp]]$RNA@data), value = TRUE, ignore.case = TRUE)
    exp <- seu_list[[sp]]$RNA@data[gene, ]
    df <- data.frame(exp = exp,
                hres = as.character(seu_list[[sp]]@meta.data$cluster_label),
                cluster = as.character(seu_list[[sp]]@meta.data$subclass_label),
                species = sp,
                stringsAsFactors = FALSE)
    df
    }) %>%
        do.call(rbind, .)


## Remove certain clusters
pdata <- pdata %>%
            filter(!is.na(cluster)) %>%
            filter(!cluster %in% "Meis2")

pdata$cluster[pdata$cluster %in% c("SMC","Peri","Endo","VLMC")] <- "Vas"
pdata$cluster[pdata$cluster %in% c("Sst","Sst Chodl")] <- "Sst"
pdata$cluster[pdata$cluster %in% c("L6 IT","L6 IT Car3")] <- "L6 IT"
pdata$cluster[pdata$hres %in% c("PVM_1","PVM_2","PVM_3")] <- "immune"
pdata$cluster[pdata$cluster %in% c("Micro-PVM")] <- "Micro"


cls_ord <- c("L2/3 IT","L5 IT","L5 ET","L5/6 NP","L6 CT","L6 IT","L6b","Lamp5","Sncg","Vip","Pvalb","Sst","Astro","OPC","Oligo","Micro","immune","Vas")
pdata$cluster <- factor(pdata$cluster, levels = cls_ord)


## Plot Type
me <- max(pdata$exp)
p <- ggplot() + 
        geom_violin(data = pdata, aes_string(x = "cluster", y = "exp", fill = "species"), scale = "width", size = 0.1, adjust = 2,trim =TRUE) + 
        theme_cowplot() + 
        RotatedAxis() + 
        scale_fill_manual(values = c("#FF420E", "#FFBB00", "#AE017E") %>% setNames(., c("Human", "Marmoset", "Mouse"))) + 
        facet_wrap(vars(species), nrow = 3, ncol = 1, strip.position = 'left') +
        scale_y_continuous(breaks = seq(0, me, by = round(me/2, digits = 1))) + 
        theme(axis.title = element_blank(), axis.text.x=element_text(size = rel(0.9)), strip.background = element_blank(), strip.text = element_text(size = rel(0.8), angle = 90), panel.spacing = unit(0.05, "in"), strip.placement = 'outside', axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), legend.position = "bottom", legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size = rel(0.5)), legend.title = element_blank())
pdf(paste0(outputdir, "Motor.HMM.", gene, ".pdf"), width = 6, height = 4)
print(p)
dev.off() 






