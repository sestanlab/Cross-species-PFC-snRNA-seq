source("../scripts/pfc.fun.R")


resolution <- "mres"
all.markers <- readRDS(file = paste0(dataDir, "Markers.hprc.", resolution, ".rds"))
all.dex <- readRDS(file = paste0(dataDir, "DEX.hprc.", resolution, ".rds")) 
all.genes <- readRDS(paste0(dataDir, "Expr_avg_ratio_by_all.rds"))[[resolution]]$avg %>% rownames()


cls.order <- c("L2_3_IT","L3_5_IT_1","L3_5_IT_2","L3_5_IT_3","L5_PT","L5_6_NP","L6_CT","L6_IT_1","L6_IT_2","L6b")
cls.infor <- data.frame(row.names = cls.order, 
                        cluster = cls.order,
                        group = rep("ExN", 10),
                        stringsAsFactors = FALSE)
all_species <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")



## Get the DEM
slim.markers <- all.markers %>%
                        subset(pct.1 >= 0.25 & pct.2 <= 0.15 & ratio_fc >= 2.5 & avg_logFC >= 1.25) %>%
                        group_by(cluster, species) %>% 
                        mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
                        filter(p_val_adj <= 1e-4) %>%
                        ungroup() %>%
                        subset(cluster %in% cls.order)
slim.dex <- all.dex %>%
                        subset(pct.1 >= 0.25 & pct.2 <= 0.15 & ratio_fc >= 2.5 & avg_logFC >= 1.25) %>%
                        group_by(cluster, speciespair) %>% 
                        mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
                        filter(p_val_adj <= 1e-4) %>%
                        ungroup() %>%
                        subset(cluster %in% cls.order)


demarkers_idx <- lapply(cls.order, function(cls) {
        cls_idx <- lapply(all_species, function(sp) {
            mar <- slim.markers$gene[slim.markers$species == sp & slim.markers$cluster == cls]
            idx <- which(slim.dex$species1 == sp & slim.dex$gene %in% mar & slim.dex$cluster == cls)
            idx
            }) %>% unlist() %>% unique()
        }) %>% unlist() %>% sort() %>% unique()

demarkers <- slim.dex[demarkers_idx, ]




## Find the Up & Down DEM
empty_df <- data.frame(row.names = all.genes)
exl_up <- lapply(all_species, function(sp) {
        sp_exl <- lapply(cls.order, function(cls) {
                    df <- empty_df
                    degenes <- demarkers$gene[demarkers$species1 == sp & demarkers$cluster == cls] %>%
                                    table() %>% .[. == 3] %>% names() %>%
                                    grep("^RP[0-9].*-", ., value = TRUE, invert = TRUE) %>%
                                    grep("^CTD-", ., value = TRUE, invert = TRUE) %>%
                                    grep("^CTB-", ., value = TRUE, invert = TRUE) %>%
                                    grep("^AC00", ., value = TRUE, invert = TRUE)  %>%
                                    grep("^AP00", ., value = TRUE, invert = TRUE)
                    df[, cls] <- 0
                    df[degenes, cls] <- 1
                    return(df)
                    }) %>% do.call(cbind, .)
        sp_exl <- sp_exl[rowSums(sp_exl) > 0, ]

        ## Order the genes
        gord <- order(apply(sp_exl, 1, function(x) min(which(x != 0))))
        sp_exl <- sp_exl[gord, ]
        sp_exl
        }) %>% setNames(., all_species)
exl_down <- lapply(all_species, function(sp) {
        sp_exl <- lapply(cls.order, function(cls) {
                    df <- empty_df
                    degenes <- demarkers$gene[demarkers$species2 == sp & demarkers$cluster == cls] %>%
                                    table() %>% .[. == 3] %>% names() %>%
                                    grep("^RP[0-9].*-", ., value = TRUE, invert = TRUE) %>%
                                    grep("^CTD-", ., value = TRUE, invert = TRUE) %>%
                                    grep("^CTB-", ., value = TRUE, invert = TRUE) %>%
                                    grep("^AC00", ., value = TRUE, invert = TRUE)  %>%
                                    grep("^AP00", ., value = TRUE, invert = TRUE)
                    df[, cls] <- 0
                    df[degenes, cls] <- 1
                    return(df)
                    }) %>% do.call(cbind, .)
        sp_exl <- sp_exl[rowSums(sp_exl) > 0, ]
        
        ## Order the genes
        gord <- order(apply(sp_exl, 1, function(x) min(which(x != 0))))
        sp_exl <- sp_exl[gord, ]
        sp_exl
        }) %>% setNames(., all_species)
save(exl_up, exl_down, file = paste0(inputdir, "HPRC.ExN.topDEM.mres.Rdata"))


