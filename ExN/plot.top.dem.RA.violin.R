source("../scripts/pfc.fun.R") 


## Plot the RA gene expression using violin plots
pfc <- readRDS(file = paste0(dataDir, "Final_PFC_HPRC.All.07102021.rds"))


cls_order <- c("L2_3_IT","L3_5_IT_1","L3_5_IT_2","L3_5_IT_3","L5_PT","L5_6_NP","L6_CT","L6_IT_1","L6_IT_2","L6b", "LAMP5_LHX6", "LAMP5_RELN", "VIP", "ADARB2", "SST_NPY", "SST", "TH", "PVALB","PVALB_CHC", "Astro","OPC","Oligo","Micro","immune","blood","Endo","PC","SMC","VLMC")


features <- c("ALDH1A1", "ALDH1A2", "ALDH1A3", "RARA", "RARB", "RARG", "RXRA", "RXRB", "RXRG", "CYP26A1", "CYP26B1", "CYP26C1", "CBLN2", "MEIS2")
features <- intersect(features, rownames(pfc))



## Do pair-wise comparison across species for each cell types

spdf2arry <- function(df, all_pairs) {
    ## Add empty cols for missing pairs (diagonal Human-Human, Chimp-Chimp, etc)
    mis_pairs <- setdiff(all_pairs, colnames(df))
    if (length(mis_pairs) > 0){
        mis_df <- matrix(0, nrow = nrow(df), ncol = length(mis_pairs), dimnames = list(rownames(df), mis_pairs))
        df <- cbind(as.matrix(df), mis_df)
    }
    df <- df[, all_pairs] %>%
            as.matrix()

    dexinfo <- lapply(1:nrow(df), function(i) {matrix(df[i, ], nrow = 4, ncol = 4, byrow = TRUE)}) %>%
            do.call(c, .) %>%
            array(., dim = c(4, 4, nrow(df)), dimnames = list(all_sps, all_sps, rownames(df)))
    return(dexinfo)
}


SummariseArray <- function(ary, all_sps) {
    genes <- dimnames(ary)[[3]]
    rSUM <- lapply(genes, function(x) rowSums(ary[, , x])) %>%
            setNames(., genes) %>%
            as.data.frame(., check.names = FALSE) %>%
            as.matrix() %>% t()
    cSUM <- lapply(genes, function(x) colSums(ary[, , x])) %>%
            setNames(., genes) %>%
            as.data.frame(., check.names = FALSE) %>%
            as.matrix() %>% t()


    ept_mat <- matrix(0, nrow = length(genes), ncol = length(all_sps), dimnames = list(genes, all_sps))

    ##------------------------------------------------
    mod3_idx <- which(apply(cSUM, 1, function(x) sum(x == 3) == 1))
    tem_mat <- cSUM[mod3_idx, ]
    tem_mat[tem_mat != 3] <- 1
    tem_mat[tem_mat == 3] <- 0
    ept_mat[mod3_idx, ] <- tem_mat


    mod1_idx <- which(apply(rSUM, 1, function(x) sum(x == 3) == 1))
    tem_mat <- rSUM[mod1_idx, ]
    tem_mat[tem_mat != 3] <- 0
    tem_mat[tem_mat == 3] <- 1
    ept_mat[mod1_idx, ] <- tem_mat


    mod2_idx <- which(apply(cSUM, 1, function(x) sum(x == 2) == 2) & 
                        apply(rSUM, 1, function(x) sum(x == 2) == 2))
    tem_mat <- rSUM[mod2_idx, ]
    tem_mat[tem_mat != 2] <- 0
    tem_mat[tem_mat == 2] <- 1
    ept_mat[mod2_idx, ] <- tem_mat
    return(ept_mat)
}



all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
pair_ord <- c("Human_Chimpanzee", "Human_Rhesus", "Human_Marmoset",
            "Chimpanzee_Human", "Chimpanzee_Rhesus", "Chimpanzee_Marmoset", 
            "Rhesus_Human", "Rhesus_Chimpanzee", "Rhesus_Marmoset",
            "Marmoset_Human", "Marmoset_Chimpanzee", "Marmoset_Rhesus")
all_pairs <- c("Human_Human", pair_ord[1:4], "Chimpanzee_Chimpanzee", pair_ord[5:8],
                "Rhesus_Rhesus", pair_ord[9:12], "Marmoset_Marmoset")


dex_res <- readRDS(file = paste0(dataDir, "DEX.hprc.mres.rds"))
slim_dex <- dex_res %>%
                    filter(pct.1 >= 0.25 & ratio_fc >= 1.5 & avg_logFC >= 0.35) %>%
                    group_by(cluster, speciespair) %>% 
                    mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
                    filter(p_val_adj <= 0.01) %>%
                    ungroup() %>%
                    subset(cluster %in% cls_order)

res_cluster <- split(slim_dex, slim_dex$cluster)
dex_sum <- lapply(res_cluster, function(x) {
                dexi <- split(x$gene, x$speciespair) %>%
                            lapply(., function(y) setNames(as.numeric(features %in% y), features)) %>%
                            as.data.frame(., check.names = FALSE) %>%
                            as.matrix() %>%
                            spdf2arry(df = ., all_pairs = all_pairs)
                dexi
                })


sp_cols <- setNames(c("#FF420E","#4CB5F5","#89DA59","#FFBB00"), c("Human", "Chimpanzee", "Rhesus", "Marmoset"))
dex_df <- lapply(names(dex_sum), function(cls) {
    subdata <- dex_sum[[cls]]
        df <- SummariseArray(ary = subdata, all_sps = all_sps) %>%
                reshape2::melt() %>%
                setNames(., c("gene", "species", "value")) %>%
                mutate(gene = as.character(gene), species = as.character(species)) %>%
                mutate(cluster = cls) %>%
                mutate(color = ifelse(value == 1, sp_cols[species], "#bdbdbd"))

    df
    }) %>% 
        do.call(rbind, .)


col_vec <- setNames(dex_df$color, paste0(dex_df$gene, "|", dex_df$species, "|", dex_df$cluster))


data <- lapply(features, function(gene) {
    print(gene)
    pdata <- data.frame(cluster = pfc@meta.data$mres,
                        species = pfc@meta.data$species,
                        exp = pfc$RNA@data[gene, ],
                        gene = gene, 
                        stringsAsFactors = FALSE) %>%
                    mutate(genespcls = paste0(gene, "|", species, "|", cluster))
    pdata
    }) %>%
        do.call(rbind, .)  %>%
        mutate(cluster = factor(as.character(cluster), levels = cls_order)) %>%
        mutate(genespcls = factor(as.character(genespcls), levels = names(col_vec))) %>%
        mutate(gene = factor(as.character(gene), levels = features))




p <- ggplot(data, aes_string(x = "cluster", y = "exp", fill = "genespcls", color = "genespcls")) + 
                    geom_violin(scale = "width", size = 0.001, adjust = 1,trim =TRUE, position = position_dodge(0.75)) + 
                    theme_cowplot() + 
                    RotatedAxis() + 
                    scale_fill_manual(values = col_vec) +
                    scale_color_manual(values = setNames(ifelse(col_vec == "#bdbdbd", NA, col_vec), names(col_vec))) +
                    facet_wrap(vars(gene), nrow = length(features), ncol = 1, strip.position = 'left', scales = "free_y") +
                    theme(axis.title = element_blank(), axis.text.x=element_text(size = rel(0.9)), axis.text.y=element_text(size = rel(0.7)), strip.background = element_blank(), strip.text = element_text(size = rel(0.7), angle = 90), panel.spacing = unit(0.05, "in"), strip.placement = 'outside', axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), legend.position = "none", legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size = rel(0.5)), legend.title = element_blank(), plot.margin = unit(c(-0.03, 0, -0.03, 0), "in"))

pdf(paste0(outputdir, "RA.violin.allcelltypes.highlight.pdf"), width = 12, height = 8)
print(p)
dev.off() 

















