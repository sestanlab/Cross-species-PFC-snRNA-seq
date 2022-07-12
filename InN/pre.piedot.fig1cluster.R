##-------------------------------------------------------------------
## Generate the data for species-pie plots
source("../scripts/pfc.fun.R")


##---------------------------------------------------------------------
## Get the species-enrichment for each cluster
resolution <- "fig1cluster"
# cluster order
homo.order <- readRDS(file = paste0(dataDir, "cluster.color.rds")) %>% rownames() %>%
                .[c(41:42, 44:86)]


avg_res <- readRDS(file = paste0(dataDir, "Expr_avg_ratio_by_species.rds"))[[resolution]]
dexres <- readRDS(file = paste0(dataDir, "DEX.hprc.", resolution, ".rds")) 
slim_dex <- dexres %>%
                    filter(pct.1 >= 0.2 & pct.2 < 0.35 & ratio_fc >= 1.5 & avg_logFC >= 0.25) %>%
                    group_by(cluster, speciespair) %>% 
                    mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
                    filter(p_val_adj <= 0.01) %>%
                    ungroup() %>%
                    subset(cluster %in% homo.order)

all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
all_genes <- rownames(avg_res$avg)
pair_ord <- c("Human_Chimpanzee", "Human_Rhesus", "Human_Marmoset",
            "Chimpanzee_Human", "Chimpanzee_Rhesus", "Chimpanzee_Marmoset", 
            "Rhesus_Human", "Rhesus_Chimpanzee", "Rhesus_Marmoset",
            "Marmoset_Human", "Marmoset_Chimpanzee", "Marmoset_Rhesus")
all_pairs <- c("Human_Human", pair_ord[1:4], "Chimpanzee_Chimpanzee", pair_ord[5:8],
                "Rhesus_Rhesus", pair_ord[9:12], "Marmoset_Marmoset")



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



## Organize the DEX results into array (sp * sp * genes) for each cluster
res_cluster <- split(slim_dex, slim_dex$cluster)
dex_sum <- lapply(res_cluster, function(x) {
                dexi <- split(x$gene, x$speciespair) %>%
                            lapply(., function(y) setNames(as.numeric(all_genes %in% y), all_genes)) %>%
                            as.data.frame(., check.names = FALSE) %>%
                            as.matrix() %>%
                            spdf2arry(df = ., all_pairs = all_pairs)
                dexi
                })



## Organize the results in the following format for pie-dot plots
## cluster, Human, Chimpanzee, Rhesus, Marmoset, gene, radius
piei <- lapply(names(dex_sum), function(cls) {
    message(paste0("Working on cluster: "), cls)
    data <- dex_sum[[cls]]
    ept_mat <- SummariseArray(ary = data, all_sps = all_sps) 
    isshare <- rowSums(ept_mat) == 0
    pre_df <- ept_mat %>%
                as.data.frame(., check.names = FALSE) %>%
                rownames_to_column("gene") %>%
                mutate(Shared = ifelse(isshare, 1, 0))


    ## Get the radius information
    sp_cls <- paste0(all_sps, "|", cls)
    cls_ratio <- avg_res$ratio[all_genes, sp_cls] %>% as.matrix()
    colnames(cls_ratio) <- all_sps

    pre_df$radius[isshare] <- apply(cls_ratio[which(isshare), ], 1, mean)
    pre_df$radius[!isshare] <- sapply(which(!isshare), function(x) {
        y <- cls_ratio[x, ] * ept_mat[x, ]
        y <- y[y!= 0]
        yy <- mean(y)
        return(yy)
        })

    pre_df$cluster <- cls
    return(pre_df)
    }) %>%
        do.call(rbind, .)%>%
        .[, c("cluster", all_sps, "Shared", "gene", "radius")]
saveRDS(piei, file = paste0(inputdir, "Pie.dot.plot.fig1cluster.InN.rds"))












                
