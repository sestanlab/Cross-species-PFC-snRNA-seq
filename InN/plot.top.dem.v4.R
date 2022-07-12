source("../scripts/pfc.fun.R")
source("../MF6_InN/inn.fun.R")


##------------------------------------------------------------------
# load data
resolution <- "fig1cluster"
all.markers <- readRDS(file = paste0(dataDir, "Markers.hprc.", resolution, ".rds"))
all.dex <- readRDS(file = paste0(dataDir, "DEX.hprc.", resolution, ".rds")) 
all_species <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")

# cluster order
homo.order <- readRDS(file = paste0(dataDir, "cluster.color.rds")) %>% rownames() %>%
                .[c(41:42, 44:86)]


load(file = "~/project/public_data/HGNC/HGNC_gset.09152021.Rdata") 


##------------------------------------------------------------------
## Get the dex & markers (DEM) genes
exp_ratio <- readRDS(file = paste0(dataDir, "Expr_avg_ratio_by_species.rds"))$fig1cluster$ratio
slim.dex <- all.dex %>%
                    filter(pct.1 >= 0.2 & pct.2 < 0.35 & ratio_fc >= 1.5 & avg_logFC >= 0.25) %>%
                    group_by(cluster, speciespair) %>% 
                    mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
                    filter(p_val_adj <= 0.01) %>%
                    ungroup()
slim.mar <- all.markers %>%
                    subset(pct.1 >= 0.2 & ratio_fc >= 1.25 & avg_logFC >= 0.25) %>%
                    group_by(cluster, species) %>% 
                    mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
                    filter(p_val_adj <= 0.01) %>%
                    ungroup()
demarkers_idx <- lapply(homo.order, function(cls) {
    cls_idx <- lapply(all_species, function(sp) {
        mar <- slim.mar$gene[slim.mar$species == sp & slim.mar$cluster == cls]
        idx <- which(slim.dex$species1 == sp & slim.dex$gene %in% mar & slim.dex$cluster == cls)
        idx
        }) %>% unlist() %>% unique()
    }) %>% unlist() %>% sort() %>% unique()
demarkers <- slim.dex[demarkers_idx, ]



dex.up.base <- lapply(all_species, function(sp) {
    sub.dem <- filter(demarkers, species1 == sp) %>%
                group_by(cluster, gene) %>%
                summarize(nhits = n(), mfc = min(ratio_fc)) %>%
                filter(nhits == 3) %>%
                ungroup() %>%
                group_by(cluster) %>%
                arrange(desc(mfc)) %>%
                ungroup() %>%
                mutate(species = sp)
    dem <- split(sub.dem, sub.dem$cluster) %>%
            .[intersect(homo.order, unique(sub.dem$cluster))] %>%
            do.call(rbind, .) %>%
            filter(!grepl("-", gene)) %>%
            filter(!grepl("^AC[0-9]", gene)) %>%
            .$gene %>%
            unique()
    print(sp)
    dem
    }) %>% 
        setNames(., all_species)



enr.list <- lapply(names(dex.up.base), function(sp) {
    cls <- paste0(rep(sp, each = length(homo.order)), "|", homo.order)
    exp_genes <- rownames(exp_ratio)[apply(exp_ratio[, cls, drop = FALSE], 1, max) >= 0.2]

    ## Further filter exp genes
    fgset <- lapply(glist, function(x) intersect(x, exp_genes))
    fgset <- fgset[sapply(fgset, length) >= 8]

    up.genes <- dex.up.base[[sp]]
    henr <- gset_enrich(universe = exp_genes, interesting_genes = up.genes, gset = fgset, min_nodeSize = 6)
    return(henr)
    }) %>%
        setNames(., names(dex.up.base))


gfs <- lapply(enr.list, function(x) x$GFname[x$FDR <= 0.1])

gf.dem <- lapply(names(gfs), function(sp) {
    topgfs <- gfs[[sp]]
    topgenes <- lapply(glist[topgfs], function(y) intersect(dex.up.base[[sp]], y)) %>% unlist() %>%
    unique()
    topgenes
    }) %>%
        setNames(., names(gfs))


save(enr.list, gfs, gf.dem, file = paste0(inputdir, "HGNC_enrich_data.Rdata"))



##-------------------------------------------------------------------
## Barplot showing the significance
all_gfs <- unlist(gfs, use.names = FALSE) %>% unique()
pdata <- lapply(all_species, function(sp) {
    df <- data.frame(row.names = all_gfs,
            gset = all_gfs,
            species = sp,
            FDR = 1,
            stringsAsFactors = FALSE)
    henr <- enr.list[[sp]]
    sh_gf <- intersect(rownames(df), rownames(henr))
    df[sh_gf, "FDR"] <- henr[sh_gf, "FDR"]
    df
    }) %>%
        do.call(rbind, .) %>%
        mutate(logFDR = -log10(FDR))

pdata$gset <- gsub("ADAM metallopeptidases with thrombospondin type 1 motif", "ADAM metallopeptidases", pdata$gset)
gf_ord <- gsub("ADAM metallopeptidases with thrombospondin type 1 motif", "ADAM metallopeptidases", all_gfs)
sp_cols <- setNames(c("#FF420E","#4CB5F5","#89DA59","#FFBB00"), c("Human", "Chimpanzee", "Rhesus", "Marmoset"))
p <- ggplot(pdata, aes_string(x = "gset", y = "logFDR", fill = "species")) +
        geom_bar(stat = "identity", size = 0.2, position = position_dodge()) +
        scale_fill_manual(values = sp_cols) +
        theme_cowplot() + 
        RotatedAxis() +
        scale_x_discrete(limits = gf_ord) +
        geom_hline(yintercept = 1, size = 0.4, linetype = "dashed") +
        theme(axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.text.x = element_text(size = rel(0.6)), axis.text.y = element_text(size = rel(1)), axis.title = element_blank(), legend.position = "none")


pdf(paste0(outputdir, "InN_DEM_enriched_pathways_allspecies.v4.pdf"), width = 5, height = 5)
print(p)
dev.off()




##-------------------------------------------------------------------
## Gene expression
source("../MF7_contact/pie.fun.R")
library(tidyr) ## gather

pie.res <- readRDS(file = paste0(inputdir, "Pie.dot.plot.fig1cluster.InN.rds"))
homo.order <- readRDS(file = paste0(dataDir, "cluster.color.rds")) %>% rownames() %>%
                .[c(41:42, 44:86)]


pie.res$radius[pie.res$radius <= 0.025] <- 0
pie.res$radius <- pie.res$radius*0.4

sp <- "Human"
genes <- gf.dem[[sp]]
sub.pie <- pie.res %>%
            subset(gene %in% genes)


p <- PlotScatterPie2_custominterval(pie.data = sub.pie, group.col = "cluster", feature.col = "gene", r.col = "radius", split.order = c("Human", "Chimpanzee", "Rhesus", "Marmoset", "Shared"), group.order = homo.order, feature.order = genes, scale.expression = FALSE, rsf = 1/0.4, x_scale = 0.5, y_scale = 0.75)
pdf(paste0(outputdir, "Pie.dot.DEM.genefamily.", sp, "..pdf"), width = 18, height = ceiling(0.3 * length(genes)))
print(p);
dev.off()       

























