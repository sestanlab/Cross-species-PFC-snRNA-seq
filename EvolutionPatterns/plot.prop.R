source("../scripts/pfc.fun.R")
library(lemon)


meta <- readRDS(file = paste0(dataDir, "PFC.filtered.meta.05082022.rds"))


group.by <- "fig1cluster"
split.by <- "species"
ind.by <- "samplename"
sp_colors <- c("#FF420E", "#4CB5F5", "#89DA59", "#FFBB00") %>% setNames(., c("Human", "Chimpanzee", "Rhesus", "Marmoset"))
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
dend <- readRDS(file = paste0("../MF1_overview/load_files/MF1_plot_dend_v3.rds"))
cls_ord <- rev(labels(dend))



meta$group2 <- meta$group
meta$group2[meta$mres %in% "Astro"] <- "Astro"
meta$group2[meta$mres %in% c("Immune", "Micro")] <- "Immune"
meta$group2[meta$mres %in% c("OPC", "Oligo")] <- "Oligo"
meta$group2[meta$mres %in% c("RB", "Endo", "SMC", "PC", "VLMC")] <- "Vas"
all_gps <- c("ExN", "InN", "Astro", "Oligo", "Immune", "Vas")
meta$fig1cluster <- gsub("OPC PDGFRA TOP2A", "OPC PDGFRA PCDH15", as.character(meta$fig1cluster))




data <- lapply(all_gps, function(gp) {
			subdata <- meta %>%
				filter(group2 == gp) %>%
                group_by(!!sym(ind.by), !!sym(group.by)) %>%
                summarise(size = n(), !!split.by := unique(!!sym(split.by))) %>%
                ungroup() %>%
                group_by(!!sym(ind.by)) %>%
                mutate(ratio = size * 100/sum(size)) %>%
                ungroup() %>%
                group_by(!!sym(split.by), !!sym(group.by)) %>%
                summarize(mRatio = mean(ratio), rse = sd(ratio)/sqrt(n())) %>%
                ungroup() %>%
                mutate(!!split.by := factor(!!sym(split.by), levels = rev(names(sp_colors))))
            return(subdata)
            }) %>%
			do.call(rbind, .) %>%
            mutate(!!group.by := factor(!!sym(group.by), levels = rev(cls_ord)))



library(ggbreak) 

p <- ggplot(data, aes_string(y = group.by, x = "mRatio", fill = split.by)) +
                  geom_bar(size = 0.01, width=.8, stat="identity", color="black", position=position_dodge(.8)) +
                  geom_errorbar(aes(xmin=mRatio-rse, xmax=mRatio+rse), width=.23, position=position_dodge(.7), size = 0.2) +
                  theme_bw() +
                  RotatedAxis() +
                  scale_fill_manual(values = sp_colors) +
                  scale_x_break(c(10, 20), scales = 0.5) +
                  theme(legend.position = "none", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(size = 0.27), axis.text.y = element_text(size = rel(0.7)))
pdf(paste0(outputdir, "SF.Celltype.proportions_v1.pdf"), width = 5, height = 10)
print(p)
dev.off()















