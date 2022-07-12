source("../scripts/pfc.fun.R")
source("../MF6_InN/inn.fun.R")


pfc <- readRDS(paste0(dataDir, file = "Final_PFC_HPRC.InN.07102021.rds"))
hmtg <- readRDS(file = paste0(dataDir, "Human_MTG_org_05072021.rds"))


pfc.cls <- c("InN SST TH GABRQ", "InN LHX6 TH STON2")
mtg.cls <- c("Inh L5-6 SST TH", "Inh L5-6 GAD1 GLP1R")


pfc.sub <- subset(pfc, fig1cluster %in% pfc.cls)
pfc.sub@meta.data$fig1cluster <- as.character(pfc.sub@meta.data$fig1cluster)

mtg.sub <- subset(hmtg, cluster %in% mtg.cls)
mtg.sub@meta.data$fig1cluster <- gsub("Inh L5-6 SST TH", "InN SST TH GABRQ", mtg.sub@meta.data$cluster) %>%
									gsub("Inh L5-6 GAD1 GLP1R", "InN LHX6 TH STON2", .)
mtg.sub@meta.data$species <- "HMTG"


seu <- merge(pfc.sub, mtg.sub)



sp_cols <- setNames(c("#FF420E","#4CB5F5","#89DA59","#FFBB00", "#FF420E"), c("Human", "Chimpanzee", "Rhesus", "Marmoset", "HMTG"))
genes <- c("HGF", "SST", "GABRQ", "STON2", "TH", "GCH1", "DDC", "SLC18A2", "DRD2", "SLC6A3", "SLC47A1", "SERPINE1")
p <- PointPlot(seu, features = genes, dot.scale = 5, cols = sp_cols, dot.min = 0.02, group.by = "fig1cluster", split.by = "species", cluster_order = pfc.cls, shape = 16, species_order = c("Human", "Chimpanzee", "Rhesus", "Marmoset", "HMTG"), scale.by = "size") +
				coord_flip() +
		        RotatedAxis() +
		        facet_wrap(vars(species), nrow = 1, ncol = 5) +
		        theme(axis.text.x=element_text(size = rel(0.6)), strip.background = element_blank(), strip.text = element_blank(), panel.spacing = unit(0.1, "in"), legend.position = "bottom", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.text.y = element_text(size = rel(0.8)), axis.title = element_blank())
pdf(paste0(outputdir, "InN_dopamine-pathway_exp_v2.pdf"), width = 4.3, height = 3.5, useDingbats = FALSE)
print(p)
dev.off()











