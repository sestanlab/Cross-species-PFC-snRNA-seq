source("../scripts/pfc.fun.R")
source("./calc.heter.fun.R")
library(tibble)


cls_ord <- c("L2_3_IT", "L3_5_IT_1","L3_5_IT_2","L3_5_IT_3", "L5_PT", "L5_6_NP", "L6_CT", "L6_IT_1", "L6_IT_2", "L6b",
        "LAMP5_LHX6", "LAMP5_RELN", "VIP", "ADARB2", "SST", "TH", "PVALB","PVALB_CHC",
        "Astro", "OPC","Oligo", "Micro","immune", "blood","Endo","PC","SMC","VLMC") ##"SST_NPY", 
all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
pdata <- lapply(cls_ord, function(cls) {
	print(cls)
	if (file.exists(paste0(inputdir, "Bootstrap_PCA_v1/Bootstrap_PCA-MeanDist_", cls, "_res.rds"))){
	df <- readRDS(file = paste0(inputdir, "Bootstrap_PCA_v1/Bootstrap_PCA-MeanDist_", cls, "_res.rds")) %>%
			filter(type == "bootstrap") %>%
			group_by(idx, species) %>%
			summarize(emp_mean = mean(meanDis)) %>%
			mutate(cluster = cls)
	return(df)	
	} else {
		return(NULL)
	}
	
	}) %>%
		do.call(rbind, .)



pdata <- pdata %>%
			mutate(species = factor(as.character(species), levels = all_sps))
pdata$group <- "NNC"
pdata$group[pdata$cluster %in% c("L2_3_IT", "L3_5_IT_1","L3_5_IT_2","L3_5_IT_3", "L5_PT", "L5_6_NP", "L6_CT", "L6_IT_1", "L6_IT_2", "L6b")] <- "ExN"
pdata$group[pdata$cluster %in% c("LAMP5_LHX6", "LAMP5_RELN", "VIP", "ADARB2", "SST_NPY", "SST", "TH", "PVALB","PVALB_CHC")] <- "InN"


cls_transfer <- c("L2-3 IT", "L3-5 IT-1", "L3-5 IT-2", "L3-5 IT-3", "L5 ET", "L5-6 NP", "L6 CT", "L6 IT-1", "L6 IT-2", "L6B", 
	"LAMP5 LHX6", "LAMP5 RELN", "VIP", "ADARB2 KCNG1", "SST", "SST HGF", "PVALB", "PVALB ChC", 
	"Astro","OPC","Oligo","Micro","Immune","RB","Endo","PC","SMC","VLMC") %>%
     setNames(., c("L2_3_IT", "L3_5_IT_1", "L3_5_IT_2", "L3_5_IT_3", "L5_PT", "L5_6_NP", "L6_CT", "L6_IT_1", "L6_IT_2","L6b", 
    "LAMP5_LHX6", "LAMP5_RELN", "VIP", "ADARB2", "SST", "TH", "PVALB","PVALB_CHC", 
    "Astro","OPC","Oligo","Micro","immune","blood","Endo","PC","SMC","VLMC"))
pdata$cluster2 <- cls_transfer[pdata$cluster]


spcols <- c("#FF420E","#4CB5F5","#89DA59","#FFBB00") %>%
            setNames(., c("Human", "Chimpanzee", "Rhesus", "Marmoset"))
## Which clusters to label "**"
source("./vis.heter.fun.R")
pval_res <- PairWiseTest(data = pdata, group.by = "cluster2", split.by = "species", value.col = "emp_mean", p.adjust.method = "fdr", fc.thre = 0, pseudo.count = 0.01)
sig_df <- pval_res$df %>%
		filter(species2 == "Marmoset") %>%
		group_by(cluster) %>% 
		summarize(nhits = sum(pval_final <= 0.01)) %>%
		filter(nhits == 3) %>%
		ungroup() %>%
		mutate(label = "**", y = ceiling(max(pdata$emp_mean) * 1.06))

p1 <- ggplot(pdata, aes(x = cluster2, y = emp_mean)) +
				geom_violin(aes(fill = species), scale = "width", size = 0.05, width = 0.8, position=position_dodge(.7)) +
            	geom_boxplot(aes(group = interaction(species, cluster2)), fill = "white", size = 0.05, width=.3, position=position_dodge(.7), color = "black", outlier.shape = NA) +
            	theme_bw() +
            	RotatedAxis() +
             	scale_fill_manual(values = spcols) +
             	scale_color_manual(values = spcols) +
            	scale_x_discrete(limits = setNames(cls_transfer, NULL)) +
            	annotate("text", x = sig_df$cluster, y = sig_df$y, label = sig_df$label) +
            	theme(legend.position = "none", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(size = 0.2),
              	plot.title = element_blank())
pdf(paste0(outputdir, "Bootstrap_4sps_PCA-MeanDist.v1.pdf"), width = 7, height = 3)
print(p1)
dev.off()











