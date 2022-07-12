## MF2B visualization
args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")
source("./network.fun.v2.R")



## Load Meta data
meta <- readRDS(paste0(dataDir, "PFC.filtered.meta.05082022.rds")) %>%
				mutate(group = gsub("Glia", "NNC", group)) %>%
				mutate(fig1cluster = as.character(fig1cluster)) %>%
				mutate(fig1cluster = gsub("OPC PDGFRA TOP2A", "OPC PDGFRA PCDH15", fig1cluster)) ## merge cOPC	 with OPC
avg_res <- readRDS(file = paste0(dataDir, "Expr_avg_ratio_by_all.rds"))[["fig1cluster"]]




## Calculate node size
cls_size <- meta %>%
                group_by(samplename, fig1cluster) %>%
                summarize(ncells = n(), group = unique(group), species = unique(species)) %>%
                ungroup() %>%
                group_by(samplename, group) %>%
                mutate(percentage = ncells * 100/sum(ncells)) %>%
                ungroup() %>%
				group_by(fig1cluster) %>%
                summarize(msize = mean(percentage)) %>%
                column_to_rownames("fig1cluster") %>%
                as.matrix() %>% .[, 1] %>%
                sqrt()



## Calculate node links
ctp <- args[1]


hvg <- readRDS(file = paste0(inputdir, "HVG_CTP_Marker_", ctp, "_metrics.rds"))$gene
avgs <- avg_res$avg[hvg, ]
colnames(avgs) <- updata_name(colnames(avgs))
cls_use <- meta %>%
			filter(group == ctp) %>%
			.$fig1cluster %>% as.factor() %>% levels()
avgs <- avgs[, cls_use]



## Organize the data together
res <- PrepNet(data = avgs, qt = 0.9)
prop_list <- CalcProps(meta, ctp = ctp)



## Visualization
res <- c(res, list(size = cls_size, props = prop_list))
PlotNet(res = res, file_name = paste0("MF2B_Pie_Net_", ctp), dot_scale = 8)




## Build pseudolist for the legend
set.seed(0)
node_meta <- data.frame(cluster = sample(letters, 6), stringsAsFactors = FALSE)
df <- data.frame(from = node_meta$cluster, to = node_meta$cluster[c(2:6, 1)], weight = 1)
cls_size <- setNames(sqrt(c(0.025, 0.25, 1, 4, 16, 32)), node_meta$cluster)
prop_list <- lapply(names(cls_size), function(x) c(Human = 0.5, Chimpanzee = 0, Rhesus = 0, Marmoset = 0)) %>%
				setNames(., names(cls_size))
res <- list(links = df,
			meta = node_meta, 
			size = cls_size, 
			props = prop_list)
source("./network.fun.v2.R")
PlotNet(res = res, file_name = paste0("MF2B_Pie_Net_legend"), dot_scale = 8)


















