## This script calculate the markers
args <- commandArgs(trailingOnly = TRUE)
source("./pfc.fun.R")
source("./zfun.wilcox.R")


## Set the parameters
group.by <- args[1]
sp <- args[2]
ctp <- args[3]
idx <- as.numeric(args[4])
seu <- readRDS(file = paste0(dataDir, "Final_PFC_HPRC.", ctp, ".rm1680.07102021.rds"))
exp_genes <- readRDS(file = paste0(dataDir, "Expgenes.hprc.rds"))
seu <- seu[exp_genes, seu@meta.data$species == sp]


allcls <- as.character(seu@meta.data[, group.by]) %>% as.factor() %>% levels()
## Subset the dataset to make sure each cluster has comparable number of cells across species
set.seed(0)
sub_cells <- lapply(allcls, function(cls) {
		cells <- colnames(seu)[seu@meta.data[, group.by] == cls]
		if (length(cells) > 5000){
			cells <- sample(cells, 5000)
		}
		cells
		}) %>% 
			unlist()

seu <- seu[, sub_cells]


ident1_list <- lapply(allcls, function(x) x) %>% setNames(., allcls)
ident2_list <- lapply(allcls, function(x) setdiff(allcls, x)) %>% setNames(., allcls)



cls.bin <- cut(1:length(ident1_list), breaks = 4) %>% as.numeric()
cls_use <- which(cls.bin == idx)


sp_mar <- allWilcox(object = seu, group.by = group.by, ident_1 = ident1_list[cls_use], ident_2 = ident2_list[cls_use], pseudocount.use = 0, adjust_method="bonferroni", nCores = 4, max.cells.per.ident = 5000)
saveRDS(sp_mar, file = paste0(dataDir, "rawMar/", "Markers.hprc.", group.by, ".", sp, ".", ctp, ".", as.character(idx), ".rds"))


