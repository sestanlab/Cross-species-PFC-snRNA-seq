source("../scripts/pfc.fun.R")


## Fast reading count matrix: assume the first row is the gene names
FastReadCounts <- function(file, showProgress = TRUE, header = FALSE, sep = "\t"){ 
    expmat <- data.table::fread(file, showProgress = showProgress, header = header, sep = sep)

    ## Save the gene names to another object and remove that column
    gene_names <- as.character(expmat[, 1][[1]])
    if (sum(duplicated(gene_names)) >= 1){
        stop("Gene names are duplicated, please check")
    }


    gcol <- colnames(expmat)[1]
    expmat[,(gcol):=NULL]


    ## Convert to Sparse marker
    to_sparse <- function(d_table){
      
      i_list <- lapply(d_table, function(x) which(x != 0))
      counts <- unlist(lapply(i_list, length), use.names = F)

      sparseMatrix(
        i = unlist(i_list, use.names = F),
        j = rep(1:ncol(d_table), counts),
        x = unlist(lapply(d_table, function(x) x[x != 0]), use.names = F),
        dims = dim(d_table),
        dimnames = list(NULL, names(d_table)))
    }

    X <- to_sparse(expmat)
    rownames(X) <- gene_names
    return(X)
}



mat <- Read10X(data.dir = "./Bhaduri_etal_2021", gene.column = 2)


dir <- "./Bhaduri_etal_2021/"
meta <- read.table(file = paste0(dir, "meta.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
		column_to_rownames("Cell.Name")
save(mat, meta, file = paste0(dir, "Bhaduri_etal_2021_RNA.Rdata"))




##---------------------------------------------------------------------------------------------------
## Do the violin plots
dir <- "./Trevino_etal_2021/"
load(file = paste0(dir, "Bhaduri_etal_2021_RNA.Rdata"))
seu <- CreateSeuratObject(mat, meta.data = NULL, assay = "RNA", min.cells = 0, min.features = 0, names.field = 1, names.delim = "_");
seu <- NormalizeData(seu, normalization.method = "LogNormalize")


fmeta <- meta[meta$CellType != "Outlier", ]
fmeta$cluster <- gsub("Dividing", "CycProg", fmeta$CellType) %>%
					gsub("^RG$", "RGC", .) %>%
					gsub("^Neuron$", "ExN", .) %>%
					gsub("^Interneuron$", "InN", .) %>%
					gsub("Microglia", "Micro", .) %>%
					gsub("Vascular", "Vas", .)

features <- c("TREM2", "FOXP2")
cls_ord <- c("CycProg", "RGC", "IPC", "ExN", "InN", "CR", "OPC", "Micro", "Vas")
pdata <- fmeta[, "cluster", drop = FALSE] %>%
            cbind(., t(seu$RNA@data[features, rownames(fmeta)])) %>%
			mutate(cluster = factor(as.character(cluster), levels = cls_ord))

plist <- lapply(features, function(gg) { 
    p <- ggplot(pdata, aes_string(x = "cluster", y = gg, fill = "cluster")) + 
            geom_violin(scale = "width", size = 0.1, adjust = 2,trim =TRUE) + 
            theme_classic() + 
            labs(y = gg) +
            theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(0.9), angle = 0, hjust = 1, vjust = 0.5), axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = rel(0.8)), axis.line.y = element_line(size = 0.2), axis.ticks.y = element_line(size = 0.2), legend.position = "none", plot.margin = unit(c(-0.03, 0, -0.03, 0), "in"))
    p
    })

plist[[length(features)]] <- plist[[length(features)]] +
            theme(axis.text.x=element_text(size = rel(0.9), angle = 45, hjust = 1, vjust = 1), axis.line.x = element_line(size = 0.2), axis.ticks.x = element_line(size = 0.2)) 
pdf(paste0("./report/", "Bhaduri_etal_FOXP2_expr_violin.pdf"), width = 4, height = 2.3)
patchwork::wrap_plots(plist, nrow = length(features), ncol = 1)
dev.off()



