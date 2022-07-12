## Load Ribosome-binding proteins dataset
source("../scripts/pfc.fun.R")
library(dplyr)


##-------------------------------------------------------------------------------------------
## Transfer ENSEMBLE ID to gene names
gtf <- rtracklayer::import("../mappability/fastagtf/hg38_final.gtf") %>%
					as.data.frame(., check.names = FALSE)


gtf$id <- extract_field(gtf$gene_id, 1, ".") %>% 
			setNames(., NULL)
df <- read.table(file = paste0("./load_files/", "Human_ribosome_binding_proteins.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "", quote = "\"")
df <- df %>%
		filter(gene_id %in% gtf$id) %>%
		mutate(gene_name = gtf$gene_name[match(gene_id, gtf$id)])
idx <- which(df$gene_name != df$gene_symbol)
saveRDS(df, file = paste0("./load_files/RBPgenes.db.rds"))



##-------------------------------------------------------------------------------------------
## Evaluate cell type differences & species differences
df <- readRDS(file = paste0("./load_files/RBPgenes.db.rds"))

inn <- readRDS(file = paste0("../data/", "Final_PFC_HPRC.InN.07102021.rds"))
inn <- subset(inn, mres %in% c("SST", "TH"))

inn@meta.data$group <- ifelse(as.character(inn$fig1cluster) == "InN SST TH GABRQ", "TH-GABRQ", "SST-bg")
inn@meta.data$group[inn@meta.data$fig1cluster == "InN LHX6 TH STON2"] <- "TH-STON2"



## Gene enriched in SST HGF GABRQ versus other SST cells
inn@meta.data$cmp_cls <- paste0(inn@meta.data$species, "|", inn@meta.data$group)
Idents(inn) <- "cmp_cls"
res1 <- FindMarkers(inn, ident.1 = "Human|TH-GABRQ", ident.2 = "Human|SST-bg", only.pos = TRUE, min.pct = 0.15, features = df$gene_name, logfc.threshold = 0.1) %>%
			mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
			filter(p_val_adj <= 0.01)
res2 <- FindMarkers(inn, ident.1 = "Human|TH-GABRQ", ident.2 = "Rhesus|TH-GABRQ", only.pos = TRUE, min.pct = 0.15, features = df$gene_name, logfc.threshold = 0.1) %>%
			mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
			filter(p_val_adj <= 0.01)
res3 <- FindMarkers(inn, ident.1 = "Human|TH-GABRQ", ident.2 = "Marmoset|TH-GABRQ", only.pos = TRUE, min.pct = 0.15, features = df$gene_name, logfc.threshold = 0.1) %>%
			mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
			filter(p_val_adj <= 0.01)
hi_genes <- intersect(rownames(res1), rownames(res2)) %>% intersect(., rownames(res3))
hi_genes


## Get average gene expression
Idents(inn) <- "cmp_cls"
avgs <- log(AverageExpression(inn, assay = "RNA", features = df$gene_name)$RNA + 1)


## Do the plot
pdata <- avgs[, c("Human|TH-GABRQ", "Human|SST-bg", "Rhesus|TH-GABRQ", "Marmoset|TH-GABRQ"), drop = FALSE] %>%
			as.data.frame() %>%
			setNames(., c("H_GABRQ", "H_bg", "R_GABRQ", "C_GABRQ")) %>%
			rownames_to_column("gene") %>%
			mutate(color = ifelse(gene %in% hi_genes, "red", "lightgrey")) %>%
			tidyr::gather(., "type", "x_avg", c("H_bg", "R_GABRQ", "C_GABRQ")) %>%
			mutate(type = factor(type, levels = c("H_bg", "R_GABRQ", "C_GABRQ"))) %>%
			mutate(label = ifelse(gene %in% hi_genes, gene, NA)) %>%
			mutate(size = ifelse(gene %in% hi_genes, 2.5, 1.3))
library(ggrastr)
set.seed(42)
p <- ggplot(pdata, aes_string(x = "x_avg", y = "H_GABRQ", color = "color", label = "label")) +
			rasterise(geom_point(aes_string(size = "size"), shape = 16), dpi = 300, scale = 1) +
			geom_text_repel(size = 4) +
			theme_bw() +
			theme(panel.grid.major = element_line(size = 0.3), panel.grid.minor = element_line(size = 0.3)) +
			scale_size_identity() +
			coord_fixed(ratio = 1, xlim = c(0, 2.8), ylim = c(0, 2.8)) + 
			geom_abline(intercept = 0, slope = 1) +
			scale_color_identity() +
			facet_wrap(vars(type), nrow = 1, ncol = 3) +
			theme(strip.background = element_blank(), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.text = element_text(size = rel(0.8)))
pdf(paste0("./report/GABRQ_specific_RBPs.pdf"), width = 10, height = 3.3)
print(p)
dev.off()






















