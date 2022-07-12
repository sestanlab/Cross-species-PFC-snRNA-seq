## Plot the average expression of FLT1-VEGFA signaling in Vas cells
source("../scripts/pfc.fun.R")



## Combinatory melt of avg and ratio matrices
## a function from this script: ../scripts/exp.fun.R
CombnMelt <- function(avg, ratio, stroke = NULL) {
    Mat2Df <- function(mat, name) {
        df <- mat %>% 
                    as.matrix() %>%
                    reshape2::melt() %>%
                    setNames(., c("features.plot", "newcls", name)) %>%
                    mutate(features.plot = as.character(features.plot), newcls = as.character(newcls))
        return(df)
    }


    ## Combine average and expression ratio
    avg.exp <- Mat2Df(mat = avg, name = "avg.exp.scaled")
    pct_exp <- Mat2Df(mat = ratio, name = "pct.exp")
    pre.plot <- dplyr::full_join(avg.exp, pct_exp, by = c("features.plot", "newcls"))


    ## Add stroke data if present
    if (!is.null(stroke)){
        pre.plot <- Mat2Df(mat = ratio, name = "stroke") %>%
                        dplyr::full_join(pre.plot, ., by = c("features.plot", "newcls")) 
    }


    ## Further format the columns
    data.plot <- pre.plot  %>%
                    mutate(species = extract_field(newcls, 1, "|"),  cluster = extract_field(newcls, 2, "|"))%>%
                    select(-newcls)

    return(data.plot)
}



res <- readRDS(file = paste0(dataDir, "Expr_avg_ratio_by_species.rds"))$fig1cluster
sel_cls <- c("cEndo CLDN5 SLC7A5", "vEndo CLDN5 IL1R1", "aEndo CLDN5 DKK2", "aaSMC ACTA2 CYP1B1", "aSMC ACTA2 CNN1", "vSMC ACTA2 CRISPLD2", "PC P2RY14 GRM8")



features <- c("CLDN5", "DKK2", "SLC7A5", "IL1R1", "FLT1", "PGF", "VEGFB", "VEGFA", "GRM8", "CYP1B1", "CNN1", "ACTA2")
avg_use <- res$avg[features, grep(paste(sel_cls, collapse = "|"), colnames(res$avg))]
ratio_use <- res$ratio[features, colnames(avg_use)]
colnames(avg_use) <- colnames(ratio_use) <- extract_field(colnames(avg_use), 1, " ")




all_df <- CombnMelt(avg = avg_scale, ratio = ratio_use, stroke = NULL) 
all_df$pct.exp <- all_df$pct.exp * 100
all_df$pct.exp[all_df$pct.exp < 5] <- NA
all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
sp_cols <- setNames(c("#FF420E","#4CB5F5","#89DA59","#FFBB00"), all_sps)

all_df$exp <- as.numeric(cut(all_df$avg.exp.scaled, 20))
all_df$color <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("lightgrey", color))(20)[value])
        }, color = sp_cols[as.character(all_df$species)], value = all_df$exp)


cls_list <- list(c1 = c("aEndo", "aSMC"), 
				c2 = c("aEndo", "aaSMC"), 
				c3 = c("cEndo", "PC"),
				c4 = c("vEndo", "vSMC"))

plist <- lapply(cls_list, function(cls) {
	avg_scale <- avg_use %>%
					as.matrix() %>%
					t() %>% scale() %>% t() %>%
					MinMax(., min = -2.5, max = 2.5)

	idx1 <- all_df$cluster == cls[1] & all_df$features.plot %in% features[1:5]
	idx2 <- all_df$cluster == cls[2] & all_df$features.plot %in% features[6:12]
	plot_data <- all_df[idx1 | idx2, ]


	p <- ggplot(plot_data, aes(x = features.plot, y = species, color = color, size = pct.exp)) +
			geom_point(shape = 16) +
			scale_radius(range = c(0, 5.5), limits = c(0, 100)) +
			scale_color_identity() +
			theme_classic() +
			coord_fixed() +
			scale_y_discrete(limits = rev(all_sps)) +
			scale_x_discrete(limits = features) +
			theme(legend.position = "right", axis.line = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.y = element_blank(), panel.grid.major = element_line(size = 0.15, color = "#4D4D4D"))
	return(p)
	})

plist[[4]] <- plist[[4]] + 
				theme(axis.text.x = element_text(size = rel(0.5), angle = 45, hjust = 1))


pdf(paste0(outputdir, "FLT1_signal.Vas.pdf"), width = 3, height = 4, useDingbats = FALSE)
patchwork::wrap_plots(plist, ncol = 1, nrow = 4, guides = "collect")
dev.off()























