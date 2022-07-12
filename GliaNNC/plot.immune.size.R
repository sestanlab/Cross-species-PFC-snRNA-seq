library(dplyr)
library(ggplot2)


meta <- readRDS("../data/PFC.filtered.meta.05082022.rds") %>%
			filter(mres == 'Immune') %>%
			mutate(fig1cluster = as.character(fig1cluster))


pdata <- meta %>%
			group_by(species, fig1cluster) %>%
			summarize(ncells = n()) %>%
			mutate(species = factor(species, levels = c("Human", "Chimpanzee", "Rhesus", "Marmoset")))
pdata$cluster <- strsplit(pdata$fig1cluster, " ", fixed= TRUE) %>%
			sapply(., function(x) x[1])


p <- ggplot(pdata, aes(x = cluster, y = ncells, fill = species, label = ncells)) +
			geom_bar(stat = "identity", size = 0.2) +
			geom_text(aes(y = ncells + 50), size = 4, angle = 45, nudge_x = 0, nudge_y = 0) +
			theme_classic() +
			lemon::coord_capped_cart(bottom='both', left='both') +
			facet_wrap(vars(species), nrow = 4, ncol = 1, strip.position = "left", scales = "fixed") +
			scale_y_continuous(limits = c(0, 900), breaks = c(0, 300, 600, 900)) +
			theme(axis.line = element_line(size = 0.2),
				axis.ticks = element_line(size = 0.2), 
				axis.title = element_blank(),
				strip.placement = 'outside', 
				strip.background = element_blank(),
				strip.text.y.left = element_text(size = 10),
				axis.text.y = element_text(size = rel(0.9), angle = 90, hjust = 0.5, vjust = 0.5),
				axis.text.x = element_text(size = rel(0.9), angle = 45, hjust = 1))
pdf("./report/Immune_cell_size.pdf", 4, 4)
print(p)
dev.off()


