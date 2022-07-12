source("../scripts/pfc.fun.R")


dir_use <- "./species_micro/"
all_files <- list.files(dir_use)


pdata <- lapply(all_files, function(ff) {
    sp <- extract_field(ff, 2, "_")
    df <- read.table(file = paste0(dir_use, ff), sep = "\t", stringsAsFactors = FALSE, header = TRUE) %>%
                as.matrix() %>%
                reshape2::melt() %>%
                setNames(., c("gene", "samplename", "exp")) %>%
                mutate(exp = as.numeric(exp)) %>%
                mutate(species = sp)

    ## Get the CPM values
    df <- df %>% 
            group_by(samplename) %>%
            mutate(libsize = sum(exp)) %>%
            ungroup() %>%
            mutate(CPM = exp * 10e6 / libsize) %>%
            mutate(logcpm = log2(CPM + 1))

    df
    }) %>%
        do.call(rbind, .) %>%
        filter(gene %in% c("FOXP2", "Foxp2"))

sp_infor <- c("Human", "Rhesus", "Marmoset", "Sheep", "Mouse", "Rat", "Hamster", "Chicken") %>%
               setNames(., c("human","macaca","marmoset","sheep","mouse","rat", "hamster","chicken"))
pdata$species <- sp_infor[pdata$species]

sp_cols <-  c("#FF420E", "#89DA59", "#FFBB00", "#8dd3c7", "#AE017E", "#9970ab", "#bf812d", "#8c510a") %>% 
        setNames(., c("Human", "Rhesus", "Marmoset", "Sheep", "Mouse", "Rat", "Hamster", "Chicken"))

p <- ggplot() + 
        geom_boxplot(data = pdata, aes_string(x = "species", y = "logcpm", fill = "species"), color = "black", width = 0.5, size = 0.2, outlier.shape = NA) + 
        theme_cowplot() + 
        RotatedAxis() + 
        scale_x_discrete(limits = setNames(sp_infor, NULL)) +
        scale_fill_manual(values = sp_cols) + 
        labs(y = "log2(CPM + 1)") +
        theme(axis.title.x = element_blank(), axis.text.x=element_text(size = rel(0.9)), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), legend.position = "bottom", legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size = rel(0.5)), legend.title = element_blank())
pdf(paste0(outputdir, "Crossspecies.micro.FOXP2.expression.pdf"), width = 6, height = 4)
print(p)
dev.off() 







