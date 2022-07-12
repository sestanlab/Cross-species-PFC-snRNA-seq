## Plot MF1
source("../scripts/pfc.fun.R")
library(circlize)


##--------------------------------------------------------------
## Update cluster names
updata_name <- function(x) {
    xx <- gsub("L5 BCL11B SYT2 POU3F1", "L5 FEZF2 BCL11B POU3F1", x) %>%
            gsub("InN SST TH GABRQ", "InN SST HGF GABRQ", .) %>%
            gsub("InN LHX6 TH STON2", "InN LHX6 HGF STON2", .) %>%
            gsub("iAstro", "Astro", .) %>%
            gsub("fAstro", "Astro", .) %>%
            gsub("pAstro", "Astro", .) %>%
            gsub("rAstro", "Astro", .) %>%
            gsub("cOPC PDGFRA TOP2A", "OPC PDGFRA PCDH15", .) %>%
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
            gsub("VLMCA8", "ABCA8",.) %>%
            gsub("L6b", "L6B",.)
    return(xx)
}



## Average expression & expression ratio
avg_res <- readRDS(file = paste0("../data/Expr_avg_ratio_by_all.rds"))[["fig1cluster"]]
avgs <- avg_res$avg
avgs <- avgs[, !grepl("cOPC", colnames(avgs))] ##cOPC was merged with OPC

exp_ratio <- avg_res$ratio
exp_ratio <- exp_ratio[, !grepl("cOPC", colnames(exp_ratio))]

colnames(avgs) <- updata_name(x = colnames(avgs))
colnames(exp_ratio) <- updata_name(x = colnames(exp_ratio))




## Meta data
meta <- readRDS(file = "../data/PFC.filtered.meta.05082022.rds") %>%
			mutate(fig1cluster = as.character(fig1cluster)) %>%
			filter(fig1cluster != "OPC PDGFRA TOP2A")




## Cluster order and dendrogram
dend <- readRDS(file = paste0("./load_files/MF1_plot_dend_v3.rds"))
allcls <- rev(labels(dend))
allmres <- c("L2-3 IT", "L6 IT-2", "L6 IT-1", "L3-5 IT-2", "L3-5 IT-1", "L3-5 IT-3", "L5 ET", "L5-6 NP", "L6 CT", "L6B", 
	"LAMP5 LHX6", "LAMP5 RELN", "VIP", "ADARB2 KCNG1", "SST NPY", "SST", "SST HGF", "PVALB", "PVALB ChC", 
	"Astro", "OPC", "Oligo", "Micro", "Immune", "Endo", "RB", "PC", "SMC", "VLMC")
mres_data <- meta %>%
            group_by(fig1cluster) %>%
            summarize(mres = unique(mres)) %>%
            mutate(fig1cluster = factor(as.character(fig1cluster), levels = allcls))
mres_sum <- mres_data[match(allcls, mres_data$fig1cluster), ] %>%
            group_by(mres) %>%
            summarize(ncls = n()) %>%
            ungroup()
mres_sum <- mres_sum[match(allmres, mres_sum$mres), ]
mres_cusum <- cumsum(mres_sum$ncls)




group_colors <- c("#C51B7D", "#4D9221", "#984ea3", "#999999")##"#08519c"
seg_idx <- c(6,11,25,26,30,34,40,47,55,64,74,76,84,86, 90,92,96,99,103,108,111)
seq_lwd <- rep(0.35, length(seg_idx))
seq_lwd[c(7, 14, 18)] <- 5




## Highligted clusters
hl_half_inter <- 0.05
hl_idx <- c(1, 43, 89, 97, 98)
hl_data <- data.frame(xl = hl_idx + hl_half_inter,
                    xr = hl_idx + (1 - hl_half_inter),
                    yb = 0, 
                    yt = 1, 
                    color = add_transparency(rep(group_colors, c(1, 1, 3, 0)), 0.75),
                    stringsAsFactors = FALSE
                    )
print(hl_data)




##-------------------------------------------------------------------
## Build the plots
pdf(paste0(outputdir, "MF1_circle_v8.pdf"), width = 22, height = 22, useDingbats = FALSE)
circos.par(start.degree = 270, cell.padding = c(0, 0, 0, 0), gap.degree=7.5, circle.margin = 0.35)
circos.initialize(xlim = c(0, length(allcls)), factors = "fasdf")




##------------------------------------------------
## Cluster labels
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.05, panel.fun = function(x, y){
    circos.rect(xleft = hl_data$xl, ybottom = hl_data$yb, 
        xright = hl_data$xr, ytop = 9.5, border=NA, col=hl_data$color) 


    text.cex <- rep(1.15, length(allcls))
    text.cex[c(hl_idx + 1)] <- 1.8
    circos.text(x = 0:(length(allcls)-1)+0.5, 
        y = rep(0,length(allcls)), font = 2, 
        labels = allcls, niceFacing = TRUE, facing="reverse.clockwise", 
        col= rep(group_colors, c(40,46,14,15)), adj=c(1, 0.5), cex=text.cex)
    ll <- rep(3.5,length(seg_idx))
    ll[c(7, 12, 16)] <- 8
})





##---------------------------------------------------------
## Rectangle: label t-type subclass
hwb <- 0.35 ## half width of the bar chart
TxData <- data.frame(xl = c(0, mres_cusum[-length(mres_cusum)]) + 0.5 - hwb, 
                    yb = rep(0, length(mres_cusum)),
                    xr = mres_cusum - 0.5 + hwb, 
                    yt = rep(1, length(mres_cusum)), 
                    color = rep(group_colors, c(10, 9, 4, 6)), 
                    stringsAsFactors = FALSE)
text_col <- rep("white", nrow(TxData))
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.05,panel.fun = function(x, y){
    circos.rect(xleft = TxData$xl, ybottom = TxData$yb, 
        xright = TxData$xr, ytop = TxData$yt, border=NA, col=TxData$color)  
    circos.text((TxData$xl + TxData$xr)/2,
                rep(0.5, length(cumsum)),
        allmres, cex = 1.3, 
        col=text_col, facing="bending.inside", niceFacing=TRUE)
})



##---------------------------------------------------------
## Rectangle: label major cluster groups
TxData <- data.frame(xl = c(0,40,86,99) + 0.5 - hwb, 
                    yb = rep(0, 4),
                    xr = c(40,86,99,114) - 0.5 + hwb, 
                    yt = rep(1, 4), 
                    color = group_colors, 
                    stringsAsFactors = FALSE)
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.05,panel.fun = function(x, y){
    circos.rect(xleft = TxData$xl, ybottom = TxData$yb, 
        xright = TxData$xr, ytop = TxData$yt, border=NA, col=TxData$color)  
    circos.text(c(20, 66, 92),
                rep(0.5,3),
        c("Glutamatergic excitatory subtypes","GABAergic inhibitory subtypes","Glial subtypes"), cex = 1.8, 
        col="white", facing="bending.inside", niceFacing=TRUE)
    circos.text(106.5,
                0.5,
        c("Non-neural subtypes"), cex = 1.8, 
        col="white", facing="bending.inside", niceFacing=TRUE)
})



##---------------------------------------------------------
## Stacked Barplot: species ratio
RatioRes <- meta %>%
                group_by(species, fig1cluster) %>%
                summarize(ncells = n()) %>%
                ungroup() %>%
                group_by(fig1cluster) %>%
                mutate(clsratio = ncells/sum(ncells)) %>%
                ungroup() 
hwb <- 0.35 ## half width of the bar chart
allsp <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
nsp <- length(allsp)
spcols <- c("#FF420E","#4CB5F5","#89DA59","#FFBB00")
RatioData <- lapply(1:length(allcls), function(x){
    subratio <- rep(0, length(allsp)) %>% setNames(., rev(allsp))
    new_ratio <- RatioRes %>% subset(fig1cluster == allcls[x]) %>% column_to_rownames("species")
    subratio[rownames(new_ratio)] <- new_ratio$clsratio

    df <- data.frame(xl = rep(x - 0.5 - hwb, nsp), 
                    yb = c(0, cumsum(subratio[1:3])), 
                    xr = rep(x - 0.5 + hwb, nsp),
                    yt = cumsum(subratio), 
                    color = rev(spcols), 
                    cluster = allcls[x], stringsAsFactors = FALSE)
    return(df)
    }) %>% do.call(rbind, .)


circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1, panel.fun = function(x, y){
    circos.rect(xleft = hl_data$xl, ybottom = hl_data$yb - 0.2, 
        xright = hl_data$xr, ytop = 1.2, border=NA, col=hl_data$color) 

    circos.rect(xleft = RatioData$xl, ybottom = RatioData$yb, 
        xright = RatioData$xr, ytop = RatioData$yt, border=NA, col=RatioData$color)
    circos.segments(seg_idx, rep(0,length(seg_idx)), seg_idx, rep(1,length(seg_idx)), lty=2, lwd=seq_lwd)
})





##------------------------------------------------
#DotPlot
mars <- list(PCDH8 = c(1, 6), OPRK1 = c(7, 11), RORB = c(12, 25), 
			FEZF2 = c(26, 26), HTR2C = c(27, 30), SYT6 = c(31, 34), CTGF = c(35, 40), 
			LAMP5 = c(41, 47), VIP = c(48, 55), KCNG1 = c(56, 64),
			SST = c(65, 74), HGF = c(75, 76), PVALB = c(77, 84), UNC5B = c(85, 86),
			AQP4 = c(87, 90), PDGFRA = c(91, 92), MOG = c(93, 96), P2RY12 = c(97, 99), PTPRC = c(97, 103),
			CLDN5 = c(104, 108), ACTA2 = c(109, 111), CEMIP = c(112, 114))
markers <- names(mars)
DotAvg_scale <- avgs[markers, allcls] %>%
                as.matrix() %>% scale() %>% as.matrix() %>% MinMax(., min = -2.5, max = 2.5)
DotAvg <- DotAvg_scale %>%
                reshape2::melt() %>%
                setNames(., c("gene", "cluster", "ExpAvg")) %>%
                mutate(gene = as.character(gene), cluster = as.character(cluster))
DotRatio <- exp_ratio[markers, allcls] %>%
                as.matrix() %>%
                reshape2::melt() %>%
                setNames(., c("gene", "cluster", "ExpRatio")) %>%
                mutate(gene = as.character(gene), cluster = as.character(cluster))
expcols <- colorRampPalette(c("lightgray", "#FF420E"))(100)
DotData <- full_join(DotAvg, DotRatio, by = c("gene", "cluster")) %>%
                group_by(gene) %>%
                mutate(color = expcols[as.numeric(cut(ExpAvg, breaks = 100))]) %>%
                ungroup() %>%
                mutate(cluster = factor(cluster, levels = allcls)) %>%
                mutate(xaxis = as.numeric(cluster) - 0.5) %>%
                mutate(gene = factor(gene, levels = markers)) %>%
                mutate(yaxis = length(markers) - as.numeric(gene) + 0.5)
SegData <- data.frame(xleft = sapply(mars, function(xx) xx[1]), 
					xright = sapply(mars, function(xx) xx[2]),
					gene = factor(names(mars), levels = names(mars)),
					stringsAsFactors = FALSE) %>%
					mutate(yaxis = length(mars) - as.numeric(gene) + 0.5) 




circos.track(ylim = c(0, length(markers)), bg.border = "black", track.height = 0.25, panel.fun = function(x, y){

    ## Plot the dots
    hg <- length(markers)
    dotSF <- 1.8
    circos.segments(SegData$xleft - 1, SegData$yaxis, 
    				SegData$xright, SegData$yaxis, lty=1, lwd=0.03)
    circos.points(DotData$xaxis, DotData$yaxis, col = DotData$color, pch=16, cex = DotData$ExpRatio * dotSF)
    start_idx <- c(0, seg_idx)
    end_idx <- c(seg_idx, length(allcls))
    ng <- length(markers)
    text_h <- ifelse(1:ng > 13, ng + 12 - 1.5*(1:ng), ng - 2.5 - 1:ng)
    text_h[5] <- 16
    text_h[c(20:22)] <- c(6, 4.5, 2.5)


    circos.rect(xleft = hl_data$xl, ybottom = hl_data$yb, 
        xright = hl_data$xr, ytop = ng, border=NA, col=hl_data$color) 


    circos.text(x = (start_idx + end_idx)/2, 
        y = text_h, font = 3, cex = 1.3, 
        labels = markers, col="black", facing="bending.inside", niceFacing = TRUE)  

    circos.segments(seg_idx, rep(0,length(seg_idx)), seg_idx, rep(length(markers),length(seg_idx)), lty=2, lwd=seq_lwd)


    ## Legend Dots
    hdw <- 1.3
    circos.text(x = 14, y = ng/2, labels = "% expressed", col="black", facing="bending.inside", niceFacing=TRUE)
    x_loc <- c(14 - 1.5*hdw, 14 - 0.5*hdw, 14 + 0.5*hdw, 14 + 1.5*hdw)
    circos.points(x = x_loc,
            y = rep(0.35 * ng, 4),
            pch=16, col="black", cex= c(0.1, 0.2, 0.5, 1) * dotSF)
    circos.text(x = x_loc,
            y = rep(0.25 * ng, 4), 
            labels = as.character(c(10, 20, 50, 100)), col="black", facing="bending.inside", niceFacing=TRUE)


    ## Legend Heatmap
    hh <- 0.035 ##half height of the heatmap legend
    circos.text(x = 21, y = ng/2, labels = "Scaled mean expr", col="black", facing="bending.inside", niceFacing=TRUE)
    x_loc <- seq(19, 23, by = 0.04)
    circos.rect(xleft = x_loc[1:length(expcols)],
                ybottom = rep((0.38 - hh) * ng, length(expcols)),
                xright = x_loc[c(1:length(expcols)) + 1], 
                ytop = rep((0.38 + hh) * ng, length(expcols)), 
                border = NA, col = expcols)
    circos.text(x = 21 + c(0 - 1.5*hdw, 0 + 1.5*hdw),
            y = rep(0.25 * ng, 2), 
            labels = as.character(c(min(DotAvg_scale), max(DotAvg_scale)) %>% round(., digits = 2)), col="black", facing="bending.inside", niceFacing=TRUE)
})




##---------------------------------------------------------
## Barplot: Cluster Ratio
hwb <- 0.4
CratioData <- data.frame(xl = 1:length(allcls) - 0.5 - hwb, 
                    yb = rep(0, length(allcls)),
                    xr = 1:length(allcls) - 0.5 + hwb, 
                    yt = as.numeric(sqrt(table(meta$fig1cluster)[allcls]/nrow(meta))), 
                    color = rep(group_colors,c(40,46,13,15)), 
                    stringsAsFactors = FALSE)
circos.track(ylim = c(0, max(CratioData$yt)), bg.border = NA, track.height = 0.08,panel.fun = function(x, y){
    circos.rect(xleft = CratioData$xl, ybottom = CratioData$yb, 
        xright = CratioData$xr, ytop = CratioData$yt, border=NA, col=CratioData$color)
    yy<-seq(0,max(CratioData$yt),length.out=5)
    yy_lab <- paste0(sapply(yy, function(mm) round((mm^2) * 100, digits = 1)), "%")
    yy_lab[1] <- 0
    circos.yaxis(at=yy, labels=yy_lab, labels.cex=0.56)
    circos.segments(seg_idx, rep(0,length(seg_idx)), seg_idx, rep(max(CratioData$yt),length(seg_idx)), lty=2, lwd=seq_lwd)
})




##---------------------------------------------------------
## Add dendrogram
h <-attr(dend, "height")
dend <- dendextend::color_branches(dend, k = 4, col = rev(group_colors))
circos.track(ylim = c(0, h), bg.border = NA, track.height = 0.075,panel.fun = function(x, y){
    circos.dendrogram(rev(dend))
})


circos.clear()

dev.off()











































