work_dir <- getwd()
inputdir <- paste0(work_dir, "/load_files/")
outputdir <- paste0(work_dir, "/report/")
dataDir <- "~/project/PFC/data/"


#Load all related packages 
library(ggplot2); 
library(ggrepel); 
library(cowplot); 
library(tibble);
library(dplyr); 
library(Matrix); 
library(parallel); 
library(viridis);
library(batchelor)
library(ComplexHeatmap)
library(circlize)
if (sum(grepl("Seurat", (.packages()))) < 1){
    library(Seurat)#if source twice, the non-developmental version could be loaded and functions may change.
}



map_gene <- function(gene_names, input_genes,ignore_case=TRUE){
    input_genes <- unique(input_genes)
  
    if (sum(grepl("\\|",gene_names))==length(gene_names)){
        if (sum(grepl("\\|",input_genes))==length(input_genes)){
              gene_id <- input_genes[input_genes %in% gene_names]
          }else{
            input_genes <- extract_field(input_genes=input_genes, 2, "|")
              gene_id <- unlist(mclapply(input_genes, function(i) ifelse(length(grep(paste0("\\|",i,"$"),gene_names, ignore.case= ignore_case))==1,gene_names[grep(paste0("\\|",i,"$"),gene_names, ignore.case= ignore_case)],"empty")))
              gene_id <- gene_id[gene_id != "empty"]
          }
    } else if(sum(grepl("\\|",gene_names))==0){
        input_genes <- extract_field(input_genes=input_genes, 2, "|")
          gene_id <- unlist(mclapply(input_genes, function(i) ifelse(length(grep(paste0("^",i,"$"),gene_names, ignore.case= ignore_case))==1,gene_names[grep(paste0("^",i,"$"),gene_names, ignore.case= ignore_case)],"empty")))
          gene_id <- gene_id[gene_id != "empty"]
    } else {
        stop("Inconsistent gene name format")
    }
    return(gene_id)
}



#Get the mito genes based on the inquiry genes
get_genes <- function(input_genes, gene_type = c("mito","ribo", "cc")[1], return_list = FALSE, revised = FALSE){
    gene_use <- list()
    if ("mito" %in% gene_type){
        mito.known <- map_gene(gene_names=input_genes, input_genes=c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB")) #Refseq annotation 103 (Macaque)
        mito.combine <- grep(pattern = "\\|MT-", x = input_genes, value = TRUE, ignore.case=TRUE) #All human mitochondria genes
        mito.single <- grep(pattern = "^MT-", x = input_genes, value = TRUE, ignore.case=TRUE) #All human mitochondria genes
        gene_use[["mito"]] <- unique(c(mito.known, mito.combine, mito.single))
    }

    if ("ribo" %in% gene_type){
        ribo.combine <- c(grep("\\|RPL", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|RPS", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|MRPL", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|MRPS", input_genes, value =TRUE, ignore.case=TRUE))
        ribo.single <- c(grep("^RPL", input_genes, value =TRUE, ignore.case=TRUE), grep("^RPS", input_genes, value =TRUE, ignore.case=TRUE), grep("^MRPL", input_genes, value =TRUE, ignore.case=TRUE), grep("^MRPS", input_genes, value =TRUE, ignore.case=TRUE))
        gene_use[["ribo"]] <- c(ribo.combine, ribo.single)
    }

    if ("cc" %in% gene_type){
        if (!revised){
            regev.s.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'s_genes.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
            regev.g2m.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'g2m_genes.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
        } else {
            regev.s.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'s_genes_revised_04202019.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
            regev.g2m.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'g2m_genes_revised_04202019.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
        }
        
        gene_use[["s"]] <- map_gene(gene_names=input_genes, input_genes=regev.s.genes,ignore_case=TRUE)
        gene_use[["g2m"]] <- map_gene(gene_names=input_genes, input_genes=regev.g2m.genes,ignore_case=TRUE)
    }

    if (return_list){
        return(gene_use)
    } else {
        all_genes <- setNames(unlist(gene_use), NULL)
        return(all_genes)
    }
}




seu_prepare <- function(counts, data = NULL, min.cells = 5, normalization.method = c("LogNormalize", "none","scran")[1], nfeatures = 2500, hvg.method = "vst", assay = "RNA") {
    #Check the normalization.method
    if (!is.null(data)){
        message("input-norm is provided and therefore normalization is not required")
        normalization.method <- "none"
    }

    if (normalization.method == "none"){
        message("normalization.method is none and therefore assume the counts has already been normlized")
    } else if (normalization.method == "LogNormalize"){
        message("normalization.method is LogNormalize and therefore will use seurat LogNormalize to transform the dataset")
    } else {
        stop("Unknown normalization.method")
    }

    inseu <- CreateSeuratObject(counts, meta.data = NULL, assay = "RNA", min.cells = min.cells, min.features = 0, names.field = 1, names.delim = "_"); #at least expressed in 5 cells #only RNA assay is supported

    #quality_gene list
    quality_genes <- get_genes(input_genes = rownames(inseu$RNA@data), gene_type = c("mito","ribo", "cc"), return_list = TRUE, revised = FALSE)
    inseu[["percent.mt"]] <- PercentageFeatureSet(inseu, pattern = NULL, features = quality_genes[["mito"]], col.name = NULL)
    inseu[["percent.ribo"]] <- PercentageFeatureSet(inseu, pattern = NULL, features = quality_genes[["ribo"]], col.name = NULL)


    #normalize the dataset if needed
    if (normalization.method == "none"){
        inseu[[assay]] <- CreateAssayObject(data = as(ifelse_check(is.null(data), inseu[["RNA"]]@data, data), "dgCMatrix"), min.cells = 0, min.features = 0)
        DefaultAssay(inseu) <- assay
    } else if (normalization.method == "LogNormalize"){
        inseu <- NormalizeData(inseu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    }
    
    if (!is.null(hvg.method)){
        inseu <- FindVariableFeatures(inseu, selection.method = hvg.method, nfeatures = nfeatures, verbose = FALSE)
    }
    return(inseu)
}




SelectHVG <- function(hvg.list, nfeatures = 2000) {

    var.features <- unname(obj = unlist(x = hvg.list)) %>%
                        table() %>%
                        sort(., decreasing = TRUE)

    tie.val <- var.features[min(nfeatures, length(x = var.features))]
    features <- names(x = var.features[which(x = var.features > tie.val)])
    if (length(x = features) > 0) {
        feature.ranks <- sapply(X = features, FUN = function(x) {
            ranks <- sapply(X = hvg.list, FUN = function(vf) {
                if (x %in% vf) {
                  return(which(x = x == vf))
                }
                return(NULL)
            })
            median(x = unlist(x = ranks))
        })
        features <- names(x = sort(x = feature.ranks))
    }
    features.tie <- var.features[which(x = var.features == tie.val)]
    tie.ranks <- sapply(X = names(x = features.tie), FUN = function(x) {
        ranks <- sapply(X = hvg.list, FUN = function(vf) {
            if (x %in% vf) {
                return(which(x = x == vf))
            }
            return(NULL)
        })
        median(x = unlist(x = ranks))
    })
    features <- c(features, names(x = head(x = sort(x = tie.ranks), 
        nfeatures - length(x = features))))
    return(features)
}



extract_field <- function(input_genes, field = 1, split_by = "_") {
    split_strings <- strsplit(input_genes, split_by, fixed=TRUE)
    if (is.numeric(field)){
        if (all(field > 0)){
            if (length(field) == 1){
                idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
                output_strings <- sapply(1:length(split_strings), function(x) split_strings[[x]][idx[x]])
            } else {
                idx <- sapply(split_strings, length) == 1
                output_strings <- sapply(1:length(split_strings), function(x) ifelse(idx[x], input_genes[x], paste(split_strings[[x]][field], collapse = split_by)))
            }
        } else {
            if (length(field) == 1) {
                field = as.integer(abs(field))
                idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
                output_strings <- sapply(1:length(split_strings), function(x) rev(split_strings[[x]])[idx[x]])
            } else {
                stop("currently doesnot support field with length >= 2")
            }
        }
    } else if (is.character(field)) {
        if (field == "rm_start"){
            idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
            output_strings <- sapply(1:length(split_strings), function(x) paste(split_strings[[x]][-1],collapse=split_by))
        } else if (field == "rm_end") {
            output_strings <- sapply(1:length(split_strings), function(x) paste(rev(rev(split_strings[[x]])[-1]),collapse=split_by))
        } else {
            stop("Currently only support rm_start, rm_end for the character field")
        }
    }
    output_strings <- output_strings %>% setNames(., input_genes)
    return(output_strings)
}


#Function in Seurat package, useful
MinMax <- function (data, min, max) {
    data2 <- data
    data2[data2 > max] <- max
    data2[data2 < min] <- min
    return(data2)
}


remove_duplicates <- function(x, return_index = FALSE){
    dp_index <- !duplicated(x) & !rev(duplicated(rev(x)))
    if(return_index){
        return(dp_index)
    } else {
        return(x[dp_index])
    }
}


ifelse_check <- function(test, yes, no){if (test){return(yes)} else{ return(no) }}


find_hvg <- function(input_data = NULL, input_list = NULL, input_meta = NULL, group.by, nfeatures = 1500, ninte_features = 2000, min.cells = 5) {
    
    if (!is.null(input_list)){
        if (class(input_list[[1]]) == "Seurat"){
            seu_list <- lapply(input_list, function(seu) {
                seu <- seu  %>%
                    NormalizeData(., normalization.method = "LogNormalize", verbose = FALSE) %>%
                    FindVariableFeatures(., selection.method = "vst", nfeatures = nfeatures, verbose = FALSE) 
                seu
            })
            rm(input_list)
        } else if (grepl("matrix", class(input_list[[1]]), ignore.case = TRUE)){
            seu_list <- lapply(input_list, function(raw) {
                seu <- seu_prepare(input_raw = raw, input_norm = NULL, min.cells = min.cells, normalization.method = "LogNormalize", nfeatures = nfeatures, hvg_method = "vst", assay_name = "RNA")
                seu
                })
        }
    }


    if (class(input_data) %in% c("Seurat")){
        seu_use <- input_data; rm(input_data)
        all_samples <- levels(as.factor(seu_use@meta.data[, group.by]))
        seu_list <- lapply(all_samples, function(sample_name) {
            seu <- seu_use[, as.character(seu_use@meta.data[, group.by]) == sample_name]  %>%
                    NormalizeData(., normalization.method = "LogNormalize", verbose = FALSE) %>%
                    FindVariableFeatures(., selection.method = "vst", nfeatures = nfeatures, verbose = FALSE) 
            seu
            }) %>% setNames(., all_samples)
    } else if (grepl("matrix", class(input_data), ignore.case = TRUE)){
        if (is.null(input_meta)){
            stop("The input data is a matrix, please provide the input_meta")
        }

        all_samples <- levels(as.factor(input_meta[, group.by]))
        seu_list <- lapply(all_samples, function(sample_name){
            sample_cells <- rownames(input_meta)[as.character(input_meta[, group.by]) == sample_name]
            inseu <- seu_prepare(input_raw = input_data[, sample_cells], input_norm = NULL, min.cells = min.cells, normalization.method = "LogNormalize", nfeatures = nfeatures, hvg_method = "vst", assay_name = "RNA")
            return(inseu)
        }) %>% setNames(., all_samples)
    }


    ##Get the integration features
    feature_use <- SelectIntegrationFeatures(object.list = seu_list, nfeatures = ninte_features)
    return(feature_use)
}       


mean.of.logs <- function (x, base = 2) {
    return(log(mean((base^x) - 1) + 1, base = base))
}



library(AUCell)
GetModuleScore <- function (assay.data, features, nbin = 24, ctrl = 100, k = FALSE, seed = 42, method = c("seurat","aucell")[2], input_dir = new_inputdir, file_name, output_dir = outputdir, rethreshold_list = NULL, cellbin.size = 8000) {
     if (is.null(x = features)) {
        stop("Missing input feature list")
    }
    features <- lapply(X = features, FUN = function(x) {
        return(intersect(x = x, y = rownames(x = assay.data)))
    })
    cluster.length <- length(x = features) #number of feature list

    if (method == "seurat"){
        set.seed(seed = seed)
        pool <- rownames(x = assay.data)
        data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
        data.avg <- data.avg[order(data.avg)]
        data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
            n = nbin, labels = FALSE, right = FALSE)
        names(x = data.cut) <- names(x = data.avg)
        ctrl.use <- vector(mode = "list", length = cluster.length)
        for (i in 1:cluster.length) {
            features.use <- features[[i]]
            for (j in 1:length(x = features.use)) {
                ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
                    data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
            }
        }
        ctrl.use <- lapply(X = ctrl.use, FUN = unique)
        ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
            ncol = ncol(x = assay.data))
        for (i in 1:length(ctrl.use)) {
            features.use <- ctrl.use[[i]]
            ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, 
               ])
        }
        features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
            ncol = ncol(x = assay.data))
        for (i in 1:cluster.length) {
            features.use <- features[[i]]
            data.use <- assay.data[features.use, , drop = FALSE]
            features.scores[i, ] <- Matrix::colMeans(x = data.use)
        }
        features.scores.use <- features.scores - ctrl.scores
        features.scores.use <- as.data.frame(x = t(x = features.scores.use))
        rownames(x = features.scores.use) <- colnames(x = assay.data)
        colnames(features.scores.use) <- names(features)

        return(features.scores.use)
    } else if (method == "aucell"){
        library(AUCell)
        if (!file.exists(paste0(input_dir, file_name, "_modulescore_auc_res.rds"))){
            #Split the cells to bins[Sometimes, the matrix is too large for rankings]
            if (ncol(assay.data) < cellbin.size){
                cellbin.size <- ceiling(ncol(assay.data)/2)
            }
            bin.ind <- ceiling(c(1:ncol(assay.data))/cellbin.size)
            max.bin <- max(bin.ind)

            auc_list <- lapply(1:max.bin, function(tembin) {
                tem_matrix <- assay.data[, bin.ind == tembin]
                tem_rankings <- AUCell_buildRankings(tem_matrix, nCores=1, plotStats=FALSE) 
                tem_auc <- AUCell_calcAUC(features, tem_rankings)#, aucMaxRank = 500)
                tem_aucmatrix <- t(as.matrix(getAUC(tem_auc)))
                rm(tem_matrix, tem_rankings)
                return(tem_aucmatrix)
                })

            hauc_matrix <- do.call(rbind, auc_list)

            if (length(auc_list) == 1){ #When the input is not a named list, the colnames will be "Geneset"
                colnames(hauc_matrix) <- names(features)
            }
            saveRDS(hauc_matrix, file = paste0(input_dir, file_name, "_modulescore_auc_res.rds"))
        } else {
            hauc_matrix <- readRDS(paste0(input_dir, file_name, "_modulescore_auc_res.rds"))
        }
        
        set.seed(seed)
        pdf(paste0(output_dir, file_name, "_modulescore_auc_auto_assignment.pdf"), paper="letter")
        par(mfrow=c(2,2)); cells_assignment <- AUCell_exploreThresholds(t(hauc_matrix), plotHist=TRUE, assign=TRUE) 
        dev.off()

        #Generate assignment matrix (rownames as cells)
        default_assign <- hauc_matrix * 0 #build an empty one
        for (gset in colnames(hauc_matrix)){
            default_assign[cells_assignment[[gset]]$assignment, gset] <- 1
        }

        outdata <- list(auc = hauc_matrix, auto = default_assign, custom = default_assign)

        if (!is.null(rethreshold_list)){
            pdf(paste0(output_dir, file_name, "_modulescore_auc_custom_assignment.pdf"), paper="letter")
            for (j in names(rethreshold_list)){
                AUCell_plotHist(t(hauc_matrix)[j,,drop = FALSE], aucThr=rethreshold_list[[j]])
                abline(v=rethreshold_list[[j]])
                cells_assignment[[j]]$assignment <- rownames(hauc_matrix)[hauc_matrix[, j] >= rethreshold_list[[j]]]
            }
            dev.off()

            ##default_assign <- hauc_matrix * 0 #build an empty one
            for (gset in colnames(hauc_matrix)){
                default_assign[, gset] <- 0
                default_assign[cells_assignment[[gset]]$assignment, gset] <- 1
            }
            outdata$custom <- default_assign
        }
        return(outdata)
    }
}



auto_levels <- function(object, split.by) { 
    all_levels <- levels(as.factor(object@meta.data[, split.by]))
    level_list <- list(c("FR", "MS", "TP", "OX"), 
                        c("Frontal", "MotorSensory", "Temporal", "Occipital"), 
                        c("Human", "Chimpanzee", "Macaque", "Marmoset"),
                        c("Human", "Chimpanzee", "Rhesus", "Marmoset"),
                        c("HSB", "PTB", "RMB", "MMB"),
                        c("E37", "E42-43", "E54", "E62-E64", "E77-78"))
    ## get the idx of reference names
    nshare <- sapply(level_list, function(x) sum(all_levels %in% x))
    idx <- which(nshare == max(nshare))[1]

    if (max(nshare) == 9){
        return(all_levels)
    } else {
        new_levels <- level_list[[idx]][level_list[[idx]] %in% all_levels]
        return(new_levels)
    }
}


gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

FeatureFig <- function(object, features, cols = c("#f5f5dc", "#31a354","#253494"), split.by = NULL, plot.scale = 0.9, ncol = NULL, file_name, output_dir = outputdir, pt.size = "auto", split.order = "auto", output.ggplot = FALSE, ngradient = 10, ...) {
    
    extra_arguments <- list(...)


    if (pt.size == "auto"){
        point.size <- round(10/sqrt(ncol(object)), digits = 1) %>% MinMax(., min = 0.1, max = 1)
    } else {
        point.size <- pt.size
    }

    ##Set the levels of factors in the supervised mode
    if (!is.null(split.by)){
        ##Get the supervised order the split.by
        new_levels <- ifelse_check("auto" %in% split.order, 
                            auto_levels(object = object, split.by = split.by),
                            split.order)
        object@meta.data[, split.by] <- factor(as.character(object@meta.data[, split.by]), levels = new_levels) 
    }


    cols <- colorRampPalette(cols)(ngradient)

    ##Get the plot list
    if (length(features) == 1){
        plot_args <- c(list(object = object, cols = cols, features = features, pt.size = point.size, combine = FALSE, split.by = split.by), extra_arguments[names(extra_arguments) %in% c("order", "reduction", "min.cutoff", "max.cutoff", "shape.by", "slot", "blend","blend.threshold", "label", "label.size", "repel", "combine", "coord.fixed", "by.col", "sort.cell")])
        plist <- do.call(FeaturePlot, plot_args)
        if (!is.null(split.by)){
            names(plist) <- paste0(features, "-", rep(new_levels, length.out = length(plist)))
        } else {
            names(plist) <- features
        }
    } else {
        plot_args <- c(list(object = object, cols = cols, features = features, pt.size = point.size, combine = FALSE, split.by = split.by), extra_arguments[names(extra_arguments) %in% c("order", "reduction", "min.cutoff", "max.cutoff", "shape.by", "slot", "blend","blend.threshold", "label", "label.size", "repel", "combine", "coord.fixed", "by.col", "sort.cell")])
        plist <- do.call(FeaturePlot, plot_args)


        ##Set the names of the plot list 
        if (!is.null(split.by)){
            new_plist <- lapply(1:(length(plist)/length(new_levels)), function(x) {
                gidx <- seq(x, length(plist), length(plist)/length(new_levels))
                plist[gidx]
                }) %>% do.call(c, .)
            plist <- new_plist; rm(new_plist)

            gene_label <- sapply(plist, function(p) ggplot_build(p)$plot$plot_env$plot$layers[[1]]$mapping$colour[2] %>% as.character())
            gene_label <- paste0(gene_label, "-", rep(new_levels, length.out = length(plist)))
        } else {
            gene_label <- sapply(plist, function(p) ggplot_build(p)$plot$plot_env$plot$layers[[1]]$mapping$colour[2] %>% as.character())
        }

        names(plist) <- gene_label
    }


    ##Polish the figures by changing the theme
    plist <- lapply(names(plist), function(x) {
        p <- plist[[x]] + 
            coord_equal(ratio = 1) + 
                theme_classic() + 
                scale_color_gradientn(colors = cols, limits = c(1,ngradient)) + 
                theme(legend.position = "bottom",
                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
                    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))+ 
                labs(title =  text_wrapper(x, width=20))
        p
        })


    ##Output the plot (either a plot object or print it)
    if (output.ggplot){
        return(plist)
    } else {
        nlevel <- ifelse(is.null(split.by), 1, length(new_levels))
        ncol <- ifelse_check(is.null(ncol), nlevel, ncol)
        nrow <- ceiling(length(plist)/ncol)
        
        jpeg(paste0(output_dir, file_name, "_feature.jpeg"), width = 10 * plot.scale * ncol, height = 11 * plot.scale * nrow, units = "in", res = 300)
        plot_grid(plotlist = plist, nrow = nrow, ncol = ncol) %>% print() %>% print()
        dev.off()
    }
}



DimFig <- function(object, file_name, group.by = NULL, split.by = NULL, cols = NULL, plot.scale = 0.9, pt.size = "auto", label.size = 4, split.order = "auto", legend.position = "right", output_dir = outputdir, output.ggplot = FALSE, ...){
    p1 <- DimFig.default(object = object, file_name = file_name, group.by = group.by, split.by = split.by, cols = cols, plot.scale = plot.scale, legend.position = legend.position, label = TRUE, pt.size = pt.size, label.size = label.size, split.order = split.order, output_dir = output_dir, output.ggplot = TRUE, ...)
    p2 <- DimFig.default(object = object, file_name = file_name, group.by = group.by, split.by = split.by, cols = cols, plot.scale = plot.scale, legend.position = "none", label = TRUE, pt.size = pt.size, label.size = label.size, split.order = split.order, output_dir = output_dir, output.ggplot = TRUE, ...)
    p3 <- DimFig.default(object = object, file_name = file_name, group.by = group.by, split.by = split.by, cols = cols, plot.scale = plot.scale, legend.position = "none", label = FALSE, pt.size = pt.size, label.size = label.size, split.order = split.order, output_dir = output_dir, output.ggplot = TRUE, ...)


    plist <- lapply(1:length(group.by), function(x) list(p1[[x]], p2[[x]], p3[[x]]))

    ##Output the plot (either a plot object or print it)
    if (output.ggplot){
        return(plist %>% do.call(c, .))
    } else {

        nlevel <- ifelse_check(is.null(split.by), 1, ifelse_check("auto" %in% split.order, auto_levels(object = object, split.by = split.by),split.order) %>% length())
        nrow <- 3
        

        max_cls <- sapply(group.by, function(x) unique(object@meta.data[, x]) %>% length()) %>% max()
        pdf_width <- ifelse_check(legend.position == "none", 10 * plot.scale * nlevel, 10 * plot.scale * nlevel + ceiling(0.15 * length(max_cls)))

        for (xx in 1:length(group.by)){
            jpeg(paste0(output_dir, file_name, paste0("_", group.by[[xx]]), ifelse(is.null(split.by), "", paste0("_", split.by)), ".jpeg"), width = pdf_width, height = 10 * plot.scale * 3, units = "in", res = 300)
            plot_grid(plotlist = plist[[xx]], nrow = 3, ncol = 1) %>% print() %>% print()
            dev.off()
        }
    }
}



DimFig.default <- function(object, file_name, group.by = NULL, split.by = NULL, cols = NULL, plot.scale = 0.9, legend.position = "right", label = TRUE, pt.size = "auto", label.size = 4, split.order = "auto", output_dir = outputdir, output.ggplot = FALSE, ...){

    extra_arguments <- list(...)


    ##Set the point size
    if (pt.size == "auto"){
        point.size <- round(10/sqrt(ncol(object)), digits = 1) %>% MinMax(., min = 0.1, max = 1)
    } else {
        point.size <- pt.size
    }


    ##Set the levels of factors in the supervised mode
    if (!is.null(split.by)){
        ##Get the supervised order the split.by
        new_levels <- ifelse_check("auto" %in% split.order, 
                            auto_levels(object = object, split.by = split.by),
                            split.order)
        object@meta.data[, split.by] <- factor(as.character(object@meta.data[, split.by]), levels = new_levels) 
    }


    plot_args <- c(list(object = object, group.by = group.by, cols = cols, pt.size = point.size, combine = FALSE, split.by = split.by, label = label, label.size = label.size), extra_arguments[names(extra_arguments) %in% c("shape.by", "reduction", "cells", "repel", "cells.highlight" , "cols.highlight", "sizes.highlight", "na.value")])
    plist <- do.call(DimPlot, plot_args)


    ##Polish the figures by changing the theme
    plist <- lapply(plist, function(p) {
        p <- p + 
            #coord_equal(ratio = 1) + 
            #    theme_classic() + 
                theme(legend.position = legend.position)#,
            #        line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
            #        axis.text.x=element_blank(),axis.text.y=element_blank(), 
            #        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))+ 
            #    labs(title =  text_wrapper(x, width=20))
        p
        })


    ##Output the plot (either a plot object or print it)
    if (output.ggplot){
        return(plist)
    } else {
        nlevel <- ifelse(is.null(split.by), 1, length(new_levels))
        #ncol <- ifelse_check(is.null(ncol), nlevel, ncol)
        nrow <- length(plist)
        

        max_cls <- sapply(group.by, function(x) unique(object@meta.data[, x]) %>% length()) %>% max()
        pdf_width <- ifelse_check(legend.position == "none", 10 * plot.scale * nlevel, 10 * plot.scale * nlevel + ceiling(0.15 * length(max_cls)))

        jpeg(paste0(output_dir, file_name, ifelse(is.null(group.by), "", paste0("_", paste(group.by, collapse = "-"))), ifelse(is.null(split.by), "", paste0("_", split.by)), ifelse(label, "_labeled", ""), ".jpeg"), width = pdf_width, height = 10 * plot.scale * nrow, units = "in", res = 300)
        plot_grid(plotlist = plist, nrow = nrow, ncol = 1) %>% print() %>% print()
        dev.off()
    }
}


DotFig <- function(object, assay = "RNA", features, cols = c("white", "red"), dot.scale = 3, dot.min = 0, group.by = "hres", file_name, output_dir = outputdir) {
    p <- DotPlot(object = object, assay = assay, features = features, cols = cols, dot.scale = dot.scale,dot.min = dot.min, group.by = group.by) + 
        RotatedAxis() + 
        theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6))


    ncls <- unique(object@meta.data[, group.by]) %>% length()
    plot_heights <- ceiling(0.2 * ncls + 2) %>% MinMax(., min = 5, max = 20)
    plot_width <- ceiling(0.6 * length(features) + 2) %>% MinMax(., min = 5, max = 20)
    jpeg(paste0(output_dir, file_name, "_markers.jpeg"), width = 4.3, height = 6, units = "in", res = 300)
    print(p)
    dev.off()
}



get_ident <- function(input_ident, ident_col = NULL, file_path, sample_name, label_names = NULL){
    if (grepl("Seurat", class(input_ident), ignore.case = TRUE)){
        tem_ident <- setNames(as.character(input_ident@meta.data[, ident_col]), colnames(input_ident))
        cells <- colnames(input_ident)
    } else {
        tem_ident <- as.character(input_ident)
        cells <- names(input_ident)
    }


    tem_file <- read.table(file_path, sep="\t", stringsAsFactors=FALSE, header=TRUE)

    tem_file$cluster <- as.character(tem_file$cluster) #change the cluster numbers to characters
    sub_file <- tem_file[tem_file$sample == sample_name, ]

    if (is.null(label_names)){
        ident_names <- sub_file$label; names(ident_names) <- sub_file$cluster
        output_ident <- ident_names[tem_ident]; names(output_ident) <- cells
    } else{
        ident_names <- setNames(sub_file[, label_names], sub_file$cluster)
        output_ident <- setNames(ident_names[tem_ident], cells)
    }
    return(output_ident)
}



#This scirpt wrap the long text to several lines, designed for the title of ggplot
text_wrapper <- function(x, width=15) {
  y <- sapply(x, function(x) paste(strsplit(x, "_", fixed=TRUE)[[1]], collapse=" "))
  return(paste(strwrap(y, width=width), collapse = "\n"))
}



make_hexbin_byregion <- function(sce, nbins = 80, dimension_reduction = "UMAP", split.by = "species"){
    all_regions <- levels(as.factor(sce@meta.data[, split.by]))
    sce_list <- list()
    for (region in all_regions){
        sub_sce <- sce[, rownames(sce@meta.data)[sce@meta.data[, split.by] == region]]
        sce_list[[region]] <- make_hexbin(sub_sce, nbins = nbins, dimension_reduction = dimension_reduction)
    }
    return(sce_list)
} 


plot_hex_byregion <- function(sce, feature, action = "mean", cols = colorRampPalette(c("lightgrey", "red"))(3), split.order = c("Human", "Chimpanzee", "Rhesus", "Marmoset"), assay = "RNA"){
    if (!feature %in% rownames(sce[[1]][[assay]])) {
        stop("Specify a valid assay.")
    }

    ##First prepare the data list
    all_regions <- ifelse_check(is.null(split.order), names(sce), split.order)
    exp_list <- lapply(all_regions, function(reg) as.numeric(GetAssayData(sce[[reg]], slot = "data")[feature, ])) %>% setNames(., all_regions)
    plot.data <- lapply(all_regions, function(reg) {
        x <- exp_list[[reg]]
        out <- sce[[reg]]@misc$hexbin[[2]]
        cID <- sce[[reg]]@misc$hexbin[[1]]
        df <- data.frame(out, 
                exp = schex:::.make_hexbin_function(x, action, cID), 
                region = reg, 
                stringsAsFactors = FALSE, check.names = FALSE)
        return(df)
        }) %>% do.call(rbind, .)
    

    maxexp <- max(plot.data$exp) %>% max(., 1)

    regp_list <- list()
    for (reg in all_regions){
        sub.data <- subset(plot.data, region == reg)
        ##cur.cols <- colorRampPalette(colors = cols)(nbks)[1:max(sub.data$scale.exp)]
        regp_list[[reg]] <- ggplot(data = sub.data) +
                    geom_hex(mapping = aes_string(x = "x", y = "y", fill = "exp"), stat = "identity") + 
                    theme_classic() + 
                    labs(title = reg) + 
                    coord_fixed() + 
                    scale_fill_gradientn(colors = cols, limits = c(0, maxexp)) + 
                    theme(legend.position = "bottom", axis.line = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), legend.title = element_blank())
    }
    ##p <- plot_grid(plotlist = regp_list, nrow = 1, ncol = length(all_regions))
    
    return(regp_list)
}




plot_hex_cbn <- function(sce, feature, action = "mean", cols = viridis(3)){
    if (!feature %in% rownames(sce$RNA)) {
        stop("Specify a valid assay.")
    }

    ##First prepare the data list
    xx <- as.numeric(GetAssayData(sce, slot = "data")[feature, ])
    out <- sce@misc$hexbin[[2]]
    cID <- sce@misc$hexbin[[1]]
    plot.data <- data.frame(out, 
                exp = schex:::.make_hexbin_function(xx, action, cID), 
                stringsAsFactors = FALSE, check.names = FALSE)

    scale_exp <- function(x) {
        oo <- (x - mean(x))/sd(x)
        return(oo)
    }
    nbks <- 30
    plot.data <- plot.data %>%
                    mutate(scale.exp = scale_exp(exp)) %>%
                    mutate(scale.exp = MinMax(exp, -2.5, 2.5)) %>%
                    mutate(scale.exp = as.numeric(cut(scale.exp, nbks)))
    plot.data$scale.exp[is.na(plot.data$scale.exp)] <- 1


    p <- ggplot(data = plot.data) +
                    geom_hex(mapping = aes_string(x = "x", y = "y", fill = "scale.exp"), stat = "identity") + 
                    theme_classic() + 
                    labs(title = feature) + 
                    coord_fixed() + 
                    scale_fill_gradientn(colors = colorRampPalette(colors = cols)(nbks)) + 
                    theme(legend.position = "bottom", axis.line = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), legend.title = element_blank())
    return(p)
}



plot_hexbin_feature_byregion <- function (sce, assay, slot, features, action, nrow = 1, ncol = NULL, colors = viridis(3), file_name = NULL, output_dir, pdf_size = c(8,5), region_order = NULL, title_color = NULL, legend = "bottom", return_rawp = FALSE) {
    if (class(sce) != "list"){
        stop("The input sce should be a list output from make_hexbin_byregion")
    }
    if (!assay %in% names(sce[[1]])) {
        stop("Specify a valid assay.")
    }
    if (!slot %in% slotNames(GetAssay(sce[[1]], assay))) {
        stop("Specify a valid assay slot")
    }


    ##First prepare the data list
    all_regions <- ifelse_check(is.null(region_order), names(sce), region_order)
    features <- intersect(features, rownames(sce[[1]][[assay]]))
    if (length(features) == 0) {
            stop("No features can be matched to the dataset")
    }
    data_flist <- lapply(all_regions, function(region) GetAssayData(sce[[region]], assay = assay, slot = slot)[features, ,drop = FALSE]) %>% setNames(., all_regions)


    feature_plist <- list()
    allp_list <- list()
    for (feature in features){
        regp_list <- list()
        max_value <- 0
        exp_list <- lapply(data_flist, function(curdata) as.numeric(curdata[feature,])) %>% setNames(., names(data_flist))
        out_list <- list()
        for (region in all_regions){
            x <- exp_list[[region]]
            out <- sce[[region]]@misc$hexbin[[2]]
            cID <- sce[[region]]@misc$hexbin[[1]]
            hh <- schex:::.make_hexbin_function(x, action, cID)
            out <- as_tibble(out)
            if (grepl("^[[:digit:]]", feature)) { ##incase the feature start with a number
                feature <- paste0("F_", feature)
            }
            feature <- gsub("-", ".", feature) ##incase the "-" will be converted to "." in tibble
            col_hh <- paste0(feature, "_", action)
            eval(parse(text = paste0("out$", col_hh, " <- hh")))
            max_value <- max(max_value, max(hh))
            out_list[[region]] <- out
        }
        max_value <- ifelse(max_value == 0, 1, max_value)
        

        regp_list <- list()
        for (region in all_regions){
            regp_list[[region]] <- schex:::.plot_hexbin(out_list[[region]], colour_by = paste0(feature, "_", action), title = paste0(feature, "  ",region), xlab = xlab, ylab = ylab) + 
                scale_fill_gradientn(limits = c(0,max_value),colours = colors) + 
                theme(legend.position = legend,line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
            if (!is.null(title_color)){
                regp_list[[region]] <- regp_list[[region]] + theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold", color = title_color[region]))
            }
        }

        allp_list[[feature]] <- regp_list
        feature_plist[[feature]] <- plot_grid(plotlist = regp_list, nrow = nrow, ncol = ifelse(is.null(ncol), length(all_regions), ncol))
    }
    
    if (return_rawp){
        return(allp_list)
    }


    if (is.null(file_name)){
        return(feature_plist)
    } else {
        jpeg(paste0(output_dir, file_name, "_hexexpression.jpeg"),pdf_size[1], pdf_size[2], units = "in", res = 300)
        for (feature in names(feature_plist)){
            print(feature_plist[[feature]])
        }
        dev.off()
    }
}



###############################################################################################

## Gene expression visualization

###############################################################################################
plot_dex_heatmap <- function(data, split.by = "species", group.by = "cluster", label_genes, font_scale = 1, height_unit = 0.05, file_name, output_dir = new_outputdir, pdf_width = 12, pdf_height = 12, show_rownames = TRUE, color_limits = NULL, row_split = NULL, label_col = NULL, heat_width = 10) {
    ## Build the plot meta data
    plot_meta <- data.frame(row.names = colnames(data), species = extract_field(colnames(data), 1, "|"), 
            cluster = extract_field(colnames(data), 2, "|"), stringsAsFactors = FALSE)
    colnames(plot_meta)[colnames(plot_meta) == "species"] <- split.by
    colnames(plot_meta)[colnames(plot_meta) == "cluster"] <- group.by

    plot_meta[, split.by] <- factor(plot_meta[, split.by], levels = c("Human", "Chimpanzee", "Rhesus", "Marmoset"))
    ##plot_meta[, group.by] <- factor(plot_meta[, split.by], levels = c("Human", "Chimpanzee", "Rhesus", "Marmoset"))
    
    ## Set the colors 
    cls.anno <- readRDS(paste0(dataDir, 'cluster.color.rds'))
    sp_colors <- c("#FF420E","#4CB5F5","#89DA59","#FFBB00") %>% setNames(., c("Human", "Chimpanzee", "Rhesus", "Marmoset"))
    cls_colors <- setNames(cls.anno$color, rownames(cls.anno))
    anno_cols <- list(sp_colors[plot_meta[, split.by] %>% as.factor() %>% levels()], 
                        cls_colors[plot_meta[, group.by] %>% as.factor() %>% levels()]) %>%
                        setNames(., c(split.by, group.by))


    if (is.null(label_genes)){
        ngenes <- nrow(data)
        pdf_heights <- max(ceiling(ngenes * height_unit), 5)
        row_fontsize <- ceiling(12/ sqrt(ngenes)) * 3 * font_scale %>% MinMax(., min = 1.5, max = 12)
        pdf(paste0(output_dir, file_name, "_heatmap.pdf"), width = pdf_width, height = pdf_heights)
        pheatmap::pheatmap(data, cluster_rows =FALSE, cluster_cols = FALSE, color = colorRampPalette(viridis(3))(30), border_color = NA, annotation_col = plot_meta, annotation_colors = anno_cols, show_rownames = show_rownames, show_colnames = FALSE, fontsize_col = 8, fontsize_row = row_fontsize, gaps_col = c(nrow(plot_meta)/4, nrow(plot_meta) * 2/4, nrow(plot_meta) * 3/4))
        dev.off()
    } else {
        column_ha <- HeatmapAnnotation(df = plot_meta, col = anno_cols, annotation_height = unit(c(0.01, 0.01), "in"))
        ## Get the range of the heatmap legends
        legend_limits <- c(-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5)
        if(!is.null(color_limits)){
            color_breaks <- color_limits
        } else {
            color_breaks <- legend_limits[legend_limits >= min(data) & legend_limits <= max(data)]
        }


        ## Match gene orders
        match_order <- function(all_genes, selected_genes) {
            selected_genes[match(all_genes[all_genes %in% names(selected_genes)], names(selected_genes))]
        }

        
        font_size <- font_scale * 5
        if (is.list(label_genes)){
            ## Convert the gene labels
            ##if (is.null(names(label_genes[[1]]))){
            ##    newlabel_genes <- lapply(label_genes, function(x) setNames(x, x))
            ##} else {
            #    newlabel_genes <- label_genes
            #}
            #print(newlabel_genes)


            htlist <- rowAnnotation(link = anno_mark(at = which(rownames(data) %in% names(label_genes[[1]])), 
                        labels = match_order(rownames(data), label_genes[[1]]), 
                        side = "left",
                        labels_gp = gpar(fontsize = font_size, col = label_col[rownames(data)[rownames(data) %in% names(label_genes[[1]])]]),
                        padding = unit(1, "mm"))) + 
                    Heatmap(data, col = colorRamp2(color_breaks, viridis(length(color_breaks))), na_col = "lightgrey",
                        name = "scaled_expr", column_title = file_name, show_column_names = FALSE, width = unit(heat_width, "in"),
                        heatmap_legend_param = list(title = "Scaled expr", color_bar = "continuous", at = color_breaks), cluster_rows = FALSE, cluster_columns = FALSE, column_split = plot_meta[, split.by], row_split = row_split, top_annotation = column_ha, row_names_gp = gpar(fontsize = 5), use_raster = TRUE, raster_quality = 3) + 
                    rowAnnotation(link = anno_mark(at = which(rownames(data) %in% names(label_genes[[2]])), 
                        labels = match_order(rownames(data), label_genes[[2]]), 
                        side = "right",
                        labels_gp = gpar(fontsize = font_size, col = label_col[rownames(data)[rownames(data) %in% names(label_genes[[2]])]]),
                        padding = unit(1, "mm")))
        } else {
            ## Convert the gene labels
            htlist <- Heatmap(data, col = colorRamp2(color_breaks, viridis(length(color_breaks))), na_col = "lightgrey",
                        name = "scaled_expr", column_title = file_name, show_column_names = FALSE, width = unit(heat_width, "in"),
                        heatmap_legend_param = list(title = "Scaled expr", color_bar = "continuous", at = color_breaks), cluster_rows = FALSE, cluster_columns = FALSE, column_split = plot_meta[, split.by], row_split = row_split, top_annotation = column_ha, row_names_gp = gpar(fontsize = 5), use_raster = TRUE, raster_quality = 3) + 
                    rowAnnotation(link = anno_mark(at = which(rownames(data) %in% names(label_genes)), 
                        labels = match_order(rownames(data), label_genes), 
                        side = "right",
                        labels_gp = gpar(fontsize = font_size, col = label_col[rownames(data)[rownames(data) %in% names(label_genes)]]), padding = unit(1, "mm")))
       
        }
        pdf(paste0(output_dir, file_name, "_heatmap.pdf"), width = pdf_width, height = pdf_height)
        draw(htlist)
        dev.off() 
    }
}



##match_order(rownames(data), newlabel_genes)

## Do label transfer on the PCA dimensions
label_transfer <- function(object, reduction = "pca", dims.use = 40, transfer_cols = "label", k = 20){
    meta_use <- object@meta.data
    data.use <- object[[reduction]]@cell.embeddings[, 1:dims.use]

    
    ##A function to get the neighoring information
    get_most <- function(x){
        return(names(sort(table(x), decreasing = TRUE))[1])
    }


    ##Assume the NA values are the inquiry cell
    anno_list <- list()

    for (anno in transfer_cols){
        ref_cells <- rownames(meta_use)[!is.na(meta_use[, anno])]
        inquiry_cells <- setdiff(rownames(meta_use), ref_cells)


        ## Get the KNN cells for all the inquiry cells
        knn_cells <- FNN::get.knnx(data = data.use[ref_cells, ,drop = FALSE], query = data.use[inquiry_cells, ,drop = FALSE], k = k) %>%
                    .$nn.index


        ref_anno <- meta_use[ref_cells, anno] %>% setNames(., ref_cells)

        ## Do the annotation transfer
        if (class(ref_anno) %in% "character"){
            new_label <- apply(knn_cells, 1, function(x) get_most(ref_anno[x])) %>% 
                            setNames(., inquiry_cells) %>%
                            c(., ref_anno)
        } else if (class(ref_anno) %in% "numeric") {
            new_label <- apply(knn_cells, 1, function(x) median(ref_anno[x])) %>% 
                            setNames(., inquiry_cells) %>%
                            c(., ref_anno)
        } else {
            stop(paste0("column ", column, " has unsupported object type"))
        }

        anno_list[[anno]] <- new_label[rownames(meta_use)]
    }
    

    new_meta <- anno_list %>%  
                    as.data.frame(., stringsAsFactors = FALSE) %>% 
                    setNames(., paste0(names(anno_list), "_new"))
    object@meta.data <- cbind(object@meta.data[, setdiff(colnames(object@meta.data), colnames(new_meta))], new_meta)
    return(object)
}



