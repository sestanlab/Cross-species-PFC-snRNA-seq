library(foreach)
###############################################################################################

## Wilcox Test Part [test on all genes]

###############################################################################################
##------------------------------------------------------------------------------------------
## Do the marker analysis use all the expressed genes.
allWilcox.default <- function(object, group.by = "hres", ident_1, ident_2 = NULL, pseudocount.use = 0, adjust_method="bonferroni", max.cells.per.ident = 5000) {
    input_identity <- as.character(object@meta.data[, group.by])
  
    #Extract the data for the DEX analysis
    data_use <- object$RNA@data


    #get the cell names for two comparing groups. 
    cells_1 <- colnames(object)[object@meta.data[, group.by] == ident_1]
    if (is.null(x = ident_2)) {
        cells_2 <- setdiff(colnames(data_use), cells_1)
    } else {
        cells_2 <- colnames(object)[object@meta.data[, group.by] %in% ident_2] 
    }
  

    ## In case cells_2 is two large
    if (max.cells.per.ident < Inf) {
        set.seed(seed = 0)
        if (length(cells_1) > max.cells.per.ident)
            cells_1 = sample(x = cells_1, size = max.cells.per.ident)
        if (length(cells_2) > max.cells.per.ident)
            cells_2 = sample(x = cells_2, size = max.cells.per.ident)
    }


    ## Get the ratio
    pct.1 <- Matrix::rowMeans(data_use[, cells_1,drop = FALSE] != 0) %>% round(., digits = 3)
    pct.2 <- Matrix::rowMeans(data_use[, cells_2,drop = FALSE] != 0) %>% round(., digits = 3)
    data.alpha <- cbind(pct.1, pct.2)


    ## Get the average expression 
    data.1 <- Matrix::rowMeans(exp(data_use[, cells_1, drop = FALSE]) - 1)
    data.2 <- Matrix::rowMeans(exp(data_use[, cells_2, drop = FALSE]) - 1)
    total.diff <- log((data.1 + pseudocount.use)/(data.2 + pseudocount.use))
    if (TRUE) {
        coldata <- data.frame(t(as.matrix(data_use[1:2,c(cells_1,cells_2)])))
        colnames(coldata) <- c("gene1","gene2")
        group <- paste(sample(x = letters, size = 6), collapse = "")
        coldata[cells_1, group] <- "Group1"
        coldata[cells_2, group] <- "Group2"
        coldata[, group] <- factor(x = coldata[, group])
        coldata$wellKey <- rownames(x = coldata)
        countdata_test <- data_use[, rownames(x = coldata), drop = FALSE]
        print(dim(countdata_test))
        
        p_val <- sapply(X = 1:nrow(countdata_test), FUN = function(x) wilcox.test(countdata_test[x, ] ~ coldata[, group])$p.value)
        names(p_val) <- rownames(countdata_test)
        to_return <- data.frame(p_val = p_val, row.names = rownames(x = countdata_test))
    }
    to_return[, "avg_logFC"] <- total.diff[rownames(x = to_return)]
    to_return <- cbind(to_return, data.alpha[rownames(x = to_return), ,drop = FALSE])
    to_return$p_val[is.na(to_return$p_val)] <- 1
    to_return$avg_logFC[is.na(to_return$avg_logFC)] <- 0
    to_return$p_val_adj = p.adjust(p = to_return$p_val, method = adjust_method, n = nrow(data_use))

    to_return <- to_return[order(to_return$p_val, -to_return$avg_logFC),]
    to_return$gene <- rownames(to_return)
    to_return$cluster <- ident_1
    return(to_return)
}



allWilcox <- function(object, group.by = "hres", ident_1, ident_2 = NULL, pseudocount.use = 0, adjust_method="bonferroni", nCores = 8, max.cells.per.ident = 5000){
    ## Make sure the ident_1 & ident_2 are list
    if (class(ident_1) != "list" || class(ident_2) != "list" ){
        stop("Please make the input ident_1 & ident_2 are list objects")
    }
    if (length(ident_1) != length(ident_2)){
        stop("Please make sure the length of ident_1 & ident_2 are the same")
    }


    ## For each element in the list, do the Wilcox test
    all_cls <- names(ident_1)


    print(Sys.time())
    cl = makeCluster(nCores, outfile="")
    doParallel::registerDoParallel(cl);
    tem_idx <- NULL

    toreturn_res <- foreach(tem_idx = 1:length(all_cls), .combine = rbind,.packages = c("stats", "Matrix", "Seurat", "dplyr"), .export = c("allWilcox.default")) %dopar% {
            cls_use <- unique(c(ident_1[[tem_idx]], ident_2[[tem_idx]]));
            seu <- object[, object@meta.data[, group.by] %in% cls_use];
            res <- allWilcox.default(object = seu, group.by = group.by, ident_1 = ident_1[[tem_idx]], ident_2 = ident_2[[tem_idx]], pseudocount.use = pseudocount.use, adjust_method="bonferroni", max.cells.per.ident = max.cells.per.ident)
            res$bg <- paste(ident_2[[tem_idx]], collapse = "|")

            print(paste0("Finish DEX for cluster: ", all_cls[tem_idx]))
            return(res)
        }
    stopCluster(cl)
    print(paste0("Finish DEX for all clusters"))
    print(Sys.time())
    return(toreturn_res)
}


