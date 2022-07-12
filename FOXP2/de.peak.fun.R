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


overlap_foxp2_peaks <- function(peaks){    
    # Process peak names chrN-ddd-ddd to dataframe to operate on
    pinfo <- data.frame(do.call('rbind', lapply(peaks, function(x) {
        setNames(strsplit(x, '-', fixed=T)[[1]], c('chr', 'start', 'end'))})),
                        row.names=peaks)

    # The fragments might be repeated and contain a '.1' suffix in the end position
    pinfo$end <- gsub('.1', '', pinfo$end, fixed=T)
    pinfo[,c('start', 'end')] <- apply(pinfo[,c('start', 'end')], 2, as.numeric)

    # Compute masks: we need matching chromosome and at least start or end inside the FOXP2 region (zoomed out x1.5)
    pinfo$chr.mask <- pinfo$chr==foxp2_chr
    pinfo$start.mask <- pinfo$start > foxp2_start & pinfo$start < foxp2_end
    pinfo$end.mask <- pinfo$end > foxp2_start & pinfo$end < foxp2_end
    pinfo$pos.mask <- apply(pinfo[c('start.mask', 'end.mask')], 1, any)
    pinfo$peak.mask <- apply(pinfo[c('chr.mask', 'pos.mask')], 1, all)

    return(peaks[pinfo$peak.mask])
}


