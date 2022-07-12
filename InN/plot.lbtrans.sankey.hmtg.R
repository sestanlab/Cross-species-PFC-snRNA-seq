source("../scripts/pfc.fun.R")
library(tidyr) ## gather
library(ggsankey) 




nameupdate <- function(x){
  x <- x %>% gsub("Inh L1 SST NMBR","Inh L1 LAMP5 NMBR",.) %>%
    gsub("Inh L5-6 GAD1 GLP1R","Inh L5-6 LHX6 GLP1R",.) %>%
    gsub("Inh L1-3 PAX6 SYT6", "Inh L1-3 VIP SYT6", .) %>%
    gsub("Inh L1-2 GAD1 MC4R", "Inh L1-2 ADARB2 MC4R",.)
  x  
}


mtg_ord <- c("Inh_L2-6_LAMP5_CA1","Inh_L1_SST_NMBR","Inh_L1-4_LAMP5_LCP2","Inh_L1-2_LAMP5_DBP","Inh_L2-6_VIP_QPCT","Inh_L3-6_VIP_HS3ST3A1","Inh_L1-3_VIP_GGH","Inh_L1-4_VIP_OPRM1","Inh_L2-4_VIP_SPAG17","Inh_L1-2_VIP_PCDH20","Inh_L2-5_VIP_TYR","Inh_L1-4_VIP_CHRNA6","Inh_L2-3_VIP_CASC6","Inh_L1-3_VIP_CCDC184","Inh_L2-4_VIP_CBLN1","Inh_L2-5_VIP_SERPINF1","Inh_L1-3_VIP_CHRM2","Inh_L1-3_VIP_ADAMTSL1","Inh_L1-2_PAX6_CDH12","Inh_L1-2_PAX6_TNFAIP8L3","Inh_L1-2_VIP_TSPAN12","Inh_L1-3_PAX6_SYT6","Inh_L1-4_VIP_PENK", "Inh_L1-2_VIP_LBH","Inh_L1_SST_CHRNA4","Inh_L1-2_GAD1_MC4R", "Inh_L3-6_SST_NPY","Inh_L1-2_SST_BAGE2","Inh_L3-5_SST_ADGRG6","Inh_L2-4_SST_FRZB","Inh_L1-3_SST_CALB1","Inh_L4-6_SST_GXYLT2","Inh_L5-6_SST_NPM1P10","Inh_L4-5_SST_STK32A","Inh_L3-6_SST_HPGD","Inh_L4-6_SST_B3GAT2","Inh_L5-6_SST_KLHDC8A","Inh_L5-6_SST_TH","Inh_L5-6_GAD1_GLP1R","Inh_L5-6_PVALB_LGR5","Inh_L2-4_PVALB_WFDC2","Inh_L5-6_SST_MIR548F2","Inh_L4-5_PVALB_MEPE","Inh_L4-6_PVALB_SULF1","Inh_L2-5_PVALB_SCUBE3") %>%
		gsub("_", " ", .) %>%
		nameupdate()
pfc_ord <- readRDS(paste0(dataDir, "PFC.cluster.order.rds"))[c(41:42, 44:86)]

mtg_mres <- c(strsplit(mtg_ord, " ", fixed = TRUE) %>% sapply(., function(x) x[3])) %>%
  toupper()
mtg_col <- sapply(mtg_mres, function(x) {
  hue <- switch(x, LAMP5 = "purple", VIP = "blue", PAX6 = "green", SNCG = "green", ADARB2 = "green", SST = "orange", PVALB = "red", LHX6 = "red")
  color = randomcoloR::randomColor(count = 1, hue = hue, luminosity = "bright")
  color
}) %>%
    setNames(., names(mtg_mres))


cls_col <- readRDS(paste0(dataDir, "cluster.color.rds"))
pfc_col <- setNames(cls_col$color, rownames(cls_col)) %>% .[pfc_ord]



FilterTransfer <- function(meta, ref_col, query_col, thre = 0.02) {
  meta$pair <- paste0(meta[, ref_col], "|", meta[, query_col])
  kp_pairs <- meta %>% 
      group_by(!!sym(ref_col), !!sym(query_col)) %>%
      summarize(nhits = n()) %>%
      ungroup() %>%
      group_by(!!sym(ref_col)) %>%
      mutate(ratio = nhits/sum(nhits)) %>%
      filter(ratio >= thre) %>%
      mutate(pair = paste0(!!sym(ref_col), "|", !!sym(query_col))) %>%
      .$pair

  meta <- meta %>%
      filter(pair %in% kp_pairs) %>%
      select(-pair)
  return(meta)
}




##seu <- readRDS(file = paste0(inputdir, "HMTG_transed_by_HDFC.transfer.rds")) 
meta <- readRDS(file = paste0(inputdir, "HDFC_transed_by_HMTG.meta.rds"))
meta <- meta[, c("fig1cluster", "predicted.id")] %>%
          mutate(predicted.id = nameupdate(predicted.id))
meta <- FilterTransfer(meta = meta, ref_col = "fig1cluster", query_col = "predicted.id", thre = 0.05)
mtg_use <- intersect(mtg_ord, meta$predicted.id)
print(length(mtg_use))



source("./inn.fun.R")
p <- plot_ggsankey(meta = meta, x_col = "fig1cluster", next_col = "predicted.id", x_ord = pfc_ord, next_ord = mtg_use, type = "sankey", space= 0.01, width= 0.1, x_cols = pfc_col, next_cols = mtg_col[mtg_use])

pdf(paste0(outputdir, "LabelTransfer_sankey_HMTG_DFC.pdf"), width = 10, height = 10)
print(p)
dev.off()



