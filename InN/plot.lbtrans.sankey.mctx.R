source("../scripts/pfc.fun.R")
library(tidyr) ## gather
library(ggsankey) 




mctx_ord <- c("Lamp5_Lhx6","Lamp5_Plch2_Dock5","Lamp5_Lsp1","Vip_Ptprt_Pkp2","Vip_Chat_Htr1f","Vip_Col15a1_Pde1a","Vip_Lmo1_Fam159b","Lamp5_Fam19a1_Tmem182","Vip_Crispld2_Htr2c","Vip_Crispld2_Kcne4","Vip_Arhgap36_Hmcn1","Vip_Gpc3_Slc18a3","Vip_Lect1_Oxtr","Vip_Lmo1_Myl1","Vip_Pygm_C1ql1","Vip_Rspo1_Itga4","Vip_Rspo4_Rxfp1_Chat","Vip_Igfbp4_Mab21l1","Vip_Igfbp6_Car10","Vip_Igfbp6_Pltp","Lamp5_Fam19a1_Pax6","Serpinf1_Aqp5_Vip","Serpinf1_Clrn1","Sncg_Gpr50","Sncg_Slc17a8","Sncg_Vip_Itih5","Sncg_Vip_Nptx2","Lamp5_Krt73","Lamp5_Ntn1_Npy2r","Sst_Chodl","Sst_Crh_4930553C11Rik_","Sst_Crhr2_Efemp1","Sst_Calb2_Necab1","Sst_Calb2_Pdlim5","Sst_Tac1_Htr1d","Sst_Hpse_Cbln4","Sst_Mme_Fam114a1","Sst_Nr2f2_Necab1","Sst_Hpse_Sema3c","Sst_Esm1","Sst_Tac2_Myh4","Sst_Rxfp1_Eya1","Sst_Rxfp1_Prdm8","Sst_Tac2_Tacstd2","Sst_Chrna2_Ptgdr","Sst_Chrna2_Glra3","Sst_Myh8_Etv1_","Sst_Nts","Pvalb_Th_Sst","Pvalb_Gabrg1","Sst_Myh8_Fibin","Sst_Tac1_Tacr3","Pvalb_Reln_Itm2a","Pvalb_Gpr149_Islr","Pvalb_Sema3e_Kank4","Pvalb_Reln_Tac1","Pvalb_Tpbg","Pvalb_Calb1_Sst","Pvalb_Akr1c18_Ntf3", "Pvalb_Vipr2") %>%
  gsub("_", " ", .)
pfc_ord <- readRDS(paste0(dataDir, "PFC.cluster.order.rds"))[c(41:42, 44:86)]

mctx_mres <- c(strsplit(mctx_ord, " ", fixed = TRUE) %>% sapply(., function(x) x[1])) %>%
  toupper()
mctx_col <- sapply(mctx_mres, function(x) {
  hue <- switch(x, LAMP5 = "purple", VIP = "blue", SERPINF1 = "green", SNCG = "green", ADARB2 = "green", SST = "orange", PVALB = "red", LHX6 = "red")
  color = randomcoloR::randomColor(count = 1, hue = hue, luminosity = "bright")
  color
}) %>%
	setNames(., mctx_ord)


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

 

meta <- readRDS(file = paste0(inputdir, "HDFC_transed_by_MCTX.meta.rds"))
meta <- meta[, c("fig1cluster", "predicted.id")]
## Filter some noise
meta <- FilterTransfer(meta = meta, ref_col = "fig1cluster", query_col = "predicted.id", thre = 0.01)
mctx_use <- intersect(mctx_ord, meta$predicted.id)



source("./inn.fun.R")
p <- plot_ggsankey(meta = meta, x_col = "fig1cluster", next_col = "predicted.id", x_ord = pfc_ord, next_ord = mctx_use, type = "sankey", space= 0.01, width= 0.1, x_cols = pfc_col, next_cols = mctx_col[mctx_use])

pdf(paste0(outputdir, "LabelTransfer_sankey_MCTX_DFC_v2.pdf"), width = 10, height = 10)
print(p)
dev.off()















