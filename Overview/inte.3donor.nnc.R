args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")
source("./robust.fun.R")


library(future)
plan("multiprocess", workers = 2)
options(future.globals.maxSize = 90*1000*1024^2)



pfc <- readRDS(file = paste0(dataDir, "Final_PFC_HPRC.All.03022022.rds"))
all_cls <- c("Astro", "Oligo", "Immune", "Vas")
sp <- args[1]
ctp <- args[2]
gp <- case_when(ctp %in% c("L2_3_IT", "L3_5_IT", "L56") ~ "ExN",
            ctp %in% c("MGE", "CGE") ~ "InN",
            ctp %in% c("Astro", "Oligo", "Immune", "Vas") ~ "NNC")


if (gp == "ExN"){
    pfc@meta.data$lres <- gsub("L3_5_IT_[0-9]", "L3_5_IT", pfc@meta.data$mres) %>%
                        gsub("L5_PT|L5_6_NP|L6_CT|L6b|L6_IT_1|L6_IT_2", "L56", .)
} else if (gp == "InN"){
	pfc@meta.data$lres <- gsub("ADARB2|LAMP5_LHX6|LAMP5_RELN|VIP", "CGE", pfc@meta.data$mres) %>%
                        gsub("PVALB|PVALB_CHC|SST|SST_NPY|TH", "MGE", .)
} else if (gp == "NNC"){
	pfc@meta.data$lres <- gsub("OPC", "Oligo", pfc@meta.data$mres) %>%
                        gsub("Micro|immune", "Immune", .) %>%
                        gsub("blood|Endo|PC|SMC|VLMC", "Vas", .)
}
pfc <- subset(pfc, species == sp & lres == ctp)
gc()


## Remove one sample and check the robustness
smps <- levels(as.factor(pfc@meta.data$samplename))

for (smp in smps){
	print(paste0("Finish integration, type (removed sample: ", smp))
	inteFile <- paste0(inputdir, "ThreeDonorInte/", "ThreeSmpInte_", sp, "_", ctp, "_rm_", smp, ".rds")
	if (!file.exists(inteFile)){
		pfc_sub <- subset(pfc, samplename != smp)
		pfc_sub <- FastIntePCA(object = pfc_sub, split.by = "samplename", nfeatures = 1000)
		saveRDS(pfc_sub, file = inteFile)
	}
}



