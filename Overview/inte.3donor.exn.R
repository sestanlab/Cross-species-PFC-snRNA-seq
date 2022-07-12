args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")
source("./robust.fun.R")


library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 45*1000*1024^2)



pfc <- readRDS(file = paste0(dataDir, "Final_PFC_HPRC.ExN.07092021.rds"))
all_cls <- c("L2_3_IT", "L3_5_IT", "L6_IT", "L56nIT")
sp <- args[1]
ctp <- args[2]
gp <- case_when(ctp %in% c("L2_3_IT", "L3_5_IT", "L6_IT", "L56nIT") ~ "ExN",
            ctp %in% c("MGE", "CGE") ~ "InN",
            ctp %in% c("Astro", "immune", "Oligo", "Vas") ~ "NNC")


if (gp == "ExN"){
    pfc@meta.data$lres <- gsub("L3_5_IT_[0-9]", "L3_5_IT", pfc@meta.data$mres) %>%
                        gsub("L6_IT_1|L6_IT_2", "L6_IT", .) %>%
                        gsub("L5_PT|L5_6_NP|L6_CT|L6b", "L56nIT", .)
}
pfc <- subset(pfc, species == sp & lres == ctp)



## Remove one sample and check the robustness
smps <- levels(as.factor(pfc@meta.data$samplename))

for (smp in smps){
	print(paste0("Finish integration, type (removed sample: ", smp))
	ff <- paste0(inputdir, "ThreeDonorInte/ThreeSmpInte_", sp, "_", ctp, "_rm_", smp, ".rds")
	if (!file.exists(ff)){
		pfc_sub <- subset(pfc, samplename != smp)
		pfc_sub <- FastIntePCA(object = pfc_sub, split.by = "samplename", nfeatures = 1000)

		saveRDS(pfc_sub, file = ff)
	}
}







