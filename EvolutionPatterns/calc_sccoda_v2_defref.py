### conda activate composition

# Setup
import sys
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
import sccoda.datasets as scd



fprefix = sys.argv[1]
ctp = fprefix.split(".")[2]
ctp = ctp.split("_")[2]
ref_dic = {"Astro" : "Astro AQP4 SLC1A2", \
			"ExN" : "L3-5 RORB TNNT2 PRRX1", \
			"InN" : "InN ADARB2 CALCR", \
			"Immune" : "Micro P2RY12 APBB1IP", \
			"Oligo" : "Oligo MOG GSN", \
			"Vas" : "SMC ACTA2 CYP1B1"}
## The references were pre-determined as the cell type with least variable presence across all samples (in the function mod.CompositionalAnalysis)



cell_counts = pd.read_csv('./load_files/'+ fprefix+'.tsv', sep='\t')
data_all = dat.from_pandas(cell_counts, covariate_columns=["samplename"])
data_all.obs["species"] = data_all.obs["samplename"].str.replace(r"[0-9]", "")




model_hc = mod.CompositionalAnalysis(data_all, formula="species", reference_cell_type = ref_dic[ctp])
sim_results = model_hc.sample_hmc()
sim_results.set_fdr(est_fdr=0.2)
intcp = sim_results.intercept_df
efft = sim_results.effect_df



intcp.to_csv('./load_files/SCCODA_v2_defref/' + fprefix + '.sccoda.intercept.tsv', sep="\t")
efft.to_csv('./load_files/SCCODA_v2_defref/' + fprefix + '.sccoda.effects.tsv', sep="\t")















