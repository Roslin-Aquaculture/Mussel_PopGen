#!/bin/python
#Name/date: TimRegan/2025.08.07

import pandas as pd

fam = pd.read_csv("SNP_v3_ELAIprep.fam", delim_whitespace=True, header=None)
fam.columns = ["FID", "IID", "PID", "MID", "Sex", "Phenotype"]

for pop in fam["FID"].unique():
    fam[fam["FID"] == pop][["FID", "IID"]].to_csv(f"{pop}.fam", sep="\t", index=False, header=False)
