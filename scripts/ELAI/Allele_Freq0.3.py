#!/bin/python
#Name/date: TimRegan/2025.08.05

import pandas as pd

gallo = pd.read_csv("gallo_freq.frq", delim_whitespace=True, comment="*", usecols=["SNP", "MAF"])
crom = pd.read_csv("cromarty_freq.frq", delim_whitespace=True, comment="*", usecols=["SNP", "MAF"])

df = gallo.merge(crom, on="SNP", suffixes=("_gallo", "_crom"))
df["DAF"] = abs(df["MAF_gallo"] - df["MAF_crom"])

informative = df[df["DAF"] > 0.3]["SNP"]
informative.to_csv("ELAI_informative.snplist", index=False, header=False)
