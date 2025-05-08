# Mussel_PopGen
Code used in analysing SNP array data representing mussel populations. 

First, ensure that all of the python scripts are in your working directory.

Ensure all conda packages are installed on whatever environment you are using (e.g. elai_env).

Change all appropriate lines in the python scripts.

The reference file "TOL_Med_SNP Pos.txt" was generated to position all SNPs along the TOL M. edulis assembly (i.e. GCF_963676685.1). Make sure it is present. 

Then run the PLINK_and_ELAI.sh bash commands.

This will output plots of Ne, pairwise Fst etc. and walk you through using ELAI.

Fimally, get summary data using the elai_visualisation.py and ELAI_EdulisLinePlot.py scripts. 
