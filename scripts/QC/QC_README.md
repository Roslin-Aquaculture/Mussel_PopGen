# Preparing AxiomAnalysis Suite output for analysis
Beginning with the .ped and .map files for your dataset, follow these commands for basic QC and SNP filtration

First ensure the correct modules are loaded:
``` bash
module load igmm/apps/plink/1.90b7.2
module load roslin/admixture/1.3.0
module load anaconda
```
Activate the appropriate conda environment and =navigate to the files
```bash
conda activate orthofinder
cd SNP_070825/SNP_v1_070825
```

Update .map file to match SNP coordinates against the TOL assembly i.e. xbMytEdul2.2 or GCF_963676685.1 
**NB Ensure this is done prior to ELAI analyses, but do not run this script if planning to perform ADMIXTURE.**
Ensure that the TOL_Med_SNP Pos.txt is also present in this directory. 
Edit the python script as is necessary.
```bash
./Update_PLINK_Map.py
```

PLINK parameters v1
```bash
plink --file SNP_v1_070825 --missing --allow-extra-chr \
 --hwe 0.000001 --geno 0.05 --out SNP_v1_070825 --make-bed
```
Example output:
#V1
##60142 variants, 391 people loaded
##Total genotyping rate is 0.89512.
##37040 variants removed due to missing genotype data (--geno)
##--hwe: 3742 variants removed due to Hardy-Weinberg exact test.
##19360 variants and 391 people pass filters and QC.

```bash
plink --bfile SNP_v1_070825  --allow-extra-chr --mind 0.30  --maf 0.04 --freq \
 --set-missing-var-ids @:# --out SNP_v1_070825 --make-bed
```
Example output:
##19360 variants, 391 people loaded 
##0 people removed due to missing genotype data (--mind).
##Total genotyping rate is 0.983692.
##11750 variants removed due to minor allele threshold(s)
##7610 variants and 391 people pass filters and QC.

```bash
cut -d' ' -f1 SNP_v1_070825.fam |sort| uniq --count
```
#These are the counts for all populations listed, ensure these are correct before continuing:
     25 Aberdeen
     65 Cromarty
     36 Ireland
     18 M.edulis_Eng
     36 M.edulis_Nordic
     74 M.gallo_Atl
     28 M.gallo_Med
     12 Mtrossulus
     29 Shetland_A
     34 Shetland_B
     34 Western_Isles

