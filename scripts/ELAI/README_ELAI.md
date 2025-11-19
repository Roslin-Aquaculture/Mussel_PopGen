# ELAI Local Ancestry Inference

This module runs ELAI on the **v3 dataset** using reference panels defined in metadata.

---
Load modules and activate conda env
```bash
module load igmm/apps/plink/1.90b7.2
module load roslin/admixture/1.3.0
module load anaconda
conda activate elai_env
```
# Ensure we have the files we need
```bash
cd SNP_070825/SNP_v3_070825_ELAI/
cp ../SNP_v3_070825/*.ped .
cp ../SNP_v3_070825/*.map .
cp ../SNP_v3_070825/TOL_ .
cp ../SNP_v3_070825/TOL_Med_SNP Pos.txt .
```
PLINK filter:
```bash
plink --file SNP_v3_070825 --geno 0.05 --mind 0.30 --maf 0.10 \
--chr 1-14 --make-bed --out SNP_v3_ELAIprep
```
Example output:

#49494 variants 351 people loaded from .bim file.

#Total genotyping rate is 0.897857.

#29744 variants removed due to missing genotype data (--geno).

#14025 variants removed due to minor allele threshold(s)

#5725 variants and 351 people pass filters and QC.

# CUstom Python script to split according to population 
(ensure to edit this file as is required)
```bash
./split_fam_by_pop.py
```

# Get allele frequencies for the two reference groups
```bash
plink --bfile SNP_v3_ELAIprep --keep Cromarty.fam     --freq --out cromarty_freq
plink --bfile SNP_v3_ELAIprep --keep M.gallo_Atl.fam --freq --out gallo_freq
```

# Calculate Difference in allele freq (ΔAF) and extract informative SNPs
filtering for ancestry-informative markers — specifically SNPs with:

∣MAFref1 − MAFref2∣ > 0.3 
```bash
./Allele_Freq0.3.py
```
or with ∣MAFref1 − MAFref2∣ > 0.2 :
```bash
./Allele_Freq0.2.py
```

# Final ELAI dataset
```bash
plink --bfile SNP_v3_ELAIprep --extract ELAI_informative.snplist \
 --make-bed --out SNP_ELAI_final
```

# Create subdatasets for each population
```bash
for pop in $(awk '{print $1}' SNP_ELAI_final.fam | sort | uniq); do
    plink --bfile SNP_ELAI_final --keep <(awk -v p=$pop '$1 == p {print $1, $2}' SNP_ELAI_final.fam) \
          --make-bed --out ${pop}_subset 
done
```

# Convert plink files to ELAI format:
Need to change all alllele to single character first for bimbam
```bash
plink --bfile Aberdeen_subset --snps-only just-acgt --biallelic-only strict --make-bed  --out Aberdeen_clean
plink --bfile Aberdeen_clean --recode bimbam --out Aberdeen_elai
plink --bfile M.gallo_Atl_subset --snps-only just-acgt --biallelic-only strict --make-bed  --out M.gallo_Atl_clean
plink --bfile M.gallo_Atl_clean --recode bimbam --out M.gallo_Atl_elai
plink --bfile Cromarty_subset --snps-only just-acgt --biallelic-only strict --make-bed  --out Cromarty_clean
plink --bfile Cromarty_clean --recode bimbam --out Cromarty_elai
plink --bfile M.edulis_Eng_subset --snps-only just-acgt --biallelic-only strict --make-bed  --out M.edulis_Eng_clean
plink --bfile M.edulis_Eng_clean --recode bimbam --out M.edulis_Eng_elai
plink --bfile M.edulis_Eng_subset --snps-only just-acgt --biallelic-only strict --make-bed  --out M.edulis_Eng_clean
plink --bfile M.edulis_Eng_clean --recode bimbam --out M.edulis_Eng_elai
plink --bfile Ireland_subset --snps-only just-acgt --biallelic-only strict --make-bed  --out Ireland_clean
plink --bfile Ireland_clean --recode bimbam --out Ireland_elai
plink --bfile Shetland_A_subset --snps-only just-acgt --biallelic-only strict --make-bed  --out Shetland_A_clean
plink --bfile Shetland_A_clean --recode bimbam --out Shetland_A_elai
plink --bfile Shetland_B_subset --snps-only just-acgt --biallelic-only strict --make-bed  --out Shetland_B_clean
plink --bfile Shetland_B_clean --recode bimbam --out Shetland_B_elai
plink --bfile Western_Isles_subset --snps-only just-acgt --biallelic-only strict --make-bed  --out Western_Isles_clean
plink --bfile Western_Isles_clean --recode bimbam --out Western_Isles_elai
```

Run ELAI for these poplationsusing M. gallo Atlantic as source population 1 and Cromarty (M. edulis) as source population 2
 ```bash
../../ELAI/elai-lin -g M.gallo_Atl_elai.recode.geno.txt -p 10 -g Cromarty_elai.recode.geno.txt -p 11 \
-g Aberdeen_elai.recode.geno.txt -p 1 \
-pos Aberdeen_elai.recode.pos.txt \
-C 2 -c 10 -s 30 -mg 50 -o ELAI_Results_Aberdeen

../../ELAI/elai-lin -g M.gallo_Atl_elai.recode.geno.txt -p 10 -g Cromarty_elai.recode.geno.txt -p 11 \
-g M.edulis_Eng_elai.recode.geno.txt -p 1 \
-pos M.edulis_Eng_elai.recode.pos.txt \
-C 2 -c 10 -s 30 -mg 50 -o ELAI_Results_M.edulis_Eng
 
../../ELAI/elai-lin -g M.gallo_Atl_elai.recode.geno.txt -p 10 -g Cromarty_elai.recode.geno.txt -p 11 \
-g Shetland_A_elai.recode.geno.txt -p 1 \
-pos Shetland_A_elai.recode.pos.txt \
-C 2 -c 10 -s 30 -mg 50 -o ELAI_Results_Shetland_A
 
../../ELAI/elai-lin -g M.gallo_Atl_elai.recode.geno.txt -p 10 -g Cromarty_elai.recode.geno.txt -p 11 \
-g Shetland_B_elai.recode.geno.txt -p 1 \
-pos Shetland_B_elai.recode.pos.txt \
-C 2 -c 10 -s 30 -mg 50 -o ELAI_Results_Shetland_B

../../ELAI/elai-lin -g M.gallo_Atl_elai.recode.geno.txt -p 10 -g Cromarty_elai.recode.geno.txt -p 11 \
-g Ireland_elai.recode.geno.txt -p 1 \
-pos Ireland_elai.recode.pos.txt \
-C 2 -c 10 -s 30 -mg 50 -o ELAI_Results_Ireland
 
../../ELAI/elai-lin -g M.gallo_Atl_elai.recode.geno.txt -p 10 -g Cromarty_elai.recode.geno.txt -p 11 \
-g Western_Isles_elai.recode.geno.txt -p 1 \
-pos Western_Isles_elai.recode.pos.txt \
-C 2 -c 10 -s 30 -mg 50 -o ELAI_Results_Western_Isles
 ```
 
This produced ELAI output files in an output folder with extensions:

 .snpinfo.txt, .admix.txt, .ps21.txt
 
 and also .em.txt and .log.txt which I did not use.

# FInally, visualise output:
Custom scripts used:
```bash
./SNP_DensityGenomePlot.py
cd output/
../elai_admix_viz.py
../ELAI_EdulisLinePlot.py
```
