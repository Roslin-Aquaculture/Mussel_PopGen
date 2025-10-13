#Add modules/programmes required
module load igmm/apps/plink/1.90b7.2
module load roslin/admixture/1.3.0
module load anaconda
conda activate orthofinder

cd /exports/cmvm/eddie/eb/groups/Regan_grp/SNP_070825/SNP_v1_070825

#Ensure that TOL_Med_SNP Pos.txt is present to map SNPs to TOL M. edulis assembly
./Update_PLINK_Map.py

#Stringent parameters v1
plink --file SNP_v1_070825 --missing --allow-extra-chr \
 --hwe 0.000001 --geno 0.05 --out SNP_v1_070825 --make-bed
#V1
##60142 variants, 391 people loaded
##Total genotyping rate is 0.89512.
##37040 variants removed due to missing genotype data (--geno)
##--hwe: 3742 variants removed due to Hardy-Weinberg exact test.
##19360 variants and 391 people pass filters and QC.

#Stringent parameters v1
plink --bfile SNP_v1_070825  --allow-extra-chr --mind 0.30  --maf 0.04 --freq \
 --set-missing-var-ids @:# --out SNP_v1_070825 --make-bed
##19360 variants, 391 people loaded 
##0 people removed due to missing genotype data (--mind).
##Total genotyping rate is 0.983692.
##11750 variants removed due to minor allele threshold(s)
##7610 variants and 391 people pass filters and QC.

cut -d' ' -f1 SNP_v1_070825.fam |sort| uniq --count
#These are the counts for all populations listed
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


plink --bfile SNP_v1_070825 --double-id --allow-extra-chr --set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out SNP_v1_070825 --make-bed
#Total genotyping rate is 0.975963.
#7610 variants and 391 people pass filters and QC.
#--indep-pairwise: Ignoring 1435 chromosome 0 variants.
#Pruning complete.  1891 of 6175 variants removed.

plink --bfile SNP_v1_070825 --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract SNP_v1_070825.prune.in --make-bed --pca --out SNP_v1_070825
#--extract: 4284 variants remaining.
#Total genotyping rate is 0.0.974869.
#4284 variants and 391 people pass filters and QC.

admixture --cv SNP_v1_070825.bed 2 > log2.out

for i in {2..5}
do
admixture --cv SNP_v1_070825.bed $i > log${i}.out
done

awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > SNP_v1_070825.cv.error
#2 0.36935
#3 0.36541
#4 0.36533
#5 0.36698
#k=4

cut -d' ' -f2 SNP_v1_070825.fam > SNP_v1_070825.txt
########################################################################################

#Stringent parameters v2
plink --file SNP_v2_070825 --missing --allow-extra-chr \
 --hwe 0.000001 --geno 0.05 --out SNP_v2_070825 --make-bed
##60142 variants, 379 people loaded
##Total genotyping rate is 0.894266.
##37477 variants removed due to missing genotype data (--geno)
##--hwe: 3628 variants removed due to Hardy-Weinberg exact test.
##19037 variants and 379 people pass filters and QC.

#Stringent parameters v2
plink --bfile SNP_v2_070825  --allow-extra-chr --mind 0.30  --maf 0.04 --freq \
 --set-missing-var-ids @:# --out SNP_v2_070825 --make-bed
##19360 variants, 379 people loaded 
##0 people removed due to missing genotype data (--mind).
##Total genotyping rate is 0.984856..
##11567 variants removed due to minor allele threshold(s)
##7470 variants and 379 people pass filters and QC.

cut -d' ' -f1 SNP_v2_070825.fam |sort| uniq --count
#These are the counts for all populations listed
     25 Aberdeen
     65 Cromarty
     36 Ireland
     18 M.edulis_Eng
     36 M.edulis_Nordic
     74 M.gallo_Atl
     28 M.gallo_Med
     29 Shetland_A
     34 Shetland_B
     34 Western_Isles

plink --bfile SNP_v2_070825 --double-id --allow-extra-chr --set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out SNP_v2_070825 --make-bed
#Total genotyping rate is 0.977312.
#7470 variants and 379 people pass filters and QC.
#--indep-pairwise: Ignoring 1405 chromosome 0 variants.
#Pruning complete.  1798 of 6065 variants removed.

plink --bfile SNP_v2_070825 --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract SNP_v2_070825.prune.in --make-bed --pca --out SNP_v2_070825
#--extract: 4267 variants remaining.
#Total genotyping rate is 0.9762.
#4267 variants and 379 people pass filters and QC.

admixture --cv SNP_v2_070825.bed 2 > log2.out

for i in {2..5}
do
admixture --cv SNP_v2_070825.bed $i > log${i}.out
done

awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > SNP_v2_070825.cv.error
#2 0.36993
#3 0.36860
#4 0.37000
#5 0.37130
#k=3 

cut -d' ' -f2 SNP_v2_070825.fam > SNP_v2_070825.txt
########################################################################################


#Stringent parameters
plink --file SNP_v3_070825 --missing --allow-extra-chr \
 --hwe 0.000001 --geno 0.05 --out SNP_v3_070825 --make-bed
##60142 variants, 351 people loaded
##Total genotyping rate is 0.89816.
##35929 variants removed due to missing genotype data (--geno)
##--hwe: 3771 variants removed due to Hardy-Weinberg exact test.
##20442 variants and 351 people pass filters and QC.

#Stringent parameters
plink --bfile SNP_v3_070825  --allow-extra-chr --mind 0.30  --maf 0.04 --freq \
 --set-missing-var-ids @:# --out SNP_v3_070825 --make-bed
##20442 variants, 351 people loaded 
##0 people removed due to missing genotype data (--mind).
##Total genotyping rate is 0.984598.
##12382 variants removed due to minor allele threshold(s)
##8060 variants and 351 people pass filters and QC.

cut -d' ' -f1 SNP_v3_070825.fam |sort| uniq --count
#These are the counts for all populations listed
     25 Aberdeen
     65 Cromarty
     36 Ireland
     18 M.edulis_Eng
     36 M.edulis_Nordic
     74 M.gallo_Atl
     29 Shetland_A
     34 Shetland_B
     34 Western_Isles

plink --bfile SNP_v3_070825 --double-id --allow-extra-chr --set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out SNP_v3_070825 --make-bed
#Total genotyping rate is 0.976875.
#7470 variants and 351 people pass filters and QC.
#--indep-pairwise: Ignoring 1498 chromosome 0 variants.
#Pruning complete.  2049 of 6562 variants removed.

plink --bfile SNP_v3_070825 --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract SNP_v3_070825.prune.in --make-bed --pca --out SNP_v3_070825
#--extract: 4513 variants remaining.
#Total genotyping rate is 0.97614.
#4513 variants and 351 people pass filters and QC.

admixture --cv SNP_v3_070825.bed 2 > log2.out

for i in {2..5}
do
admixture --cv SNP_v3_070825.bed $i > log${i}.out
done

awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > SNP_v3_070825.cv.error
#2 0.37457
#3 0.37462
#4 0.37590
#5 0.38107
#k=2

cut -d' ' -f2 SNP_v3_070825.fam > SNP_v3_070825.txt



##########################################

#Calculate pairwise Fst using PLINK
for pop1 in $(awk '{print $1}' SNP_v3_070825.fam | sort -u); do
    for pop2 in $(awk '{print $1}' SNP_v3_070825.fam | sort -u); do
        if [ "$pop1" != "$pop2" ]; then
            echo "Calculating Fst between $pop1 and $pop2"
            # Generate a properly formatted within file (FID IID GroupID)
            awk -v p1="$pop1" -v p2="$pop2" '
                $1 == p1 {print $1, $2, 1}
                $1 == p2 {print $1, $2, 2}
            ' SNP_v3_070825.fam > temp_within.txt
            # Debug: Check the first few lines
            head temp_within.txt
            # Run PLINK Fst
            plink --bfile SNP_v3_070825 --fst --within temp_within.txt --allow-extra-chr --out fst_${pop1}_${pop2}
            # Clean up
            rm temp_within.txt
        fi
    done
done
 
# Create subdatasets for each population
for pop in $(awk '{print $1}' SNP_v3_070825.fam | sort | uniq); do
    plink --bfile SNP_v3_070825 --keep <(awk -v p=$pop '$1 == p {print $1, $2}' SNP_v3_070825.fam) \
          --make-bed --out ${pop}_subset --allow-extra-chr
done
 
#Calculate O(HOM) (Observed homozygous genotypes), E(HOM) (Expected homozygous genotypes) and F (Inbreeding coefficient)
for pop in $(awk '{print $1}' SNP_v3_070825.fam | sort | uniq); do
  echo "Processing population: $pop"
  plink --bfile "$pop"_subset --allow-extra-chr --het --out $pop
done
  
#Calculate Effective Population Size (Ne)
for pop in $(awk '{print $1}' SNP_v3_070825.fam | sort | uniq); do
  echo "Processing population: $pop"
  plink --bfile "$pop"_subset --allow-extra-chr --r2 gz --ld-window-kb 1000 --out $pop
done
  
#Calculate Linkage Disequilibrium (LD) Decay
for pop in $(awk '{print $1}' SNP_v3_070825.fam | sort | uniq); do
  echo "Processing population: $pop"
  plink --bfile "$pop"_subset --allow-extra-chr --r2 --ld-window-r2 0 --ld-window 99999 --ld-window-kb 1000 --out $pop
  awk '{print $2, $5, $7}' $pop.ld > $pop.ld_summary
done
 
#Calculate Fst: Creates clusters according to Family ID:
plink --bfile SNP_v3_070825 --fst --family --allow-extra-chr --out fst_results
 
#Updated this script, see below.
./Summary_Plots.py
#######################################################
##########################################################
## Moving onto ELAI

conda deactivate
conda activate elai_env

#Repeated with new parameters for ELAI

module load igmm/apps/plink/1.90b7.2
module load roslin/admixture/1.3.0
module load anaconda
conda activate elai_env

cd /exports/cmvm/eddie/eb/groups/Regan_grp/SNP_070825/SNP_v3_070825_ELAI2/
cp ../SNP_v3_070825/*.ped .
cp ../SNP_v3_070825/*.map .
cp ../SNP_v3_070825/TOL_ .
cp ../SNP_v3_070825/TOL_Med_SNP\ Pos.txt .

plink --file SNP_v3_070825 --geno 0.05 --mind 0.30 --maf 0.10 \
--chr 1-14 --make-bed --out SNP_v3_ELAIprep
#49494 variants 351 people loaded from .bim file.
#Total genotyping rate is 0.897857.
#29744 variants removed due to missing genotype data (--geno).
#14025 variants removed due to minor allele threshold(s)
#5725 variants and 351 people pass filters and QC.

#Relaxed
plink --file SNP_v3_070825 --geno 0.1 --mind 0.30 --maf 0.04 \
--chr 1-14 --make-bed --out SNP_v3_ELAIprep
#49494 variants 351 people loaded from .bim file.
#Total genotyping rate is 0.897857.
#20355 variants removed due to missing genotype data (--geno).
#11432 variants removed due to minor allele threshold(s)
#17707 variants and 351 people pass filters and QC.


./split_fam_by_pop.py

# Get allele frequencies for the two reference groups
plink --bfile SNP_v3_ELAIprep --keep Cromarty.fam     --freq --out cromarty_freq
plink --bfile SNP_v3_ELAIprep --keep M.gallo_Atl.fam --freq --out gallo_freq

#Calculate Difference in allele freq (ΔAF) and extract informative SNPs
#filtering for ancestry-informative markers — specifically SNPs with:
# ∣MAFref1 − MAFref2∣ > 0.3 
./Allele_Freq0.3.py
#1338 SNPs kept stringent
#2719 SNPs kept relaxed
./Allele_Freq0.2.py
#4919 SNPs retained (realxed)

#Final ELAI dataset
plink --bfile SNP_v3_ELAIprep --extract ELAI_informative.snplist \
 --make-bed --out SNP_ELAI_final
##Stringent:
#5725 variants and 351 people loaded
#1338 variants and 351 people pass filters and QC.
##Relaxed:
#17707 variants and 351 people loaded
#2719 variants and 351 people pass filters and QC.
##MAF difference 0.2:
#17707 variants loaded
#4919 variants and 351 people pass filters and QC.


# Create subdatasets for each population
for pop in $(awk '{print $1}' SNP_ELAI_final.fam | sort | uniq); do
    plink --bfile SNP_ELAI_final --keep <(awk -v p=$pop '$1 == p {print $1, $2}' SNP_ELAI_final.fam) \
          --make-bed --out ${pop}_subset 
done

#Convert plink files to ELAI format:
#Need to change all alllele to single character first for bimbam
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

#Run ELAI for these poplationsusing M. gallo Atlantic as source population 1 and Cromarty (M. edulis) as source population 2
 
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
 
 
#This produced ELAI output files in an output folder with extensions:
# .snpinfo.txt, .admix.txt, .ps21.txt and also .em.txt and .log.txt which I did not use.

./SNP_DensityGenomePlot.py
cd output/
../elai_admix_viz.py
../elai_visualisation.py
../ELAI_EdulisLinePlot.py
