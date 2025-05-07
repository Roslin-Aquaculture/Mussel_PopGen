"""
Ensure that all relevant Python scripts have been added to the directory and updated accordingly (see notes)
Update_PLINK_Map.py
Summary_Plots.py
elai_admix_viz.py
ELAI_EdulisLinePlot.py
"""


#2025/04/22: Wild population from Aberdeen added.
module load igmm/apps/vcftools/0.1.16
module load igmm/apps/plink/1.90b7.2
module load anaconda
conda activate orthofinder
#THis conda env has the appropriate viz packages I need, orthofinder is not actually necessary though

#Move to dir with PLINK output files from Ambre's analyses using Axiom
cd /exports/cmvm/eddie/eb/groups/Regan_grp/ELAI/SNP_v3_100225
 
#Ensure TOL_Med_SNP Pos.txt is present. Otherwise same as previous script, 
#but changed to ensure  SNP_v3_260325 was being examined.
./Update_PLINK_Map.py
 
cp SNP_v3_260325_updated.map SNP_v3_260325.map
#Still needed to replace all instances of "---" with "0" for position or else I got an error. 

plink --file SNP_v3_260325 --missing --allow-extra-chr --no-fid --no-parents \
--no-sex --no-pheno --hwe 0.000001 --geno 0.1 --out SNP_v3_260325 --make-bed
# 7468 variants and 358 mussels passed.

plink --bfile SNP_v3_260325  --allow-extra-chr --mind 0.15  --maf 0.01 --freq \
 --set-missing-var-ids @:# --out SNP_v3_260325 --make-bed
# 3063 variants removed due to MAF, 4405 variants and 358 mussels left
 
#Update population_map.txt
#Verify that our .fam file has the FID column set to population labels:
awk 'NR==FNR {pop[$1]=$2; next} {sub(/\r$/, "", pop[$2]); printf "%s %s %s %s %s %s\n", pop[$2], $2, $3, $4, $5, $6}' population_map.txt SNP_v3_260325.fam > updated.fam
  
#Replace old .fam
mv updated.fam SNP_v3_260325.fam
 
#Calculate pairwise Fst using PLINK
for pop1 in $(awk '{print $1}' SNP_v3_260325.fam | sort -u); do
    for pop2 in $(awk '{print $1}' SNP_v3_260325.fam | sort -u); do
        if [ "$pop1" != "$pop2" ]; then
            echo "Calculating Fst between $pop1 and $pop2"
            # Generate a properly formatted within file (FID IID GroupID)
            awk -v p1="$pop1" -v p2="$pop2" '
                $1 == p1 {print $1, $2, 1}
                $1 == p2 {print $1, $2, 2}
            ' SNP_v3_260325.fam > temp_within.txt
            # Debug: Check the first few lines
            head temp_within.txt
            # Run PLINK Fst
            plink --bfile SNP_v3_260325 --fst --within temp_within.txt --allow-extra-chr --out fst_${pop1}_${pop2}
            # Clean up
            rm temp_within.txt
        fi
    done
done
 
# Create subdatasets for each population
for pop in $(awk '{print $1}' SNP_v3_260325.fam | sort | uniq); do
    plink --bfile SNP_v3_260325 --keep <(awk -v p=$pop '$1 == p {print $1, $2}' SNP_v3_260325.fam) \
          --make-bed --out ${pop}_subset --allow-extra-chr
done
 
#Calculate O(HOM) (Observed homozygous genotypes), E(HOM) (Expected homozygous genotypes) and F (Inbreeding coefficient)
for pop in $(awk '{print $1}' SNP_v3_260325.fam | sort | uniq); do
  echo "Processing population: $pop"
  plink --bfile "$pop"_subset --allow-extra-chr --het --out $pop
done
  
#Calculate Effective Population Size (Ne)
for pop in $(awk '{print $1}' SNP_v3_260325.fam | sort | uniq); do
  echo "Processing population: $pop"
  plink --bfile "$pop"_subset --allow-extra-chr --r2 gz --ld-window-kb 1000 --out $pop
done
  
#Calculate Linkage Disequilibrium (LD) Decay
for pop in $(awk '{print $1}' SNP_v3_260325.fam | sort | uniq); do
  echo "Processing population: $pop"
  plink --bfile "$pop"_subset --allow-extra-chr --r2 --ld-window-r2 0 --ld-window 99999 --ld-window-kb 1000 --out $pop
  awk '{print $2, $5, $7}' $pop.ld > $pop.ld_summary
done
 
#Calculate Fst: Creates clusters according to Family ID:
plink --bfile SNP_v3_260325 --fst --family --allow-extra-chr --out fst_results
 
#Updated this script, see below.
./Summary_Plots.py

###########################
## Moving onto ELAI
conda deactivate
conda activate elai_env

 
#Convert plink files to ELAI format:
plink --bfile M.gallo_Atl_subset --recode bimbam --allow-extra-chr 0 --out M.gallo_Atl_elai
plink --bfile Cromarty_subset --recode bimbam --allow-extra-chr 0 --out Cromarty_elai
plink --bfile Shetland_subset --recode bimbam --allow-extra-chr 0 --out Shetland_elai
plink --bfile Western_Isles_subset --recode bimbam --allow-extra-chr 0 --out Western_Isles_elai
plink --bfile Irl_Atl_subset --recode bimbam --allow-extra-chr 0 --out Irl_Atl_elai
plink --bfile M.edulis_Nor_subset --recode bimbam --allow-extra-chr 0 --out M.edulis_Nor_elai
plink --bfile M.edulis_Eng_subset --recode bimbam --allow-extra-chr 0 --out M.edulis_Eng_elai
plink --bfile Aberdeen_subset --recode bimbam --allow-extra-chr 0 --out Aberdeen_elai

#Run ELAI for these poplationsusing M. gallo Atlantic as source population 1 and Cromarty (M. edulis) as source population 2
 
../../../ELAI/elai-lin -g M.gallo_Atl_elai.recode.geno.txt -p 10 -g Cromarty_elai.recode.geno.txt -p 11 \
-g Western_Isles_elai.recode.geno.txt -p 1 \
-pos Western_Isles_elai.recode.pos.txt \
-C 2 -c 10 -s 30 -mg 50 -o ELAI_Results_WesIsl
 
../../../ELAI/elai-lin -g M.gallo_Atl_elai.recode.geno.txt -p 10 -g Cromarty_elai.recode.geno.txt -p 11 \
-g Shetland_elai.recode.geno.txt -p 1 \
-pos Shetland_elai.recode.pos.txt \
-C 2 -c 10 -s 30 -mg 50 -o ELAI_Results_Shetland
 
../../../ELAI/elai-lin -g M.gallo_Atl_elai.recode.geno.txt -p 10 -g Cromarty_elai.recode.geno.txt -p 11 \
-g Irl_Atl_elai.recode.geno.txt -p 1 \
-pos Irl_Atl_elai.recode.pos.txt \
-C 2 -c 10 -s 30 -mg 50 -o ELAI_Results_Irl_Atl
 
../../../ELAI/elai-lin -g M.gallo_Atl_elai.recode.geno.txt -p 10 -g Cromarty_elai.recode.geno.txt -p 11 \
-g M.edulis_Eng_elai.recode.geno.txt -p 1 \
-pos M.edulis_Eng_elai.recode.pos.txt \
-C 2 -c 10 -s 30 -mg 50 -o ELAI_Results_M.edulis_Eng
 
../../../ELAI/elai-lin -g M.gallo_Atl_elai.recode.geno.txt -p 10 -g Cromarty_elai.recode.geno.txt -p 11 \
-g M.edulis_Nor_elai.recode.geno.txt -p 1 \
-pos M.edulis_Nor_elai.recode.pos.txt \
-C 2 -c 10 -s 30 -mg 50 -o ELAI_Results_M.edulis_Nor
 
../../../ELAI/elai-lin -g M.gallo_Atl_elai.recode.geno.txt -p 10 -g Cromarty_elai.recode.geno.txt -p 11 \
-g Aberdeen_elai.recode.geno.txt -p 1 \
-pos Aberdeen_elai.recode.pos.txt \
-C 2 -c 10 -s 30 -mg 50 -o ELAI_Results_Aberdeen

#This produced ELAI output files in an output folder with extensions:
# .snpinfo.txt, .admix.txt, .ps21.txt and also .em.txt and .log.txt which I did not use.

cd output/
../elai_admix_viz.py
../ELAI_EdulisLinePlot.py
