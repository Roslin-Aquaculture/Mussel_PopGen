Includes *fixed* convergence diagnostic script and triangle plot script.

# entropy Analysis

This module runs entropy on the **v3 dataset only** and provides:

- MCMC setup  
- Convergence assessment  
- Q-matrix extraction  
- Triangle plot for hybridisation interpretation  

---

# Running entropy
Load modules, activate conda env and navigate to data
```bash
module load igmm/apps/plink/1.90b7.2
module load roslin/admixture/1.3.0
module load anaconda
module load R
conda activate elai_env

cd SNP_110825/SNP_v3_110825
```

# Create merge list
```bash
ls *_clean.bed | sed 's/.bed//' | grep -v Aberdeen_clean > mergelist.txt
```

# Merge into one dataset
```bash
plink --bfile Aberdeen_clean --merge-list mergelist.txt --make-bed --out mussel_all
```
# 1. Use plink to get a VCF, then convert with a simple Python script.
```bash
plink --bfile mussel_all --recode vcf --out mussel_all
```

# Needed to get R scripts from source:
```bash
wget https://bitbucket.org/buerklelab/mixedploidy-entropy/get/246ccf1003c4.zip
unzip 246ccf1003c4.zip
mv -r 246ccf1003c4 entropy
```
# Install further required packages:
```bash
conda install popgen-entropy
conda install -c conda-forge r-vcfr
conda install conda-forge::hdf5
conda install bioconda::bioconductor-rhdf5
```

# Create required ploidy filoe from vcf
```bash
grep -m 1 "^#CHROM" mussel_all.vcf | cut -f10- | tr '\t' '\n' | sed 's/.*/2/' > ploidy_inds.txt
```

# Custom R script
run from the directory that contains inputdataformat.R (from the entropy distribution)
```bash
Rscript entropy/auxfiles/inputdataformat.R mussel_all.vcf mussel_all 2
```

# 4. Run entropy for each K value, Where:
Input data file (-i): mpgl (or another accepted format). This is your genotype/read-count matrix in the format entropy expects.

Ploidy (-n): either a single value (e.g. 2) or a text file listing ploidy per individual/locus as required.

Optional: an initial admixture proportions file (-q) can be supplied to initialise q values.

MCMC settings: -k (clusters), -l (total MCMC steps), -b (burn-in),

  -t (thin/store every t steps), -o (output .hdf5).
  
 -D 1 requests DIC/WAIC if compiled into the binary; this helps compare K.
 
```bash
cd /exports/cmvm/eddie/eb/groups/Regan_grp/SNP_110825/SNP_v3_110825
```

# full runs (longer) after you validate diagnostics above — recommended
adjust -l (iterations), -b (burn-in), -t (thin) as desired.
```bash
for K in 1 2 3 4; do
  for CH in 1 2 3; do
    entropy -i mussel_all.mpgl -m 1 -n 2 -k ${K} \
      -q qk${K}inds.txt \
      -l 5000 -b 2000 -t 10 -D 1 \
      -o Full_entropy/mussel_k${K}_chain${CH}.hdf5 \
      -r $RANDOM
  done
done

cd Full_entropy/
```

# Measure convergence with this custom R script:
produces converge_K2_convergence.txt summary files
```bash
Rscript convergence_diag_fixed.R converge_K2 2 \
mussel_k2_chain1.hdf5 mussel_k2_chain2.hdf5 mussel_k2_chain3.hdf5

Rscript convergence_diag_fixed.R converge_K3 3 \
mussel_k3_chain1.hdf5 mussel_k3_chain2.hdf5 mussel_k3_chain3.hdf5

Rscript convergence_diag_fixed.R converge_K4 4 \
mussel_k4_chain1.hdf5 mussel_k4_chain2.hdf5 mussel_k4_chain3.hdf5
```
# Example output:
K2

Potential scale reduction factors:

     Point est. Upper C.I.
     
[1,]       1.01       1.01

K3

Potential scale reduction factors:

     Point est. Upper C.I.
     
[1,]       1.15       1.17

K4

Potential scale reduction factors:

     Point est. Upper C.I.
     
[1,]        1.2       1.24


# Plot the admixture, custom R script
e.g. for K=2
```bash
Rscript plot_admix_from_hdf5.R inds_mussel_all.txt K2 2 mussel_k2_chain1.hdf5 mussel_k2_chain2.hdf5 mussel_k2_chain3.hdf5
```

# Triangle plot custom R script:
```bash
Rscript triangle_plot_entropy.R triangle_k2_mussel.pdf \
  inds_mussel_all.txt \
  mussel_k2_chain1.hdf5 mussel_k2_chain2.hdf5 mussel_k2_chain3.hdf5
```
