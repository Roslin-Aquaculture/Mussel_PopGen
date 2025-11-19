# PCA and ADMIXTURE Pipeline

This module performs LD-pruned PCA and ADMIXTURE analyses for each dataset version (v1–v3).
A complete example is shown for **v1**, and identical commands must be run for **v2** and **v3**.

---

# 1. Load software

```bash
module load igmm/apps/plink/1.90b7.2
module load roslin/admixture/1.3.0
```

# 2. LD-pruning (example: v1)
```bash
plink --bfile SNP_v1_110825 --double-id --allow-extra-chr \
 --set-missing-var-ids @:# \
 --indep-pairwise 50 10 0.1 \
 --out SNP_v1_110825
```

# 3. Extract pruned SNPs and compute PCA
```bash
plink --bfile SNP_v1_110825 --double-id --allow-extra-chr \
 --extract SNP_v1_110825.prune.in \
 --make-bed --pca \
 --out SNP_v1_110825
```

# Outputs:
SNP_v1_110825.eigenvec
SNP_v1_110825.eigenval

# 4. ADMIXTURE (example: v1)
```bash
for K in {2..6}; do
    admixture --cv SNP_v1_110825.bed $K > log${K}.out
done
```

# Extract CV error:
```bash
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > SNP_v1_110825.cv.error
```
# 5. Label preparation
```bash
cut -d' ' -f2 SNP_v1_110825.fam > SNP_v1_110825.txt
```

# 6. Repeat for v2 and v3
Use exactly the same commands, substituting:
SNP_v2_110825
SNP_v3_110825
