# Mussel_PopGen
This repository provides a complete, reproducible workflow for analysing population structure, introgression, and local ancestry in *Mytilus* spp. using the 60k MytiSNP-array  ([See Nascimento-Schulze et al.](https://doi.org/10.1111/eva.13552)).

This workflow includes QC, PCA, ADMIXTURE, FST, LD, entropy, ELAI, and a triangle plot for admixture of K=2, together with documentation and example commands.

Here, we are workkng on datasets which have already been produced from raw .cel files using AxiomSuite.   

## Dataset Versions

Three filtered datasets are used throughout this study:

| Version | Description | Typical Use |
|--------|-------------|--------------|
| **v1** | Full dataset (all individuals, including *M. trossulus* and *M. gallo mediterranean*) | PCA, ADMIXTURE, initial exploratory analyses |
| **v2** | *M. trossulus* removed | PCA, ADMIXTURE, FST, LD |
| **v3** | *M. trossulus* + *M. gallo mediterranean* removed | **Primary dataset** for entropy, ELAI, introgression, LD, FST |

Each module includes a full worked example on **v1**, and you must repeat the same pipeline for **v2** and **v3**.

---

# Modules

## 1. QC and Filtering
See: `scripts/QC/QC_README.md`  
Includes SNP call filtering, missingness thresholds, removal of sex-linked markers where applicable, and generation of v1–v3 datasets.

## 2. PCA and ADMIXTURE Population Structure
See: `scripts/PCA_ADMIXTURE/PCA_ADMIXTURE.md`  
Includes full worked example for v1 and instructions to repeat for v2 and v3.

## 3. Entropy (Structure-like Bayesian Clustering)
See: `scripts/entropy/README_entropy.md`  
Includes:  
- Running entropy on v3  
- Convergence diagnostics using deviance  
- Extracting Q matrices  
- Plotting triangle plots (updated robust script included)

## 4. ELAI (Local Ancestry)
See: `scripts/ELAI/README_ELAI.md`  
Includes:  
- Input preparation  
- Reference definition  
- Running ELAI 
- Visualisation
- 

---

# Requirements

- PLINK ≥ 1.9  
- ADMIXTURE ≥ 1.3  
- entropy (v1.3 or later)  
- ELAI (latest commit)  
- R ≥ 4.0 with packages: `hdf5r`, `ggplot2`, `dplyr`, `gridExtra`, `reshape2`, `RColorBrewer`

---

# Citation
If using this repository, please cite the accompanying publication and the tools used: PLINK, ADMIXTURE, entropy, and ELAI.

# Extra notes for working with other datasets
Ensure that all of the python scripts are in your working directory.
Ensure all conda packages are installed on whatever environment you are using (e.g. elai_env).
Change all appropriate lines in the custom python and R scripts.
The reference file "TOL_Med_SNP Pos.txt" was generated to position all SNPs along the TOL *M. edulis* assembly (i.e. GCF_963676685.1). Make sure it is present. 


