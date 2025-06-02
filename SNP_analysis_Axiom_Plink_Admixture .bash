#AXIOM

#Axiom is a genotyping solution that uses microarray technology to analyze genetic variations across the genome. It is widely used in various fields, including agricultural genomics, human genetics, and environmental genomics.

#Use axiom analysis suite software to process the raw data (.cel files) Perform quality control (QC) steps to filter out low-quality dta and ensure accuracy. This includes checks for call rates. 
#I carried out the first QC selection on Axiom. ie: I check for the call rate cut off  ≥ 80 (following Clemence recommendation) (The call rate is a quality metric that measures the proportion of successfully genotyped SNPs or samples. A call rate cutoff of ≥ 80% means that only SNPs or samples with a call rate of 80% or higher are considered acceptable for further analysis) 

		# Axiom SNP array analysis

#Analysis configration = 	Axiom_Myt_v1_96orMore.r1 (Default) 

#Sample QC =  DQC > 0
			#QC call rate > 80
			#Percent of passing sample > 80
			#Avererage call rate for passing > 80
#SNP QC = 	cr-cutoff > 80 
#"cr-cutoff: ≥ 80" refers to the call rate cutoff. The call rate is a quality metric that measures the proportion of successfully genotyped SNPs or samples. A call rate cutoff of ≥ 80% means that only SNPs or samples with a call rate of 80% or higher are considered acceptable for further analysis.


		# PLINK 

# Once the files have been imported to the HPC, QC can be done using PLINK.
#Data missingess within individuals and variants can be explored using the command:

plink --bfile Filename --missing --allow-extra-chr --no-fid --no-parents --no-sex --no-pheno -- hwe 0.000001 --out Filename --make-bed

# The PLINK command processes a binary genotype file set (Filename) by allowing the inclusion of non-standard chromosomes and generating output files with the same prefix (Filename). It calculates missing genotype statistics (--missing) and ignores family IDs, parental IDs, sex, and phenotype information (--no-fid, --no-parents, --no-sex, --no-pheno). Finally, it creates a new binary file set (--make-bed). This command simplifies the dataset by excluding unnecessary metadata fields and computes missing data statistics, preparing the data for subsequent genetic analysis.
#Missing Data Statistics: It calculates missing genotype statistics and creates .lmiss and .imiss files.
#Handling Extra Chromosomes: It allows the processing of non-standard chromosomes.
#Ignoring Specific Fields: It ignores family ID, parental IDs, sex, and phenotype information in the input data.
#Output Files: It generates a new set of binary PLINK files (newFilename.bed, newFilename.bim, newFilename.fam) and the missing data statistics files (newFilename.lmiss, newFilename.imiss).


plink --bfile Filename --allow-extra-chr --out Filename --mind 0.15 --geno 0.1 --maf 0.01 --freq --make-bed

# The PLINK command processes a binary genotype file set (Filename) by allowing the inclusion of non-standard chromosomes and generating output files with the prefix testQC. It filters the data to exclude individuals with more than 15% missing genotypes (--mind 0.15), SNPs with more than 5% missing data (--geno 0.05), and SNPs with a minor allele frequency below 2% (--maf 0.02). Additionally, it calculates allele frequencies (--freq) and creates a new binary file set (--make-bed). This preprocessing ensures high-quality data for downstream genetic analysis by retaining only reliable SNPs and samples.

plink --bfile Filename --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out Filename --make-bed

#The PLINK command processes a binary genotype file set (Filename) by allowing the inclusion of non-standard chromosomes and generating output files with the same prefix (Filename). It sets both individual IDs and family IDs to the same value (--double-id), assigns missing variant IDs using a specified format (--set-missing-var-ids @:#), and performs linkage disequilibrium (LD) pruning with a sliding window approach (--indep-pairwise 50 10 0.1). Finally, it creates a new binary file set (--make-bed). This command prepares the dataset by ensuring proper ID formatting, handling missing variant IDs, and pruning SNPs in high LD, which is useful for downstream analyses requiring independent SNPs.

plink --bfile Filename --double-id --allow-extra-chr --set-missing-var-ids @:# --extract Filename.prune.in --make-bed --pca --out Filename

#The PLINK command processes a binary genotype file set (Filename) by allowing the inclusion of non-standard chromosomes and generating output files with the same prefix (Filename). It sets both individual IDs and family IDs to the same value (--double-id), assigns missing variant IDs using a specified format (--set-missing-var-ids @:#), and extracts a pruned set of SNPs listed in outputfile.prune.in (--extract). It then performs principal component analysis (PCA) on the pruned dataset (--pca) and creates a new binary file set (--make-bed). This command prepares the dataset for PCA by ensuring proper ID formatting, handling missing variant IDs, extracting independent SNPs, and performing PCA to identify population structure or genetic variation patterns.

		#ADMIXTURE 


awk '{$1=0;print $0}' (Filename).bim > (Filename).bim.tmp
mv (Filename).bim.tmp (Filename).bim

admixture --cv (Filename).bed 2 > log2.out

# Two files will be produced
# The .Q file contains the cluster assignments for each individual
# The .P file contains the population allele frequencies for each SNP included in the analysis

# Contrary to Structure analysis, admixutre calculates the best number of K clusters representing the population being analysed, with no prior assumption. 
# Therefore, in this example, a set of k clusters, generally from 2 to 10, will be analysed in a loop and the output directed to a log.out file. However, more this value can be adjusted based on preivous knowledge of the population set being analysed in the study. must test the best fitting K value for the population being analysed.  

for i in {2..10}
do
admixture --cv (Filename).bed $i > log${i}.out
done

# The best fitting number of K clusters will be the one with the lowest number of cross-validation error
# This value can be extracted from the log.out generated files using the following command:

awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > (Filename).cv.error


#to get the sample id list 
awk '{print $1}' "Filename.ped" > Filename.txt