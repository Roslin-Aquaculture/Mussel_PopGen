#!/bin/python
#Name/date: TimRegan/2025.01.10
#File: Summary_Plots.py
 
# Script Breakdown
##Heterozygosity Analysis: Processes all .het files. Calculates observed and expected heterozygosity. Creates a summary CSV and a bar plot comparing observed vs. expected heterozygosity
##Effective Population Size (Ne): Reads .ld files for each population. Computes Ne. Summarizes and plots the average Ne for each population.
##LD Decay: Bins SNP distances (10kb-1Mb), creates R2 for each bin and summary CSV with a line plot showing LD decay for each population.
##Final Summary: Combines heterozygosity and Ne data into a single table. Saves it as final_summary.csv.
##Output Files:
##CSV Files: heterozygosity_summary.csv ne_summary.csv ld_decay_summary.csv final_summary.csv
##Plots: heterozygosity_plot.png ne_plot.png ld_decay_plot.png
 
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob
 
# Parameters
recombination_rate = 2 / 100  # 2 cM/Mb converted to proportion
ld_distance_bins = np.linspace(0.01, 1, 20)  # LD decay bins (0.01 Mb to 1 Mb)
 
# Initialize summary dataframes
heterozygosity_summary = []
ne_summary = []
ld_decay_summary = []
pairwise_fst = []
 
# Process each population
#populations = [f.split('.')[0] for f in glob.glob("*.het")]  # Extract population names
populations = [f.rsplit('.', 1)[0] for f in glob.glob("*.het")]  # Extract population names by removing the last '.ext'
 
# Get all population pairs
population_pairs = [(pop1, pop2) for i, pop1 in enumerate(populations) for pop2 in populations[i+1:]]
 
# Process each Fst file
for pop1, pop2 in population_pairs:
    fst_file = f"fst_{pop1}_{pop2}.fst"
     
    try:
        # Load the Fst summary file
        fst_data = pd.read_csv(fst_file, delim_whitespace=True)
 
        # Calculate the mean Fst value
        mean_fst = fst_data['FST'].mean()  # Adjust column name if needed
 
        # Append result
        pairwise_fst.append({"Pop1": pop1, "Pop2": pop2, "Mean_Fst": mean_fst})
 
    except FileNotFoundError:
        print(f"Warning: {fst_file} not found.")
    except Exception as e:
        print(f"Error processing {fst_file}: {e}")
 
# Convert to DataFrame
fst_df = pd.DataFrame(pairwise_fst)
 
# Make the matrix symmetric by adding reversed pairs
reversed_pairs = fst_df.rename(columns={"Pop1": "Pop2", "Pop2": "Pop1"})
fst_df = pd.concat([fst_df, reversed_pairs])
 
# Pivot to create an Fst matrix
fst_matrix = fst_df.pivot(index="Pop1", columns="Pop2", values="Mean_Fst")
 
# Make matrix symmetric
fst_matrix = fst_matrix.combine_first(fst_matrix.T)
 
# Fill diagonal with zeros (Fst within same population = 0)
np.fill_diagonal(fst_matrix.values, 0)
 
# Save Fst matrix
fst_matrix.to_csv("fst_heatmap_matrix.csv")
 
# Create a heatmap of pairwise Fst values
plt.figure(figsize=(10, 8))
sns.heatmap(fst_matrix, annot=True, cmap="coolwarm", fmt=".3f", cbar_kws={"label": "Mean Fst"})
plt.title("Pairwise Fst Heatmap")
plt.tight_layout()
plt.savefig("fst_heatmap.png")
plt.show()
 
for pop in populations:
    # Load heterozygosity data
    het_file = f"{pop}.het"
    het_data = pd.read_csv(het_file, delim_whitespace=True)
    obs_het = 1 - (het_data['O(HOM)'].sum() / het_data['N(NM)'].sum())
    exp_het = 1 - (het_data['E(HOM)'].sum() / het_data['N(NM)'].sum())
    F_IS = (exp_het - obs_het) / exp_het  # Wright’s F_IS = (He - Ho) / He
    heterozygosity_summary.append({
        'Population': pop,
        'Obs_Het': obs_het,
        'Exp_Het': exp_het,
        'F_IS': F_IS
    })
 
    # Load LD data
    ld_file = f"{pop}.ld"
    ld_data = pd.read_csv(ld_file, delim_whitespace=True)
    ld_data['Distance'] = abs(ld_data['BP_B'] - ld_data['BP_A']) / 1e6  # Convert to Mb
    ld_data = ld_data[(ld_data['Distance'] > 0.01) & (ld_data['Distance'] < 1)]  # Filter distances
    ld_data['Ne'] = 1 / (4 * recombination_rate * ld_data['R2'])
    mean_ne = ld_data['Ne'].mean()
    ne_summary.append({'Population': pop, 'Mean_Ne': mean_ne})
 
    # Summarize LD decay
    ld_data['Bin'] = pd.cut(ld_data['Distance'], ld_distance_bins)
    decay_summary = ld_data.groupby('Bin')['R2'].mean().reset_index()
    decay_summary['Population'] = pop
    ld_decay_summary.append(decay_summary)
 
# Convert summaries to dataframes
heterozygosity_df = pd.DataFrame(heterozygosity_summary)
ne_df = pd.DataFrame(ne_summary)
ld_decay_df = pd.concat(ld_decay_summary)
 
# Save summaries as CSV
heterozygosity_df.to_csv("heterozygosity_summary.csv", index=False)
ne_df.to_csv("ne_summary.csv", index=False)
ld_decay_df.to_csv("ld_decay_summary.csv", index=False)
 
# Merge heterozygosity and Ne data for a final summary
final_summary = pd.merge(heterozygosity_df, ne_df, on='Population')
final_summary.to_csv("final_summary.csv", index=False)
 
# Print final summary
print("Final Summary:")
print(final_summary)
 
# Plot 1: Observed vs Expected Heterozygosity
plt.figure(figsize=(10, 6))
sns.barplot(data=heterozygosity_df.melt(id_vars='Population', var_name='Metric', value_name='Heterozygosity'),
            x='Population', y='Heterozygosity', hue='Metric')
plt.xticks(rotation=45)
plt.title("Heterozygosity Across Populations")
plt.tight_layout()
plt.savefig("heterozygosity_plot.png")
plt.show()

# Plot 1b: F_IS per population
plt.figure(figsize=(10, 6))
sns.barplot(data=heterozygosity_df, x='Population', y='F_IS', color='gray')
plt.xticks(rotation=45)
plt.ylabel("F_IS (Inbreeding Coefficient)")
plt.title("Wright’s F_IS Across Populations")
plt.tight_layout()
plt.savefig("FIS_plot.png")
plt.show()

# Plot 2: Effective Population Size (Ne)
plt.figure(figsize=(10, 6))
sns.barplot(data=ne_df, x='Population', y='Mean_Ne')
plt.xticks(rotation=45)
plt.title("Effective Population Size (Ne) Across Populations")
plt.ylabel("Mean Ne")
plt.tight_layout()
plt.savefig("ne_plot.png")
plt.show()
 
 
# Convert 'Bin' to numeric (using midpoints)
ld_decay_df['Bin'] = ld_decay_df['Bin'].apply(lambda interval: interval.mid)
 
# Plot with seaborn
plt.figure(figsize=(12, 8))
sns.lineplot(data=ld_decay_df, x='Bin', y='R2', hue='Population', marker='o')
 
# Set X-axis limits based on actual data range (or adjust manually if needed)
plt.xlim(ld_decay_df['Bin'].min(), ld_decay_df['Bin'].max())
 
# Update X-axis ticks to show appropriate bin midpoints
plt.xticks(
    ticks=ld_decay_df['Bin'].unique(),
    labels=[f"{b:.2f}" for b in ld_decay_df['Bin'].unique()],
    rotation=45
)
 
# Update axis labels and title
plt.xlabel("Distance (Mb)")
plt.ylabel("LD (R^2)")
plt.title("LD Decay Across Populations")
 
# Adjust legend
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Population")
 
# Tight layout for better spacing
plt.tight_layout()
 
# Save and show the plot
plt.savefig("ld_decay_plot.png")
plt.show()
