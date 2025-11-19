#!/bin/python
# -*- coding: utf-8 -*-
# Name/date: TimRegan/2025.08.07
# File: ELAI_EdulisLinePlot.py

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.patches as patches

pop_colors = {
    #"M.edulis_Nor": "red",
    "Aberdeen": "teal",
    "Ireland": "gold",
    "M.edulis_Eng": "darkslategrey",
    "Shetland_A": "magenta",
    "Shetland_B": "purple",
    "Western_Isles": "lightgreen"
}


def load_elai_data(ps21_file, snpinfo_file, num_individuals, num_snps, num_sources=2):
    """
    Load and reshape ELAI ancestry dosage data.
    """
    # Load SNP info file
    snpinfo = pd.read_csv(snpinfo_file, sep="\t", header=0, dtype=str)
    snpinfo.columns = [col.strip() for col in snpinfo.columns]
    snpinfo = snpinfo[["rs", "pos", "chr"]].copy()
    snpinfo.columns = ["SNP_ID", "Pos", "Chr"]

    snpinfo["Pos"] = pd.to_numeric(snpinfo["Pos"], errors="coerce")
    snpinfo["Chr"] = snpinfo["Chr"].str.strip()
    snpinfo = snpinfo.dropna(subset=["Chr", "Pos"])
    snpinfo["Chr"] = snpinfo["Chr"].astype(str)

    print(f"Filtered SNP info: {snpinfo.shape[0]} SNPs")

    # Load ELAI ancestry data
    ancestry_data = np.loadtxt(ps21_file)
    print(f"Loaded ancestry matrix: {ancestry_data.shape}")

    ancestry_data = ancestry_data.reshape(num_individuals, num_snps, num_sources)
    ancestry_avg = ancestry_data.mean(axis=0)  # average across individuals
    total_dosage = ancestry_avg.sum(axis=1)
    ancestry_prop = ancestry_avg / total_dosage[:, np.newaxis]
    ancestry_df = pd.DataFrame(ancestry_prop, columns=["Ancestry_1", "Ancestry_2"])
    ancestry_df = ancestry_df.iloc[: len(snpinfo)]

    # Merge SNP info with ancestry
    combined_data = pd.concat([snpinfo, ancestry_df], axis=1)
    return combined_data

def plot_edulis_ancestry(data, title, output_file):
    """
    Plot only M. edulis ancestry (Ancestry_2) with alternating colour lines per chromosome,
    and overlay mean +/- standard deviation across all SNPs.
    """
    data = data[~data["Chr"].isin(["0", 0])]
    data["Chr"] = pd.to_numeric(data["Chr"], errors="coerce")
    data = data.dropna(subset=["Chr"])
    data["Chr"] = data["Chr"].astype(int)
    data = data.sort_values(["Chr", "Pos"])

    # Compute genome-wide position offsets
    chr_offsets = {}
    cumulative_length = 0
    chr_ranges = {}  # We will use this to store the start and end genomic positions of each chromosome.
    for chrom in sorted(data["Chr"].unique()):
        chr_data = data[data["Chr"] == chrom]
        min_pos = chr_data["Pos"].min()
        max_pos = chr_data["Pos"].max()
        chr_offsets[chrom] = cumulative_length
        cumulative_length += max_pos - min_pos + 1
        chr_ranges[chrom] = (chr_offsets[chrom], chr_offsets[chrom] + max_pos - min_pos)

    data["Genome_Pos"] = data.apply(lambda row: row["Pos"] + chr_offsets[row["Chr"]], axis=1)

    # Compute mean and SSD across the genome
    mean_dosage = data["Ancestry_2"].mean()
    std_4dosage = 4*(data["Ancestry_2"].std())
    std_6dosage = 6 * data["Ancestry_2"].std()

    # Plot
    fig, ax = plt.subplots(figsize=(16, 3))

   # Alternate background color
    for i, (chrom, (start, end)) in enumerate(chr_ranges.items()):
        color = "#f0f0f0" if i % 2 == 0 else "#ffffff"
        ax.axvspan(start, end, facecolor=color, zorder=0)

    for i, chrom in enumerate(sorted(data["Chr"].unique())):
        chr_data = data[data["Chr"] == chrom]
        color = "blue" if i % 2 == 0 else "black"
        ax.plot(chr_data["Genome_Pos"], chr_data["Ancestry_2"], color=color, linewidth=0.7)

        # Label chromosomes
        mid = chr_offsets[chrom] + (chr_data["Pos"].max() - chr_data["Pos"].min()) / 2
        ax.text(mid, 1.03, f"Chr{chrom}", ha='center', va='bottom', fontsize=10, transform=ax.get_xaxis_transform())

    # Mean and ±SSD lines
    #ax.axhline(y=mean_dosage, color="red", linestyle="-", linewidth=1, label="Mean Edulis Dosage")
    #ax.axhline(y=mean_dosage + std_4dosage, color="grey", linestyle="--", linewidth=0.8, label="+/- 4*SD")
    ax.axhline(y=mean_dosage, color="red", linestyle="-", linewidth=1)
    ax.axhline(y=mean_dosage + std_4dosage, color="grey", linestyle="--", linewidth=0.8)
    ax.axhline(y=mean_dosage - std_4dosage, color="grey", linestyle="--", linewidth=0.8)

    # Axes and labels
    ax.set_xlabel("Genome-wide Position (Gb)")
    ax.set_ylabel("$\it{M. edulis}$ Dosage")
    #ax.set_ylim(0, 1)
    # Set custom y-axis limits based on mean and 6*SD
    ax.set_ylim(max(0, mean_dosage - std_6dosage), min(1, mean_dosage + std_6dosage))

    buffer = (data["Genome_Pos"].max() - data["Genome_Pos"].min()) * 0.02
    ax.set_xlim(data["Genome_Pos"].min() - buffer, data["Genome_Pos"].max() + buffer)

    ax.set_title(f"{title}", fontsize=11, pad=25)
    ax.legend(loc='upper right', frameon=False)

    plt.tight_layout()
    plt.savefig(f"{output_file}", dpi=300, bbox_inches="tight")
    plt.close()

def create_standalone_legend(output_file="ancestry_legend.png"):
    """Create a standalone legend plot."""
    fig, ax = plt.subplots(figsize=(3, 3))
    
    # Create proxy artists for the legend
    legend_elements = [
        patches.Patch(facecolor=pop_colors[pop], label=pop)
        for pop in pop_colors.keys()
    ]
    
    # Create legend
    ax.legend(handles=legend_elements, frameon=False)
    
    # Remove axes
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()

def plot_all_populations(all_data, output_file="All_Populations_Medulis_ancestry.png"):
    """Plot all populations' M. edulis ancestry on the same figure."""
    fig, ax = plt.subplots(figsize=(16, 6))
    
    # First population for chromosome structure
    first_pop_data = next(iter(all_data.values()))
    
    # Compute genome-wide positions (same for all populations)
    data = first_pop_data.copy()
    data = data[~data["Chr"].isin(["0", "0.0", "0"])]
    data["Chr"] = pd.to_numeric(data["Chr"], errors="coerce")
    data = data.dropna(subset=["Chr"])
    data["Chr"] = data["Chr"].astype(int)
    data = data.sort_values(["Chr", "Pos"])

    # Calculate chromosome offsets
    chr_offsets = {}
    cumulative_length = 0
    chr_ranges = {}
    for chrom in sorted(data["Chr"].unique()):
        chr_data = data[data["Chr"] == chrom]
        min_pos = chr_data["Pos"].min()
        max_pos = chr_data["Pos"].max()
        chr_offsets[str(chrom)] = cumulative_length
        cumulative_length += max_pos - min_pos + 1
        chr_ranges[chrom] = (chr_offsets[str(chrom)], chr_offsets[str(chrom)] + max_pos - min_pos)

    # Calculate genome positions for reference data
    data["Genome_Pos"] = data.apply(
        lambda row: row["Pos"] + chr_offsets[str(row["Chr"])], 
        axis=1
    )

    # Alternate background color
    for i, (chrom, (start, end)) in enumerate(chr_ranges.items()):
        color = "#f0f0f0" if i % 2 == 0 else "#ffffff"
        ax.axvspan(start, end, facecolor=color, zorder=0)

    # Plot each population
    for pop_name, pop_data in all_data.items():
        pop_data = pop_data.copy()
        pop_data = pop_data[pop_data["Chr"].isin(chr_offsets.keys())]
        
        pop_data["Genome_Pos"] = pop_data.apply(
            lambda row: row["Pos"] + chr_offsets[str(row["Chr"])], 
            axis=1
        )
        
        ax.plot(pop_data["Genome_Pos"], pop_data["Ancestry_2"], 
                color=pop_colors[pop_name], 
                linewidth=0.7, 
                alpha=0.7)

    # Add chromosome labels
    for chrom in sorted(data["Chr"].unique()):
        chr_data = data[data["Chr"] == chrom]
        mid = chr_offsets[str(chrom)] + (chr_data["Pos"].max() - chr_data["Pos"].min()) / 2
        ax.text(mid, 1.03, f"Chr{chrom}", ha='center', va='bottom', fontsize=10, transform=ax.get_xaxis_transform())

    # Configure plot
    ax.set_xlabel("Genome-wide Position (Gb)")
    ax.set_ylabel("$\it{M. edulis}$ Dosage")
    ax.set_ylim(0, 1)
    
    buffer = (data["Genome_Pos"].max() - data["Genome_Pos"].min()) * 0.02
    ax.set_xlim(data["Genome_Pos"].min() - buffer, data["Genome_Pos"].max() + buffer)

    ax.set_title("All Populations - $\it{M. edulis}$ Ancestry", fontsize=11, pad=25)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    
    # Create separate legend
    create_standalone_legend()

def main():
    num_snps = 6802
    all_population_data = {}  # Store data for combined plot

    populations = {
        "Western_Isles": {"file_prefix": "ELAI_Results_Western_Isles", "n_ind": 34},
        "Shetland_B": {"file_prefix": "ELAI_Results_Shetland_B", "n_ind": 34},
        "Shetland_A": {"file_prefix": "ELAI_Results_Shetland_A", "n_ind": 29},
        #"M.edulis_Nordic": {"file_prefix": "ELAI_Results_M.edulis_Nordic", "n_ind": 36},
        "M.edulis_Eng": {"file_prefix": "ELAI_Results_M.edulis_Eng", "n_ind": 18},
        "Ireland": {"file_prefix": "ELAI_Results_Ireland", "n_ind": 36},
        "Aberdeen": {"file_prefix": "ELAI_Results_Aberdeen", "n_ind": 25},
    }

    for pop_name, pop_info in populations.items():
        ps21_file = f"{pop_info['file_prefix']}.ps21.txt"
        snpinfo_file = f"{pop_info['file_prefix']}.snpinfo.txt"
        data = load_elai_data(ps21_file, snpinfo_file, pop_info["n_ind"], num_snps)

        if data.empty:
            print(f"ERROR: No data loaded for {pop_name}. Skipping...")
            continue

        output_file = f"{pop_name}_Medulis_ancestry_plot.png"
        plot_edulis_ancestry(data, f"{pop_name} - M. edulis Ancestry", output_file)
        
        # Store data for combined plot
        all_population_data[pop_name] = data

    # Create combined plot if we have data
    if all_population_data:
        plot_all_populations(all_population_data)
        create_standalone_legend() 
    else:
        print("No valid population data found for combined plot.")

if __name__ == "__main__":
    main()
