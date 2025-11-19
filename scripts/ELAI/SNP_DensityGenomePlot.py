#!/bin/python
# -*- coding: utf-8 -*-
# Name/date: TimRegan/2025.05.06
# File: SNP_DensityGenomePlot.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def load_snp_info(snpinfo_file, has_header=True, column_order=("SNP_ID", "Chr", "Pos")):
    """
    Load and clean SNP info file. Can handle both headered and headerless input.
    """
    sep = "\t" if has_header else r"\s+"
    header = 0 if has_header else None
    snpinfo = pd.read_csv(snpinfo_file, sep=sep, header=header, dtype=str)

    if has_header:
        snpinfo.columns = [col.strip() for col in snpinfo.columns]
        snpinfo = snpinfo.rename(columns={
            "Marker ID": "SNP_ID",
            "Chromosome": "Chr",
            "Physical Position": "Pos"
        })
    else:
        snpinfo.columns = column_order

    snpinfo["Pos"] = pd.to_numeric(snpinfo["Pos"], errors="coerce")
    snpinfo["Chr"] = pd.to_numeric(snpinfo["Chr"], errors="coerce")
    snpinfo = snpinfo.dropna(subset=["Chr", "Pos"])
    snpinfo["Chr"] = snpinfo["Chr"].astype(int)

    return snpinfo

def compute_genome_bins(snpinfo, bin_size=1_000_000):
    """
    Bin SNPs into 1 Mb windows across genome, return dataframe with bin counts and positions.
    """
    snpinfo = snpinfo[~snpinfo["Chr"].isin(["0", 0])]
    snpinfo["Chr"] = pd.to_numeric(snpinfo["Chr"], errors="coerce")
    snpinfo = snpinfo.dropna(subset=["Chr"])
    snpinfo["Chr"] = snpinfo["Chr"].astype(int)
    snpinfo = snpinfo.sort_values(["Chr", "Pos"])

    bins = []
    chr_offsets = {}
    cumulative = 0

    for chrom in sorted(snpinfo["Chr"].unique()):
        chr_snps = snpinfo[snpinfo["Chr"] == chrom]
        chr_min = chr_snps["Pos"].min()
        chr_max = chr_snps["Pos"].max()
        chr_offsets[chrom] = cumulative
        chrom_length = chr_max - chr_min + 1
        num_bins = int(np.ceil(chrom_length / bin_size))

        for i in range(num_bins):
            start = chr_min + i * bin_size
            end = min(start + bin_size, chr_max + 1)
            count = chr_snps[(chr_snps["Pos"] >= start) & (chr_snps["Pos"] < end)].shape[0]
            genome_pos = cumulative + start
            bins.append((chrom, start, end, count, genome_pos))

        cumulative += chr_max - chr_min + 1

    return pd.DataFrame(bins, columns=["Chr", "Start", "End", "SNP_Count", "Genome_Pos"]), chr_offsets

def plot_snp_density(density_df, chr_offsets, title, output_file):
    """
    Plot SNP density heatmap across genome using 1 Mb bins.
    """
    plt.figure(figsize=(16, 2))  # Reduced height
    fig, ax = plt.subplots(figsize=(16, 2))

    cmap = plt.cm.viridis
    norm = mcolors.Normalize(vmin=0, vmax=density_df["SNP_Count"].max())

    for _, row in density_df.iterrows():
        ax.axvspan(row["Genome_Pos"], row["Genome_Pos"] + 1_000_000,
                   color=cmap(norm(row["SNP_Count"])), linewidth=0)

    # Chromosome boundaries
    for chrom, offset in chr_offsets.items():
        ax.axvline(x=offset, color='grey', linestyle='--', linewidth=0.5)
        chr_mid = offset + (density_df[density_df["Chr"] == chrom]["Start"].max() / 2)
        ax.text(chr_mid, 1.02, f"Chr{chrom}", ha='center', va='bottom', fontsize=8, transform=ax.get_xaxis_transform())

    ax.set_xlim(density_df["Genome_Pos"].min(), density_df["Genome_Pos"].max())
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title(f"{title} - SNP Density Heatmap (1 Mb bins)", fontsize=11, pad=20)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation="vertical", pad=0.01)
    cbar.set_label("SNP Count per 1Mb", rotation=270, labelpad=15)

    plt.tight_layout()
    plt.savefig(f"{output_file}_snp_density.png", dpi=300, bbox_inches="tight")
    plt.show()

def main():
    full_snp_file = "TOL_Med_SNP Pos.txt"
    filtered_snp_file = "Ireland_elai.recode.pos.txt"
    title_base = "xbMytEdul2.2"
    output_prefix = "TOL_SNP_Density"

    # Load full SNPs (has header)
    print("Loading full SNP list...")
    full_snps = load_snp_info(full_snp_file, has_header=True)
    full_density_df, full_chr_offsets = compute_genome_bins(full_snps)
    plot_snp_density(full_density_df, full_chr_offsets,
                     f"{title_base} - Full SNPs", f"{output_prefix}_full")

    # Load filtered SNPs (no header, space-separated, reordered columns)
    print("Loading filtered SNP list...")
    filtered_snps = load_snp_info(
        filtered_snp_file,
        has_header=False,
        column_order=("SNP_ID", "Pos", "Chr")
    )
    filtered_density_df, filtered_chr_offsets = compute_genome_bins(filtered_snps)
    plot_snp_density(filtered_density_df, filtered_chr_offsets,
                     f"{title_base} - Filtered SNPs", f"{output_prefix}_filtered")



if __name__ == "__main__":
    main()
