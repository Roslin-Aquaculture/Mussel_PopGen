#!/bin/python
# Name: elai_admix_viz.py
# Purpose: To visualise the global ancestry for each pop on one graph (split violin)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns  # For violin plots

# Function to Load Admixture Data
def load_admix_data(admix_file, population_name):
    """
    Load admixture proportion data from .admix.txt file.
    """
    admix_data = np.loadtxt(admix_file)  # Load as NumPy array
    df = pd.DataFrame(admix_data, columns=["$\it{M. galloprovincialis}$ (Atlantic)", "$\it{M. edulis}$ (Cromarty)"])
    df["Population"] = population_name 
    return df

# Function to Plot Split Violin
def plot_split_violin(admix_files, titles, output_file):
    """
    Creates a split violin plot to show the distribution of admixture proportions.
    """
    all_data = []
    for admix_file, title in zip(admix_files, titles):
        df = load_admix_data(admix_file, title)
        all_data.append(df)

    df_melted = pd.concat(all_data).melt(id_vars=["Population"], var_name="Ancestry", value_name="Proportion")

    fig, ax = plt.subplots(figsize=(10, 6))  # Wider to accommodate legend

    # Muted orange & blue
    colors = ["#e69f00", "#4c72b0"]
    
    sns.violinplot(
        x="Population", y="Proportion", hue="Ancestry", data=df_melted,
        split=True, palette=colors, inner="quartile", saturation=1, ax=ax
    )
    plt.setp(ax.collections, alpha=0.7)

    ax.set_title("Global Admixture Proportions")
    ax.set_ylabel("Ancestry Proportion")
    ax.set_ylim(0, 1)

    # Legend outside plot
    ax.legend(title="Ancestry", bbox_to_anchor=(1.02, 0.5), loc="center left", frameon=False)

    # Adjust layout for legend
    fig.subplots_adjust(right=0.85)
    plt.savefig(f"{output_file}_split_violin.png", dpi=300, bbox_inches="tight")
    plt.show()

# Main Function
def main():
    admix_files = [
        #"ELAI_Results_M.edulis_Nor.admix.txt",
        "ELAI_Results_M.edulis_Eng.admix.txt",
        "ELAI_Results_Ireland.admix.txt",
        "ELAI_Results_Western_Isles.admix.txt",
        "ELAI_Results_Shetland_A.admix.txt",
        "ELAI_Results_Shetland_B.admix.txt",
        "ELAI_Results_Aberdeen.admix.txt"
    ]
    
    #titles = ["Nordic $\it{M. edulis}$", "English $\it{M. edulis}$", "Ireland", "Western Isles", "Shetland", "Aberdeen"]
    titles = ["English $\it{M. edulis}$", "Ireland", "Western Isles", "Shetland_A", "Shetland_B", "Aberdeen"]
    plot_split_violin(admix_files, titles, "global_admixture")

if __name__ == "__main__":
    main()
