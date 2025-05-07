#! bin/python
#Name: elai_admix_viz.py
#Purpose: To visualise the global ancestry for each pop on one graph
 
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
    df = pd.DataFrame(admix_data, columns=["M. gallo_Atl", "Cromarty (M. edulis)"])
     
    # Add population column for hue
    df["Population"] = population_name 
    return df
 
# Function to Plot Split Violin
def plot_split_violin(admix_files, titles, output_file):
    """
    Creates a split violin plot to show the distribution of admixture proportions.
    """
    all_data = []  # Store all population data
 
    # Load each population's data
    for admix_file, title in zip(admix_files, titles):
        df = load_admix_data(admix_file, title)
        all_data.append(df)
 
    # Combine into one DataFrame
    df_melted = pd.concat(all_data).melt(id_vars=["Population"], var_name="Ancestry", value_name="Proportion")
 
    plt.figure(figsize=(8, 6))
 
    # Define Muted Colors (Matching Per-Chromosome Graphs)
    colors = ["#e69f00", "#4c72b0"]  # Muted orange & blue 
 
    # Create split violin plot
    ax = sns.violinplot(x="Population", y="Proportion", hue="Ancestry", data=df_melted,
                        split=True, palette=colors, inner="quartile", saturation=1)
    plt.setp(ax.collections, alpha=0.7)
 
    plt.title("Global Admixture Proportions")
    plt.ylabel("Ancestry Proportion")
    plt.ylim(0, 1)  # Full scale (0-100%)
 
    plt.tight_layout()
    plt.savefig(f"{output_file}_split_violin.png")
    plt.show()
 
# Main Function to Run for All Populations
def main():
    admix_files = [
        "ELAI_Results_WesIsl.admix.txt",
        "ELAI_Results_Shetland.admix.txt",
        "ELAI_Results_Irl_Atl.admix.txt",
        "ELAI_Results_M.edulis_Eng.admix.txt",
        "ELAI_Results_M.edulis_Nor.admix.txt",
        "ELAI_Results_Aberdeen.admix.txt"
    ]
     
    titles = ["Western Isles", "Shetland", "Ireland", "M. edulis Eng", "M. edulis Nor", "Aberdeen"]
 
    plot_split_violin(admix_files, titles, "global_admixture")
 
if __name__ == "__main__":
    main()
