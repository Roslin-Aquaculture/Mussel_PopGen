#!/bin/python
#Name/date: TimRegan/2025.08.05

import pandas as pd

def update_map_file(map_file, ref_file, output_file):
    # Load the PLINK map file
    map_cols = ["Chromosome", "Marker ID", "Genetic Distance", "Physical Position"]
    map_df = pd.read_csv(map_file, sep="\t", comment='#', names=map_cols, dtype={"Chromosome": str})
    
    # Load the reference file
    ref_df = pd.read_csv(ref_file, sep="\t", dtype={"Chromosome": str})
    
    # Merge to update chromosome and physical position
    updated_map_df = map_df.merge(ref_df, on="Marker ID", how="left", suffixes=("", "_new"))
    
    # Replace old chromosome and position with new values where available
    updated_map_df["Chromosome"] = updated_map_df["Chromosome_new"].fillna(updated_map_df["Chromosome"])
    updated_map_df["Physical Position"] = updated_map_df["Physical Position_new"].fillna(updated_map_df["Physical Position"])
    
    # Drop extra columns
    updated_map_df = updated_map_df[map_cols]

    # Replace any "---" with "0" in Physical position column
    updated_map_df["Physical Position"] = updated_map_df["Physical Position"].replace("---", "0")
    
    # Write the updated file
    updated_map_df.to_csv(output_file, sep="\t", index=False, header=False)
    
    print(f"Updated map file saved to {output_file}")

# Example usage
map_file = "SNP_v3_110825.map"  # Replace with actual path
ref_file = "TOL_Med_SNP Pos.txt"  # Replace with actual path
output_file = "SNP_v3_110825.map"
update_map_file(map_file, ref_file, output_file)
