#! bin/python
#Name/date: TimRegan/2025.04.22

"""
This code relies on having a reference file which maps out the position of each SNP used in the array. It should have the following format:
Marker ID	Chromosome	Physical Position
AX-603035157	1	60282288
AX-604275399	1	11604362
...etc.
"""

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
    
    # Write the updated file
    updated_map_df.to_csv(output_file, sep="\t", index=False, header=False)
    
    print(f"Updated map file saved to {output_file}")

# Example usage
map_file = "SNP_v3_260325.map"  # Replace with actual path to .map file used
output_file = "SNP_v3_260325_updated.map"
ref_file = "TOL_Med_SNP Pos.txt"  # Replace with actual path to position of reference file. 

update_map_file(map_file, ref_file, output_file)
