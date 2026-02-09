import os
import pandas as pd

def aggregate_ri():
    """
    Aggregates RI events from Human (nucleus), Mouse (total), and Arabidopsis (nucleus)
    into a consolidated data frame.
    """
    # Base directory for Part4_Splicing_Analysis
    base_dir = "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis"
    
    # Define species and their respective fraction subdirectories
    species_fractions = {
        "human": "nucleus",
        "mouse": "total",
        "arabidopsis": "nucleus"
    }

    # Consolidated data list
    consolidated_data = []

    for species, fraction in species_fractions.items():
        print(f"Processing {species} ({fraction})...")
        
        # Results directory for Step 01 output
        step1_out = os.path.join(base_dir, "results", species, fraction, "step01_data_prep")
        
        # Files to process
        target_files = {
            "Dependent": "UFM1_dependent.tsv",
            "Independent": "UFM1_independent.tsv"
        }

        for dependency, filename in target_files.items():
            filepath = os.path.join(step1_out, filename)
            
            if not os.path.exists(filepath):
                print(f"  [WARNING] File not found: {filepath}")
                continue

            # Load data
            df = pd.read_csv(filepath, sep='\t')

            # Filter for RI events
            ri_df = df[df['EventType'] == 'RI'].copy()

            if ri_df.empty:
                print(f"  [INFO] No RI events in {filename}")
                continue

            # Column extraction and sign-flip
            # User requirement: PSI positive = Retention in Control (WT)
            # rMATS IncLevelDifference = PSI_Anisomycin - PSI_WT
            # Therefore, we use -IncLevelDifference
            
            # Select relevant columns
            # Column mapping:
            # Species -> species
            # geneSymbol -> gene name (fallback to GeneID if symbol is missing/NA)
            # -dPSI_num.WT -> PSI (Sign flipped)
            
            # Note: dPSI_num.WT is used in prepare_rmats_data.R for filtering
            # and corresponds to IncLevelDifference.WT (which is case - ctrl)
            
            # Process gene name (Symbol fallback to ID)
            ri_df['GeneName'] = ri_df['geneSymbol'].fillna(ri_df['GeneID'])
            
            # Calculate final PSI (Percentage, Sign flipped)
            # Using dPSI_num.WT which is the filtered IncLevelDifference value
            ri_df['PSI_Percentage'] = -ri_df['dPSI_num.WT'] * 100

            # Extract final columns
            extracted = ri_df[['GeneID', 'GeneName', 'PSI_Percentage']].copy()

            extracted['Species'] = species
            extracted['Dependency'] = dependency
            
            consolidated_data.append(extracted)
            print(f"  [SUCCESS] Extracted {len(extracted)} RI events from {filename}")

    if not consolidated_data:
        print("[ERROR] No data aggregated.")
        return

    # Combine all dataframes
    final_df = pd.concat(consolidated_data, ignore_index=True)

    # Reorder columns
    final_df = final_df[['Species', 'GeneID', 'GeneName', 'PSI_Percentage', 'Dependency']]


    # Output path
    output_path = os.path.join(base_dir, "results", "RI_Consolidated_Aggregation.tsv")
    final_df.to_csv(output_path, sep='\t', index=False)
    
    print(f"\n[COMPLETE] Consolidated data saved to: {output_path}")
    print(f"Total RI events: {len(final_df)}")
    print(final_df.head())

if __name__ == "__main__":
    aggregate_ri()
