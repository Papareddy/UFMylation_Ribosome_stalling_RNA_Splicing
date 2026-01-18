
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main():
    parser = argparse.ArgumentParser(description="Plot protein impact fractions.")
    parser.add_argument("--per_event_file", required=True, help="Path to per_event_compact_for_plotting.tsv")
    parser.add_argument("--output_file", required=True, help="Path to save the output plot (e.g., plot.png)")
    args = parser.parse_args()

    # Read data
    perE = pd.read_csv(args.per_event_file, sep="\t")

    impact_order = [
        "Frameshift_NMD_likely",
        "Frameshift_no_NMD",
        "Start_Stop_disruption",
        "Inframe_CDS_change",
        "CDS_neutral",
        "UTR_only"
    ]

    perE = perE[perE['protein_impact_class'].isin(impact_order)]
    perE['protein_impact_class'] = pd.Categorical(perE['protein_impact_class'], categories=impact_order, ordered=True)

    # Fractions
    tab = pd.crosstab(perE['dataset'], perE['protein_impact_class'])
    frac = tab.div(tab.sum(axis=1), axis=0)

    # Colors per class
    cols = {
        "Frameshift_NMD_likely": "#b2182b",
        "Frameshift_no_NMD": "#ef8a62",
        "Start_Stop_disruption": "#fddbc7",
        "Inframe_CDS_change": "#d1e5f0",
        "CDS_neutral": "#67a9cf",
        "UTR_only": "#2166ac"
    }
    plot_cols = [cols[c] for c in frac.columns]

    # Check for Directional Datasets
    datasets = set(frac.index)
    datasets = set(frac.index)
    is_directional = all(x in datasets for x in ["UFM1_dependent_dPSI_positive", "UFM1_dependent_dPSI_negative", "UFM1_independent_dPSI_positive", "UFM1_independent_dPSI_negative"])
    
    if is_directional:
        # Create 2 subplots: Inclusion and Exclusion
        fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
        
        # Subplot 1: Inclusion -> dPSI_positive
        ax = axes[0]
        p_inc = frac.loc["UFM1_independent_dPSI_positive"]
        l_inc = frac.loc["UFM1_dependent_dPSI_positive"]
        
        # Order by Preserved Inc
        idx_inc = p_inc.sort_values(ascending=False).index
        p_inc = p_inc[idx_inc]
        l_inc = l_inc[idx_inc]
        pc_inc = [cols[c] for c in idx_inc]
        
        ax.bar(l_inc.index, l_inc, color=pc_inc, edgecolor=None)
        ax.bar(p_inc.index, -p_inc, color=[adjust_alpha(c, 0.5) for c in pc_inc], edgecolor=None)
        ax.set_title("Inclusion Events (UFM1_Dep vs UFM1_Indep)")
        ax.axhline(0, color='black', linewidth=1)
        ax.set_xticklabels(l_inc.index, rotation=45, ha='right')
        ax.set_ylabel("Fraction of events")

        # Subplot 2: Exclusion -> dPSI_negative
        ax = axes[1]
        p_exc = frac.loc["UFM1_independent_dPSI_negative"]
        l_exc = frac.loc["UFM1_dependent_dPSI_negative"]
        
        # Order by Preserved Exc
        idx_exc = p_exc.sort_values(ascending=False).index
        p_exc = p_exc[idx_exc]
        l_exc = l_exc[idx_exc]
        pc_exc = [cols[c] for c in idx_exc]
        
        ax.bar(l_exc.index, l_exc, color=pc_exc, edgecolor=None)
        ax.bar(p_exc.index, -p_exc, color=[adjust_alpha(c, 0.5) for c in pc_exc], edgecolor=None)
        ax.set_title("Exclusion Events (UFM1_Dep vs UFM1_Indep)")
        ax.axhline(0, color='black', linewidth=1)
        ax.set_xticklabels(l_exc.index, rotation=45, ha='right')

        plt.suptitle("Protein Impact Class Fractions (Directional)", fontsize=14)
        plt.tight_layout()
        plt.savefig(args.output_file, dpi=300)
        print(f"Plot saved to {args.output_file} (Directional)")
        
    else:
        # Standard Plot (UFM1_Dep vs UFM1_Indep)
        if "UFM1_independent" in frac.index and "UFM1_dependent" in frac.index:
            # Order by preserved
            frac = frac.loc[:, frac.loc['UFM1_independent'].sort_values(ascending=False).index]
            preserved = frac.loc["UFM1_independent"]
            lost = frac.loc["UFM1_dependent"]
            plot_cols = [cols[c] for c in frac.columns]
            
            plt.figure(figsize=(10, 6))
            plt.bar(lost.index, lost, color=plot_cols, edgecolor=None)
            plt.bar(preserved.index, -preserved, color=[adjust_alpha(c, 0.5) for c in plot_cols], edgecolor=None)

            plt.xticks(rotation=45, ha='right')
            plt.ylabel("Fraction of events")
            plt.title("Protein Impact Class Fractions (UFM1_Dep vs. UFM1_Indep)")
            plt.axhline(0, color='black', linewidth=1)
            plt.tight_layout()
            plt.savefig(args.output_file, dpi=300)
            print(f"Plot saved to {args.output_file}")
        else:
            print(f"[WARN] Expected sets (lost, preserved) or directional sets not found in {datasets}. Skipping plot.")

def adjust_alpha(hex_color, alpha):
    """Adjust alpha of a hex color string."""
    rgb = tuple(int(hex_color.lstrip('#')[i:i+2], 16) / 255.0 for i in (0, 2, 4))
    return (*rgb, alpha)

if __name__ == "__main__":
    main()
