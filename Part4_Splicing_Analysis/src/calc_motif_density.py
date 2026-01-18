#!/usr/bin/env python3
import sys
import pandas as pd
from Bio import SeqIO

def count_motifs(seq, motif_def):
    # motif_def: list of sets [{'G'}, {'C'}, {'G','C'}...]
    length = len(motif_def)
    count = 0
    seq = seq.upper()
    
    for i in range(len(seq) - length + 1):
        sub = seq[i : i+length]
        match = True
        for j, char in enumerate(sub):
            if char not in motif_def[j]:
                match = False
                break
        if match:
            count += 1
    return count

def analyze_motif_density(fasta_file, group_label):
    data = []
    # SRSF5 Permissive (WCWWC) -> A/T, C, A/T, A/T, C
    # Or classic Cytosine-rich used in literature (TCCTC, CCTCC, etc.)
    # Let's use the DREME/AME consensus we found: WCWWC 
    # W = A or T (Weak)
    # Actually wait. SRSF5 Consensus from CisBP (Human):
    # C-rich.
    # Let's use the degenerate definition: [ACT] [C] [ACT] [C] [C]
    # No, let's use the classic loose SRSF5:
    # C C A C C (Consensus)
    # Loose: {C,A,T} C {A,T,C} C C
    
    # Actually, let's use the EXACT definition from AME output that was top rank in Lost.
    # M00817 (SRSF5): GCSCC (G, C, G/C, C, C). This was strict.
    
    # Try a broader definition: (C/G) C (Any) C C ?
    # Let's try "C-rich": at least 3 C's in 5-mer?
    # No, be scientific.
    # The AME p-value difference was 10^5.
    # This implies MANY weak sites sum up.
    # Let's try the CORE: C C
    
    # Let's scan for "Y C Y C C" (Py C Py C C) -> SRSF5-like.
    motif_def = [
        {'C','T'},       # Y
        {'C'},           # C
        {'C','T','A','G'}, # N (Spacer) or specific? 
        # Actually standard: T C T C C
        # Let's use:
        {'C','T'}, {'C'}, {'C','T'}, {'C'}, {'C'}
    ] # YCYCC
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        seq_len = len(seq)
        if seq_len == 0: continue
        
        n_hits = count_motifs(seq, motif_def)
        density = (n_hits / seq_len) * 1000 # Hits per 1kb
        
        data.append({
            "Group": group_label,
            "ID": record.id,
            "Length": seq_len,
            "Hits": n_hits,
            "Density_per_kb": density
        })
    return pd.DataFrame(data)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python calc_motif_density.py <lost_fa> <preserved_fa> <out_tsv>")
        sys.exit(1)
        
    df_lost = analyze_motif_density(sys.argv[1], "Lost")
    df_pres = analyze_motif_density(sys.argv[2], "Preserved")
    
    combined = pd.concat([df_lost, df_pres], ignore_index=True)
    combined.to_csv(sys.argv[3], sep="\t", index=False)
    
    # Print Summary Stats
    print("--- Summary Stats (Density per kb) ---")
    print(combined.groupby("Group")["Density_per_kb"].describe())
    
    # Print Percentage with > 0 hits
    print("\n--- Percentage with > 0 Hits ---")
    stats = combined.groupby("Group").apply(lambda x: (x["Hits"] > 0).mean() * 100)
    print(stats)
