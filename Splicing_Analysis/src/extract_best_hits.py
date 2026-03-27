#!/usr/bin/env python3
import sys
import os
import random
from Bio import SeqIO

def get_hamming_dist(seq, motif_consensus):
    # motif_consensus: list of sets of allowed chars
    # seq: string
    dist = 0
    for i, char in enumerate(seq):
        allowed = motif_consensus[i]
        if char not in allowed:
            dist += 1
    return dist

def extract_best_hit(fasta_in, fasta_out, motif_def):
    # motif_def: list of sets, e.g. [{'G'}, {'C'}, {'G','C'}, {'C'}, {'C'}]
    motif_len = len(motif_def)
    
    with open(fasta_out, 'w') as out:
        for record in SeqIO.parse(fasta_in, "fasta"):
            seq = str(record.seq).upper()
            if len(seq) < motif_len: continue
            
            best_seq = ""
            min_dist = motif_len + 1
            
            # Simple sliding window
            for i in range(len(seq) - motif_len + 1):
                sub = seq[i : i+motif_len]
                dist = get_hamming_dist(sub, motif_def)
                if dist < min_dist:
                    min_dist = dist
                    best_seq = sub
                if dist == 0: break # Found perfect match
            
            if best_seq and min_dist <= 1:
                # Good match (0-1 mismatch)
                out.write(f">{record.id}_dist{min_dist}\n{best_seq}\n")
            else:
                # No Good Match -> Inject Random Noise
                # This ensures the Logo 'Bits' drop (entropy increases) without creating an 'N' block.
                # Generate a random 5-mer
                noise = "".join(random.choices(['A','C','G','T'], k=5))
                out.write(f">{record.id}_no_hit\n{noise}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_best_hits.py <input.fa> <output.fa>")
        sys.exit(1)
        
    # SRSF5 Consensus: G C S C C 
    # S = G or C
    motif_gcscc = [{'G'}, {'C'}, {'G','C'}, {'C'}, {'C'}]
    
    extract_best_hit(sys.argv[1], sys.argv[2], motif_gcscc)
