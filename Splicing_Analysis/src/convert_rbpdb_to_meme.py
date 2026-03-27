#!/usr/bin/env python3
import os
import glob
import numpy as np
import argparse

def parse_pfm_file(fpath):
    """
    Reads a PFM file. Expects 4 rows (A, C, G, T) with space-separated counts.
    Returns ID and a list of probability rows (for MEME output).
    """
    matrix = []
    try:
        with open(fpath) as f:
            for line in f:
                parts = line.strip().split()
                if not parts: continue
                # Parse as floats/ints
                matrix.append([float(x) for x in parts])
    except Exception as e:
        print(f"[WARN] Failed to read {fpath}: {e}")
        return None, None
        
    if len(matrix) != 4:
        print(f"[WARN] {fpath} does not have 4 rows (A, C, G, T). Skipping.")
        return None, None
        
    # Convert to Numpy for easier transposition and normalization
    mat = np.array(matrix) # Shape (4, Width)
    
    # Check for empty
    if mat.size == 0: return None, None
    
    # Transpose to (Width, 4) -> Rows are positions, Cols are A C G T
    mat_t = mat.T
    
    # Normalize to probabilities
    # Add pseudocount? MEME usually handles raw counts or probs. 
    # letter-probability matrix says "alength= 4 w= ... nsites= ..."
    # P(x) = count(x) / sum(counts)
    
    probs = []
    nsites = 0 # Max sum across columns?
    
    for row in mat_t:
        total = sum(row)
        if total == 0:
            # Uniform if no counts? Or skip?
            prob_row = [0.25, 0.25, 0.25, 0.25]
        else:
            prob_row = [x / total for x in row]
            nsites = max(nsites, total)
        probs.append(prob_row)
        
    fname = os.path.basename(fpath)
    motif_id = os.path.splitext(fname)[0]
    
    return motif_id, probs, nsites

def write_meme(pfm_dir, out_file):
    print(f"Converting PFMs from {pfm_dir} to {out_file}...")
    files = glob.glob(os.path.join(pfm_dir, "*.pfm"))
    files.sort()
    
    if not files:
        print(f"[WARN] No .pfm files found in {pfm_dir}")
        return

    with open(out_file, "w") as out:
        # Header
        out.write("MEME version 4\n\n")
        out.write("ALPHABET= ACGT\n\n")
        out.write("strands: + -\n\n")
        out.write("Background letter frequencies\n")
        out.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")
        
        count = 0
        for fpath in files:
            mid, probs, nsites = parse_pfm_file(fpath)
            if mid is None: continue
            
            width = len(probs)
            out.write(f"MOTIF {mid} {mid}\n")
            out.write(f"letter-probability matrix: alength= 4 w= {width} nsites= {int(nsites)} E= 0\n")
            for row in probs:
                out.write(f"{row[0]:.6f} {row[1]:.6f} {row[2]:.6f} {row[3]:.6f}\n")
            out.write("\n")
            count += 1
            
        print(f"Wrote {count} motifs to {out_file}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--human_pfm_dir", help="Path to Human PFM directory")
    parser.add_argument("--mouse_pfm_dir", help="Path to Mouse PFM directory")
    parser.add_argument("--out_dir", required=True)
    args = parser.parse_args()
    
    os.makedirs(args.out_dir, exist_ok=True)
    
    if args.human_pfm_dir:
        write_meme(args.human_pfm_dir, os.path.join(args.out_dir, "RBPDB_Human.meme"))
        
    if args.mouse_pfm_dir:
        write_meme(args.mouse_pfm_dir, os.path.join(args.out_dir, "RBPDB_Mouse.meme"))

if __name__ == "__main__":
    main()
