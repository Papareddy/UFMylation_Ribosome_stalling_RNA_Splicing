import sys
import os

def convert_cisbp_to_meme(input_path, output_file):
    """
    Converts CisBP-RNA motifs to MEME format.
    Handles both:
    1. A single 'PWM.txt' file (bulk format type A).
    2. A directory containing 'RBP_Information_all_motifs.txt' and a 'pwms_all_motifs' subdirectory (bulk format type B).
    """
    
    motifs = {} # ID -> {name: "", matrix: []}
    
    # CASE 1: Input is a Directory (New Bulk Format, e.g., Mouse)
    if os.path.isdir(input_path):
        print(f"Detected directory input: {input_path}")
        
        # 1. Load RBP Info Map
        info_file = os.path.join(input_path, "RBP_Information_all_motifs.txt")
        id_to_name = {}
        if os.path.exists(info_file):
            with open(info_file, 'r') as f:
                header = f.readline().strip().split('\t')
                # Find indices
                try:
                    idx_motif_id = header.index("Motif_ID")
                    idx_rbp_name = header.index("RBP_Name")
                except ValueError:
                    print("Error: Could not find Motif_ID or RBP_Name in RBP_Information_all_motifs.txt header.")
                    sys.exit(1)
                    
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) > max(idx_motif_id, idx_rbp_name):
                        m_id = parts[idx_motif_id]
                        r_name = parts[idx_rbp_name]
                        if m_id != "." and m_id != "":
                            id_to_name[m_id] = r_name
        else:
            print(f"Warning: {info_file} not found. Using IDs as names.")

        # 2. Iterate over PWM files
        pwm_dir = os.path.join(input_path, "pwms_all_motifs")
        if not os.path.exists(pwm_dir):
             print(f"Error: {pwm_dir} not found.")
             sys.exit(1)
             
        files = [f for f in os.listdir(pwm_dir) if f.endswith(".txt")]
        print(f"Found {len(files)} PWM files.")
        
        for fname in files:
            motif_id = fname.replace(".txt", "")
            rbp_name = id_to_name.get(motif_id, motif_id)
            
            matrix = []
            with open(os.path.join(pwm_dir, fname), 'r') as f:
                # Header: Pos A C G U
                header = f.readline() 
                for line in f:
                    if not line.strip(): continue
                    parts = line.strip().split('\t')
                    # Format: Pos ProbA ProbC ProbG ProbU
                    # MEME expects: ProbA ProbC ProbG ProbT
                    # Parts[1]=A, Parts[2]=C, Parts[3]=G, Parts[4]=U
                    row = f"{parts[1]}\t{parts[2]}\t{parts[3]}\t{parts[4]}"
                    matrix.append(row)
            
            if matrix:
                motifs[motif_id] = {"name": rbp_name, "matrix": matrix}

    # CASE 2: Input is a single File (Old Bulk Format)
    else:
        print(f"Detected file input: {input_path}")
        
        motifs_list = [] # Temporarily use a list to match old logic
        current_motif = {}
        matrix_lines = []
        reading_matrix = False
        
        with open(input_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    if current_motif and matrix_lines:
                        current_motif['matrix'] = matrix_lines
                        motifs_list.append(current_motif)
                        current_motif = {}
                        matrix_lines = []
                        reading_matrix = False
                    continue
                
                parts = line.split('\t')
                key = parts[0]
                
                if key == "RBP":
                    if current_motif: 
                         pass 
                    current_motif = {}
                elif key == "RBP Name":
                    current_motif['alt_name'] = parts[1] if len(parts) > 1 else "Unknown"
                elif key == "Motif":
                    current_motif['id'] = parts[1] if len(parts) > 1 else "Unknown"
                elif key == "Pos":
                    reading_matrix = True
                    continue
                elif reading_matrix:
                    if len(parts) >= 5:
                        probs = parts[1:5]
                        matrix_lines.append(" ".join(probs))
                        
            if current_motif and matrix_lines:
                current_motif['matrix'] = matrix_lines
                motifs_list.append(current_motif)
        
        # Convert list to dictionary for common output
        for m in motifs_list:
            motif_id = m.get('id', 'Unknown')
            alt_name = m.get('alt_name', 'Unknown')
            motifs[motif_id] = {"name": alt_name, "matrix": m['matrix']}


    # Write Output (Common)
    with open(output_file, 'w') as out:
        out.write("MEME version 4\n\n")
        out.write("ALPHABET= ACGU\n\n")
        out.write("strands: + -\n\n")
        out.write("Background letter frequencies\n")
        out.write("A 0.25 C 0.25 G 0.25 U 0.25\n\n")
        
        for m_id, data in motifs.items():
            out.write(f"MOTIF {m_id} {data['name']}\n")
            out.write(f"letter-probability matrix: alength= 4 w= {len(data['matrix'])} nsites= 20 E= 0\n")
            for row in data['matrix']:
                out.write(f"{row}\n")
            out.write("\n")
            
    print(f"Converted {len(motifs)} motifs.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python cisbp_to_meme.py <input_path> <output_meme>")
        print("  <input_path>: Can be 'PWM.txt' file OR the top-level CisBP directory.")
        sys.exit(1)
        
    convert_cisbp_to_meme(sys.argv[1], sys.argv[2])
