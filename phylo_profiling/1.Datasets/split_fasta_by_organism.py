import os
import sys

def split_fasta(input_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    file_handles = {}
    
    try:
        with open(input_file, 'r') as f:
            current_handle = None
            for line in f:
                if line.startswith('>'):
                    # Header format: >organism_id.accession ...
                    # We want the part before the first dot
                    header = line[1:].strip()
                    organism_id = header.split('.')[0]
                    
                    if organism_id not in file_handles:
                        out_path = os.path.join(output_dir, f"{organism_id}.fasta")
                        file_handles[organism_id] = open(out_path, 'w')
                    
                    current_handle = file_handles[organism_id]
                
                if current_handle:
                    current_handle.write(line)
                    
    finally:
        for handle in file_handles.values():
            handle.close()
    
    print(f"Successfully split {input_file} into {len(file_handles)} files in {output_dir}")

if __name__ == "__main__":
    input_fasta = "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/phylo_profiling/1.Datasets/eukaryotes.fasta"
    output_directory = "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/phylo_profiling/1.Datasets/species_fasta"
    
    split_fasta(input_fasta, output_directory)
