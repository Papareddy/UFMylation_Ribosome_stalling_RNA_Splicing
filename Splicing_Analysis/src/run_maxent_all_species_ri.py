#!/usr/bin/env python3
import os
import subprocess
import shutil

def run_cmd(cmd):
    full_cmd = ["mamba", "run", "-n", "splicing-functional"] + cmd
    print(f"[EXEC] {' '.join(full_cmd)}")
    subprocess.check_call(full_cmd)

def main():
    project_root = "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis"
    script_dir = os.path.join(project_root, "src")
    results_dir = os.path.join(project_root, "results")
    consolidated_dir = os.path.join(results_dir, "MaxEnt_RI_Consolidated")
    os.makedirs(consolidated_dir, exist_ok=True)

    species_config = {
        "human": {
            "lost": os.path.join(results_dir, "human/nucleus/step01_data_prep/UFM1_dependent.tsv"),
            "preserved": os.path.join(results_dir, "human/nucleus/step01_data_prep/UFM1_independent.tsv"),
            "genome": os.path.join(project_root, "data/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
            "gtf": os.path.join(project_root, "data/human/pcg_gencode.v45.annotation.gtf.gz")
        },
        "mouse": {
            "lost": os.path.join(results_dir, "mouse/total/step01_data_prep/UFM1_dependent.tsv"),
            "preserved": os.path.join(results_dir, "mouse/total/step01_data_prep/UFM1_independent.tsv"),
            "genome": os.path.join(project_root, "data/mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa"),
            "gtf": os.path.join(project_root, "data/mouse/Mus_musculus.GRCm39.112.gtf.gz")
        },
        "arabidopsis": {
            "lost": os.path.join(results_dir, "arabidopsis/nucleus/step01_data_prep/UFM1_dependent.tsv"),
            "preserved": os.path.join(results_dir, "arabidopsis/nucleus/step01_data_prep/UFM1_independent.tsv"),
            "genome": os.path.join(project_root, "data/arabidopsis/Arabidopsis_thaliana_TAIR10.dna.primary_assembly.fa"),
            "gtf": os.path.join(project_root, "data/arabidopsis/no_plastid_no_rRNA.Arabidopsis_thaliana.TAIR10.56.gtf")
        }
    }

    for species, paths in species_config.items():
        print(f"\n=== Processing {species} ===")
        species_outdir = os.path.join(consolidated_dir, species)
        os.makedirs(species_outdir, exist_ok=True)

        # 1. Run RI events (Lost/Preserved) using extract_ri_motifs.py
        cmd_ri = [
            "python3", os.path.join(script_dir, "extract_ri_motifs.py"),
            "--lost", paths["lost"],
            "--preserved", paths["preserved"],
            "--genome_fasta", paths["genome"],
            "--outdir", species_outdir
        ]
        try:
            run_cmd(cmd_ri)
        except Exception as e:
            print(f"[ERROR] RI Extraction failed for {species}: {e}")

        # 2. Process Constitutive Introns
        print(f"[INFO] Generating and scoring Constitutive Introns for {species}...")
        
        # Generate constitutive BEDs
        const_bed_base = os.path.join(species_outdir, "constitutive")
        cmd_const_gen = [
            "python3", os.path.join(script_dir, "generate_constitutive_introns.py"),
            "--gtf", paths["gtf"],
            "--exclude_files", paths["lost"], paths["preserved"],
            "--out_bed", f"{const_bed_base}.bed",
            "--n_sample", "5000"
        ]
        try:
            run_cmd(cmd_const_gen)
            
            # Extract and Score 5' and 3' splice sites for constitutive
            for kind in ["5ss", "3ss"]:
                bed = f"{const_bed_base}.{kind}.bed"
                fa = f"{const_bed_base}.{kind}.fa"
                
                # Extract FASTA
                cmd_getfasta = [
                    "bedtools", "getfasta", "-s", "-fi", paths["genome"], 
                    "-bed", bed, "-fo", fa, "-name"
                ]
                run_cmd(cmd_getfasta)
                
                # Score using the same logic (we'll call a helper if we had one, but we'll use a hack)
                # Actually, extract_ri_motifs.py has a score_ss function.
                # Let's add a small helper script or just run a python snippet.
                score_cmd = [
                    "python3", "-c", 
                    f"import sys; sys.path.append('{script_dir}'); import extract_ri_motifs; extract_ri_motifs.score_ss('{fa}', '{kind}')"
                ]
                run_cmd(score_cmd)
                
                # Copy to consolidated folder
                dst = os.path.join(consolidated_dir, f"{species}_constitutive_{kind}.scores.tsv")
                src = fa + ".scores"
                if os.path.exists(src):
                    shutil.copy(src, dst)
                    print(f"[INFO] Copied {src} to {dst}")
                    
        except Exception as e:
            print(f"[ERROR] Constitutive processing failed for {species}: {e}")

        # 3. Consolidate results for Lost/Preserved (from step 1)
        for grp in ["UFM1_dependent", "UFM1_independent"]:
            for kind in ["5ss", "3ss"]:
                src = os.path.join(species_outdir, f"{grp}.{kind}.fa.scores")
                if os.path.exists(src):
                    dst = os.path.join(consolidated_dir, f"{species}_{grp}_{kind}.scores.tsv")
                    shutil.copy(src, dst)
                    print(f"[INFO] Copied {src} to {dst}")

    print(f"\n[DONE] Consolidated results in {consolidated_dir}")

if __name__ == "__main__":
    main()
