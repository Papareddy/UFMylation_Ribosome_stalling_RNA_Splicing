
import os
import subprocess
import glob

# Configuration
ROOT_DIR = "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis"
MOTIF_DB_MOUSE = os.path.join(ROOT_DIR, "data/motifs/CisBP_Mouse_All.meme")
MOTIF_DB_HUMAN = os.path.join(ROOT_DIR, "data/motifs/CisBP_Human_All.meme")

# List of corrupted directories (identified previously)
# I will use a glob to find them again just to be robust, checking for missing header in ame.tsv
def find_corrupted_ame_dirs():
    corrupted_dirs = []
    # Search recursively for ame.tsv
    # Limit depth to avoid scanning everything if possible, but full scan is fine
    for ame_file in glob.glob(os.path.join(ROOT_DIR, "**", "ame.tsv"), recursive=True):
        if "results" in ame_file and "coverage_test" not in ame_file:
             # Skip old results if strictly focusing on coverage_test, but user wanted to fix "results" too if relevant.
             # Actually, let's fix everything we find corrupted.
             pass
        
        try:
            with open(ame_file, 'r') as f:
                first_line = f.readline()
                if not first_line.startswith("rank"):
                    corrupted_dirs.append(os.path.dirname(ame_file))
        except Exception as e:
            print(f"[WARN] Could not read {ame_file}: {e}")
    return corrupted_dirs

def fix_ame_run(ame_dir):
    print(f"\n[FIX] Fixing: {ame_dir}")
    parent_dir = os.path.dirname(ame_dir)
    dir_name = os.path.basename(ame_dir)
    
    # Parse directory name: ame_{REGION}_UFM1_{TYPE}_vs_constitutive
    # e.g. ame_3ss_UFM1_dependent_vs_constitutive
    try:
        parts = dir_name.split('_')
        # name parts: ame, region, UFM1, type, vs, constitutive
        # region can be "3ss", "5ss", "intron", "exon"
        # However, for exon, it might be ame_exon_...
        
        region = parts[1] # 3ss, 5ss, intron, exon
        
        # Determine Type: dependent or independent
        # If "dependent", it corresponds to "lost" (RI) or specific logic for SE
        analysis_type_str = ""
        if "dependent" in dir_name and "independent" not in dir_name:
            analysis_type_str = "dependent"
        elif "independent" in dir_name:
            analysis_type_str = "independent"
        else:
            print(f"[ERROR] Could not determine type (dep/indep) from {dir_name}")
            return

        # Determine Species
        if "/human/" in ame_dir:
            species = "human"
            motif_db = MOTIF_DB_HUMAN
        else:
            species = "mouse"
            motif_db = MOTIF_DB_MOUSE
            
        # Determine Input Files based on Directory Type (RI vs SE)
        # RI paths: .../RI_DeepDive_.../
        # SE paths: .../SE_DeepDive_.../
        
        is_ri = "RI_DeepDive" in parent_dir
        is_se = "SE_DeepDive" in parent_dir
        
        target_fa = ""
        control_fa = ""
        
        if is_ri:
            # RI Logic:
            # Dependent -> lost.{region}.fa
            # Independent -> preserved.{region}.fa
            # Control -> constitutive_introns.{region}.fa
            
            if analysis_type_str == "dependent":
                target_fa = os.path.join(parent_dir, f"lost.{region}.fa")
            else:
                target_fa = os.path.join(parent_dir, f"preserved.{region}.fa")
                
            control_fa = os.path.join(parent_dir, f"constitutive_introns.{region}.fa")
            
        elif is_se:
            # SE Logic:
            # Dependent -> included.{region}.fa (or skipped?) 
            # Wait, SE analysis usually compares "Included" vs "Excluded" or "Regulated" vs "Background"
            # In analyze_se_vs_constitutive.py:
            # "UFM1_dependent" maps to args.se_fa (Sensitive)
            # "UFM1_independent" maps to args.control_fa (Control/Insensitive)
            # BUT both are compared against "constitutive_exons".
            # The filenames in parent dir:
            # sensitive.exon.fa (Dependent)
            # insensitive.exon.fa (Independent)
            # constitutive_exons.exon.fa
            
            # Let's check files in parent if possible, or assume logic.
            # Assuming 'sensitive' = dependent, 'insensitive' = independent
            
            # region for SE is usually 'exon'
            
            if analysis_type_str == "dependent":
                target_fa = os.path.join(parent_dir, f"sensitive.{region}.fa")
            else:
                target_fa = os.path.join(parent_dir, f"insensitive.{region}.fa")
            
            control_fa = os.path.join(parent_dir, f"constitutive_exons.{region}.fa")
            
        else:
            print("[WARN] Unknown analysis type (not RI/SE DeepDive). Skipping.")
            return

        # Verification
        if not os.path.exists(target_fa):
            print(f"[ERROR] Target FASTA not found: {target_fa}")
            # Try fallback for RI: maybe "intron" instead of "3ss" if 3ss missing?
            # But the dir name is ame_3ss...
            return
        if not os.path.exists(control_fa):
            print(f"[ERROR] Control FASTA not found: {control_fa}")
            return
            
        # Run AME
        # cmd: mamba run -n meme_env ame --evalue-report-threshold 1000 --control {control} --oc {out} {target} {db}
        # Note: Added --verbose 1 to reduce output
        cmd = [
            "mamba", "run", "-n", "meme_env",
            "ame", "--verbose", "1",
            "--evalue-report-threshold", "1000",
            "--control", control_fa,
            "--oc", ame_dir,
            target_fa,
            motif_db
        ]
        
        print(f"[EXEC] {' '.join(cmd)}")
        subprocess.check_call(cmd)
        print("[SUCCESS] AME re-run completed.")
        
    except Exception as e:
        print(f"[ERROR] Failed to fix {ame_dir}: {e}")

def main():
    print("Scanning for corrupted ame.tsv files...")
    corrupted = find_corrupted_ame_dirs()
    print(f"Found {len(corrupted)} corrupted directories.")
    
    for d in corrupted:
        fix_ame_run(d)

if __name__ == "__main__":
    main()
