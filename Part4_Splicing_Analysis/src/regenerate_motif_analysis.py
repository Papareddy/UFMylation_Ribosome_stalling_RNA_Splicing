
import os
import subprocess
import glob
import sys

# Configuration
ROOT_DIR = "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis"
SCRIPT_DIR = os.path.join(ROOT_DIR, "src")

# Paths
DATA_DIR = os.path.join(ROOT_DIR, "data")
GENOMES_DIR = os.path.join(DATA_DIR, "genomes")

# Mappings
SPECIES_MAP = {
    "mouse": {
        "fasta": os.path.join(DATA_DIR, "mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa"),
        "gtf": os.path.join(DATA_DIR, "mouse/Mus_musculus.GRCm39.112.gtf.gz"),
        "motif_db": os.path.join(DATA_DIR, "motifs/CisBP_Mouse_All.meme")
    },
    "human": {
        "fasta": os.path.join(DATA_DIR, "human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
        "gtf": os.path.join(DATA_DIR, "human/pcg_gencode.v45.annotation.gtf.gz"),
        "motif_db": os.path.join(DATA_DIR, "motifs/CisBP_Human_All.meme")
    }
}
# Note: Check if Human/Human.fa exists or if it is distinct. 
# Usually I download genomes to Mouse/Mouse.fa etc.

def run_cmd(cmd, desc):
    print(f"\n[EXEC] {desc}")
    print(" ".join(cmd))
    subprocess.check_call(cmd)

def regenerate_analysis(run_outdir):
    print(f"Scanning {run_outdir}...")
    
    # Identify Species
    species = "mouse" if "/mouse/" in run_outdir else "human"
    conf = SPECIES_MAP[species]
    
    # Check Step 1 Exists
    step1_out = os.path.join(run_outdir, "step01_data_prep")
    if not os.path.exists(step1_out):
        print(f"[WARN] Step 1 missing in {run_outdir}. Skipping.")
        return

    step10_out = os.path.join(run_outdir, "step10_motif_analysis")
    if not os.path.exists(step10_out):
        print(f"[WARN] Step 10 missing in {run_outdir}. Skipping.")
        return

    # Find DeepDive folders
    # RI_DeepDive_comparison
    # SE_DeepDive_comparison
    
    for dd_dir in glob.glob(os.path.join(step10_out, "*_DeepDive_*")):
        dirname = os.path.basename(dd_dir)
        print(f"Found Analysis Dir: {dirname}")
        
        # Parse: {TYPE}_DeepDive_{COMPARISON}
        parts = dirname.split("_DeepDive_")
        if len(parts) != 2: continue
        
        etype = parts[0] # RI or SE
        comparison = parts[1] # combined, dPSI_positive, etc.
        
        # Identify Inputs
        if comparison == "combined":
             l_name, p_name = "UFM1_dependent.tsv", "UFM1_independent.tsv"
        elif comparison == "dPSI_positive":
             l_name, p_name = "UFM1_dependent_dPSI_positive.tsv", "UFM1_independent_dPSI_positive.tsv"
        elif comparison == "dPSI_negative":
             l_name, p_name = "UFM1_dependent_dPSI_negative.tsv", "UFM1_independent_dPSI_negative.tsv"
        else:
             print(f"Unknown comparison: {comparison}. Skipping.")
             continue
             
        l_path = os.path.join(step1_out, l_name)
        p_path = os.path.join(step1_out, p_name)
        
        if not os.path.exists(l_path):
            print(f"Input missing: {l_path}. Skipping.")
            continue
            
        # 1. Extraction
        if etype == "RI":
            extract_script = os.path.join(SCRIPT_DIR, "extract_ri_motifs.py")
            analyze_script = os.path.join(SCRIPT_DIR, "analyze_ri_vs_constitutive.py")
        elif etype == "SE":
            extract_script = os.path.join(SCRIPT_DIR, "extract_se_motifs.py")
            # Assuming analyze_se_vs_constitutive.py exists and matches args
            analyze_script = os.path.join(SCRIPT_DIR, "analyze_se_vs_constitutive.py")
        else:
            continue
            
        # Run Extraction
        # Note: extract_ri_motifs args: --lost, --preserved, --genome_fasta, --outdir
        cmd_ext = [
            "mamba", "run", "-n", "splicing-functional", "python3",
            extract_script,
            "--lost", l_path,
            "--preserved", p_path,
            "--genome_fasta", conf["fasta"],
            "--outdir", dd_dir
        ]
        try:
            run_cmd(cmd_ext, f"Extracting {dirname}")
        except Exception as e:
            print(f"[ERROR] Extraction failed: {e}")
            continue

        # 2. Analysis (AME + MaxEnt + Plots)
        # analyze_ri_vs_constitutive args: --gtf, --genome_fasta, --lost_tsv, --preserved_tsv, --lost_fa/lost_prefix?, --outdir, --script_dir, --motif_db
        
        # Wait, args vary by script.
        # Check analyze_ri... args
        # --lost_fa, --preserved_fa (inputs, relative names or absolute?)
        # For RI: extract produces "lost.intron.fa" etc.
        # Analysis script args: --lost_fa (path to lost.intron.fa)
        
        # Verify RI args
        if etype == "RI":
            lost_fa = os.path.join(dd_dir, "lost.intron.fa")
            pres_fa = os.path.join(dd_dir, "preserved.intron.fa") # Check output naming of extract script.
            # extraction produces: {outdir}/lost.intron.bed/fa
            # yes.
            
            cmd_ana = [
                "mamba", "run", "-n", "splicing-functional", "python3",
                analyze_script,
                "--gtf", conf["gtf"],
                "--genome_fasta", conf["fasta"],
                "--lost_tsv", l_path,
                "--preserved_tsv", p_path,
                "--lost_fa", lost_fa,
                "--preserved_fa", pres_fa,
                "--outdir", dd_dir,
                "--script_dir", SCRIPT_DIR,
                "--motif_db", conf["motif_db"]
            ]
        else:
            # SE Args
            # analyze_se_vs_constitutive args: --lost_prefix, --preserved_prefix (instead of _fa?)
            # Need to check run_pipeline.py for SE invocation.
            # Lines 538: --lost_prefix {step10}/lost
            lost_prefix = os.path.join(dd_dir, "lost")
            pres_prefix = os.path.join(dd_dir, "preserved")
            
            cmd_ana = [
                "mamba", "run", "-n", "splicing-functional", "python3",
                analyze_script,
                "--gtf", conf["gtf"],
                "--genome_fasta", conf['fasta'],
                "--lost_tsv", l_path,
                "--preserved_tsv", p_path,
                "--lost_prefix", lost_prefix,
                "--preserved_prefix", pres_prefix,
                "--outdir", dd_dir,
                "--script_dir", SCRIPT_DIR,
                "--motif_db", conf["motif_db"]
            ]
            
        try:
            run_cmd(cmd_ana, f"Analyzing {dirname}")
        except Exception as e:
            print(f"[ERROR] Analysis failed: {e}")

def main():
    # Scan known directories
    # coverage_test/mouse_RI_20/mouse/total
    # coverage_test/human_RI_20/human/nucleus (if exists)
    
    # Or just scan recursively for step10...?
    # Better to be targeted to avoid scanning huge results if not desired.
    # List:
    targets = [
        "coverage_test/mouse_RI_20/mouse/total",
        "coverage_test/human_RI_20/human/nucleus",
        "coverage_test/mouse_SE_20/mouse/total", # If exists
        # "results/mouse/total", # If fixing everything
        # "results/human/nucleus"
    ]
    
    for t in targets:
        full_p = os.path.join(ROOT_DIR, t)
        if os.path.exists(full_p):
            regenerate_analysis(full_p)
        else:
            print(f"[WARN] Target not found: {full_p}")

if __name__ == "__main__":
    main()
