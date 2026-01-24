import subprocess
import os
import sys
import argparse
import shutil
import glob
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description="Master Pipeline for Figure 4 Splicing Analysis", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--species", required=True, choices=['human', 'mouse', 'arabidopsis'], help="Species to analyze.")
    parser.add_argument("--fraction", choices=["nucleus", "cytosol", "total"], default="nucleus", help="Cellular fraction (default: nucleus)")
    parser.add_argument("--outdir", default="results", help="Main output directory")
    
    # NEW FLAG: Explicitly skip alignment even if FASTA is found
    parser.add_argument("--no-fasta", action="store_true", help="Skip protein alignment even if FASTA is found in data directory.")
    
    parser.add_argument("--dpsi", default="0.1", help="dPSI threshold (NOTE: Currently unused - hardcoded in R script as 0.2 for SE, 0.1 for other event types)")
    parser.add_argument("--min-reads", default="20", help="Minimum reads threshold (applied to sum of replicates per condition; both conditions must pass).")
    parser.add_argument("--normalize", default="log2ratio", help="Normalization method for frame shift density.")
    parser.add_argument("--nperm", default="1000", help="Number of permutations for statistics.")
    parser.add_argument("--nbins", default="5", help="Number of bins for density plots.")
    parser.add_argument("--anchored_window", type=int, default=100, help="Window size (bp) for Start/Stop codon anchored density analysis.")
    parser.add_argument("--pool", action="store_true", help="Pool data before processing.")
    parser.add_argument("--event_types", nargs='+', default=["SE"], help="List of splicing event types to include. Options: SE, A3SS, A5SS, MXE, RI.")
    parser.add_argument("--fdr", type=float, default=0.05, help="FDR threshold for rMATS filtering.")
    parser.add_argument("--fdr_domain", type=float, default=0.01, help="FDR threshold for domain enrichment significance.")
    parser.add_argument("--background", choices=['genome', 'rmats'], default='genome', help="Background for enrichment analysis: 'genome' (all genes) or 'rmats' (only genes tested in rMATS).")
    parser.add_argument("--direction", action="store_true", help="Enable direction-based splitting (dPSI_positive vs dPSI_negative) in AAfeatures.sh.")
    parser.add_argument("--direct", action="store_true", help="Enable directional DeepDive analysis in Step 8 (dPSI_positive and dPSI_negative). Default: only combined and ESE_Analysis.")
    parser.add_argument("--run-mirna", action="store_true", help="Run miRNA analysis in Step 11 (default: skipped).")
    parser.add_argument("--cache-dir", default="data/cache", help="Directory to store persistent cache (TxDb, Biomart results).")
    parser.add_argument("--motif-db", help="Path to MEME Motif DB for AME analysis in Step 10. (default: CisBP for the chosen species)", default=argparse.SUPPRESS)
    parser.add_argument("--motifs", nargs='+', default=None, help="List of specific motifs to plot in RNA Map (Step 12).")
    
    # Execution Control
    parser.add_argument("--start-step", type=int, default=1, help="Start pipeline from this step (1-9).")
    parser.add_argument("--steps", nargs='+', help="Run only specific steps (e.g. 1-7 or 3 4). Overrides --start-step.")
    
    args = parser.parse_args()
    if not hasattr(args, 'motif_db'): args.motif_db = None

    # Parse --steps range if present
    if args.steps:
        expanded_steps = []
        for s in args.steps:
            if '-' in s:
                try:
                    start, end = map(int, s.split('-'))
                    expanded_steps.extend(range(start, end + 1))
                except ValueError:
                    print(f"[ERROR] Invalid range format: {s}")
                    sys.exit(1)
            else:
                try:
                    expanded_steps.append(int(s))
                except ValueError:
                    print(f"[ERROR] Invalid step number: {s}")
                    sys.exit(1)
        args.steps = sorted(list(set(expanded_steps)))
    return args

LOG_DIR = "logs"

def _run_and_log(cmd, log_name):
    """Helper to run shell commands and log status to file and console."""
    global LOG_DIR
    os.makedirs(LOG_DIR, exist_ok=True)
    log_file = os.path.join(LOG_DIR, f"{log_name}.log")
    print(f"[EXEC] {' '.join(cmd)} (Log: {log_file})")
    
    with open(log_file, "w") as f:
        # Use Popen to capture output in real-time
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        
        # Read stdout line by line and write to both
        for line in process.stdout:
            sys.stdout.write(line)
            f.write(line)
            
        process.wait()
        
    if process.returncode != 0:
        print("\n" + "!"*60)
        print(f"       PIPELINE FAILED AT STEP: {log_name}")
        print("!"*60 + "\n")
        print(f"[ERROR] Check detailed log: {log_file}")
        sys.exit(1)

def run_splice_impact(species, indir, outdir, gtf_file, cache_dir, direction):
    """Helper for Step 2: SpliceImpactR Integration."""
    script_path = os.path.join(os.path.dirname(__file__), "src", "get_splice_impact_features.R")
    
    lost_file = os.path.join(indir, "UFM1_dependent.tsv")
    preserved_file = os.path.join(indir, "UFM1_independent.tsv")
    
    cmd = ["Rscript", script_path,
           f"--lost={lost_file}",
           f"--preserved={preserved_file}",
           f"--outdir={outdir}",
           f"--cache_dir={cache_dir}",
           f"--species={species}"]
           
    if gtf_file:
        cmd.append(f"--gtf={gtf_file}")
    
    if direction:
        cmd.append("--direction=TRUE") # Passed but might not be used depending on R script logic
        
    _run_and_log(cmd, "step02_splice_impact_r")

def main():
    args = get_args()
    
    # Set global LOG_DIR to be inside outdir
    global LOG_DIR
    LOG_DIR = os.path.join(args.outdir, "logs")
    os.makedirs(LOG_DIR, exist_ok=True)
    
    project_root = os.path.dirname(os.path.abspath(__file__))
    script_dir = os.path.join(project_root, "src")
    
    # --- File Detection & Setup ---
    data_dir = os.path.join(project_root, "data", args.species)
    # Input Directories
    # Input Directories
    if args.species == 'mouse':
        print("[INFO] Species is Mouse: Hardwiring fraction to 'total'")
        fraction = "total"
        frac_dir = os.path.join(data_dir, fraction)
    elif args.fraction:
        fraction = args.fraction
        frac_dir = os.path.join(data_dir, fraction)
    else:
        # Auto-detect: prefer nucleus if it exists, otherwise check root
        if os.path.isdir(os.path.join(data_dir, "nucleus")):
            fraction = "nucleus"
            frac_dir = os.path.join(data_dir, fraction)
        else:
             fraction = "" # No fraction
             frac_dir = data_dir
    
    # Dynamic detection of WT and UFM directories
    wt_candidates = glob.glob(os.path.join(frac_dir, "[Ww][Tt]*"))
    ufm_candidates = glob.glob(os.path.join(frac_dir, "[Uu][Ff][Mm]*"))
    
    wt_dir = wt_candidates[0] if wt_candidates else os.path.join(frac_dir, "wt")
    ufm_dir = ufm_candidates[0] if ufm_candidates else os.path.join(frac_dir, "ufm")
    
    if not wt_candidates:
        print(f"[WARN] Could not auto-detect WT directory (starting with WT) in {frac_dir}. Defaulting to 'wt'")
    else:
        print(f"[INFO] Detected WT directory: {os.path.basename(wt_dir)}")
        
    if not ufm_candidates:
        print(f"[WARN] Could not auto-detect UFM directory (starting with UFM) in {frac_dir}. Defaulting to 'ufm'")
    else:
        print(f"[INFO] Detected UFM directory: {os.path.basename(ufm_dir)}")

    # Locate GTF
    # Search recursively because sometimes files are in subdirs (e.g. nucleus/)
    gtf_list = glob.glob(os.path.join(data_dir, "**", "*.gtf*"), recursive=True)
    gtf_file = gtf_list[0] if gtf_list else None
    if not gtf_file:
        print(f"[WARN] GTF file not found in {data_dir}. Some steps may fail.")
    else:
        print(f"[INFO] Found GTF: {os.path.basename(gtf_file)}")

    # Locate DNA FASTA
    dna_list = glob.glob(os.path.join(data_dir, "**", "*.dna.*fa"), recursive=True)
    if not dna_list:
        dna_list = glob.glob(os.path.join(data_dir, "**", "*.fa"), recursive=True)
        # Filter out protein fasta if mixed? usually protein is .pep or pc_translations
        dna_list = [f for f in dna_list if "pc_translations" not in f and "pep" not in f]
        
    dna_fasta_file = dna_list[0] if dna_list else None
    if dna_fasta_file: print(f"[INFO] Found DNA FASTA: {os.path.basename(dna_fasta_file)}")

    # Locate Protein FASTA
    prot_list = glob.glob(os.path.join(data_dir, "**", "*.pep.*fa*"), recursive=True)
    if not prot_list:
        prot_list = glob.glob(os.path.join(data_dir, "**", "*pc_translations*fa*"), recursive=True)
        
    # Ensure we don't pick up the index file (.fai)
    if prot_list:
        prot_list = [f for f in prot_list if not f.endswith('.fai')]
        
    prot_fasta_file = prot_list[0] if prot_list else None
    
    if prot_fasta_file: print(f"[INFO] Found Protein FASTA: {os.path.basename(prot_fasta_file)}")

    # --- Set Default Motif DB (CisBP) if not provided ---
    if not args.motif_db:
        motif_map = {
            "human": os.path.join(project_root, "data/motifs/CisBP_Human_All.meme"),
            "mouse": os.path.join(project_root, "data/motifs/CisBP_Mouse_All.meme"),
            "arabidopsis": os.path.join(project_root, "data/motifs/CisBP_Arabidopsis_All.meme")
        }
        
        default_db = motif_map.get(args.species.lower())
        if default_db and os.path.exists(default_db):
             print(f"[INFO] Using Default Motif DB: {os.path.basename(default_db)}")
             args.motif_db = default_db
        else:
             if default_db:
                 print(f"[WARN] Default Motif DB not found at {default_db}. AME steps may be skipped.")
             else:
                 print(f"[WARN] No default motif DB for species '{args.species}'. Provide --motif_db.")
    else:
        print(f"[INFO] Using User-Provided Motif DB: {args.motif_db}")

    # Output Directories
    run_outdir = os.path.join(args.outdir, args.species, fraction)
    
    step1_out = os.path.join(run_outdir, "step01_data_prep")
    step2_out = os.path.join(run_outdir, "step02_domain_enrichment")
    step3_out = os.path.join(run_outdir, "step03_protein_attributes")
    step4_out = os.path.join(run_outdir, "step04_protein_sequence_impact")
    step4_annotated = os.path.join(step4_out, "annotated")  # Functional impact annotations
    step5_out = os.path.join(run_outdir, "step05_frameshift_density")
    step6_out = os.path.join(run_outdir, "step06_aa_features")
    step7_out = os.path.join(run_outdir, "step07_biophysical_properties")
    step8_out = os.path.join(run_outdir, "step08_motif_analysis")
    step9_out = os.path.join(run_outdir, "step09_mechanism_investigation")
    step10_out = os.path.join(run_outdir, "step10_rna_maps")
    step11_out = os.path.join(run_outdir, "step11_genomic_associations")
    
    for d in [step1_out, step2_out, step3_out, step4_out, step5_out, step6_out, step7_out, step8_out, step9_out, step10_out, step11_out]:
        os.makedirs(d, exist_ok=True)
    
    def should_run(step_num):
        if args.steps:
            return step_num in args.steps
        return step_num >= args.start_step

    # --- Execution Ladder ---

    # Step 1: R Preparation
    if should_run(1):
        print("[INFO] === Step 1: Data Preparation ===")
        # Note: wt_dir and ufm_dir detection above might fail if data dir strcuture is different.
        # But we assume standard structure.
        event_types_str = ",".join(args.event_types)
        _run_and_log(["mamba", "run", "-n", "splicing-functional", "Rscript", os.path.join(script_dir, "prepare_rmats_data.R"), 
                        f"--wt_dir={wt_dir}", f"--ufm_dir={ufm_dir}", f"--outdir={step1_out}", 
                        f"--fdr={args.fdr}", f"--dpsi={args.dpsi}", f"--min_reads={args.min_reads}",
                        f"--event_types={event_types_str}"], "step01_data_prep")

        # Automatically export BED/RDS for downstream use (Moved from old Step 14)
        print("[INFO] === Step 1: Exporting Events (BED/RDS) ===")
        cmd_export = ["mamba", "run", "-n", "splicing-functional", "Rscript",
                    os.path.join(script_dir, "export_events.R"),
                    "--dependent", os.path.join(step1_out, "UFM1_dependent.tsv"),
                    "--independent", os.path.join(step1_out, "UFM1_independent.tsv"),
                    "--outdir", step1_out] # Save directly in Step 1 output
        _run_and_log(cmd_export, "step01_export_events")

    # Step 2: Domain Enrichment (R Integration)
    if should_run(2):
        print("[INFO] === Step 2: Extracting Domains & Protein Attributes ===")
        # Note: R script outputs to 'output_dir' (run_outdir here). We move files after generation.
        run_splice_impact(args.species, step1_out, run_outdir, gtf_file, args.cache_dir, args.direction)
        
        # Move files to respective step folders for organization
        domain_tsv_src = os.path.join(run_outdir, 'domain_enrichment.tsv')
        domain_tsv_dst = os.path.join(step2_out, 'domain_enrichment.tsv')
        if os.path.exists(domain_tsv_src):
            os.makedirs(step2_out, exist_ok=True)
            shutil.move(domain_tsv_src, domain_tsv_dst)
            print(f"[INFO] Moved domain_enrichment.tsv to {step2_out}")

        biophys_tsv_src = os.path.join(run_outdir, 'biophysical_enrichment.tsv')
        biophys_tsv_dst = os.path.join(step3_out, 'biophysical_enrichment.tsv')
        if os.path.exists(biophys_tsv_src):
            os.makedirs(step3_out, exist_ok=True)
            shutil.move(biophys_tsv_src, biophys_tsv_dst)
            print(f"[INFO] Moved biophysical_enrichment.tsv to {step3_out}")
    
    # Step 2 Plotting: Domain Enrichment Volcano Plot
    if should_run(2):
        domain_input = os.path.join(step2_out, 'domain_enrichment.tsv')
        if os.path.exists(domain_input):
            print("[INFO] === Step 2: Domain Enrichment Plotting ===")
            enrich_cmd = ["mamba", "run", "-n", "splicing-functional", "python3",
                          os.path.join(script_dir, "analyze_domain_enrichment.py"),
                          f"--enrichment={domain_input}",
                          f"--fdr={args.fdr_domain}",
                          f"--outfile={os.path.join(step2_out, 'Volcano_Domain_Enrichment_Comparison.png')}"]
            _run_and_log(enrich_cmd, "step02_domain_enrichment")

    # Step 3: Protein Attribute Enrichment Plotting
    if should_run(3):
        print("[INFO] === Step 3: Protein Attribute Enrichment Plotting ===")
        biophys_input = os.path.join(step3_out, 'biophysical_enrichment.tsv')
        
        if os.path.exists(biophys_input):
            bio_cmd = ["mamba", "run", "-n", "splicing-functional", "python3",
                       os.path.join(script_dir, "analyze_domain_enrichment.py"),
                       f"--enrichment={biophys_input}",
                       f"--fdr={args.fdr_domain}",
                       "--protein-attributes",
                       f"--outfile={os.path.join(step3_out, 'Protein_Attributes_Enrichment.png')}"]
            _run_and_log(bio_cmd, "step03_protein_attributes")


    # Step 4: Functional Impact Classification (Part 1: Generate Annotations)
    if should_run(4):
        print("[INFO] === Step 4: Functional Impact Classification ===")
        os.makedirs(step4_annotated, exist_ok=True)  # Create annotated subdirectory
        _run_and_log(["mamba", "run", "-n", "splicing-functional", "python3", 
                        os.path.join(script_dir, "splicing_functional_impat.py"), 
                        "--lost_table", os.path.join(step1_out, "UFM1_dependent.tsv"), 
                        "--preserved_table", os.path.join(step1_out, "UFM1_independent.tsv"), 
                        "--gtf", gtf_file, "--outdir", step4_annotated, "--advanced"], "step04_functional_impact")

    # Step 4 (Part 2): Protein Primary Sequence Impact  
    if should_run(4):
        print("[INFO] === Step 4: Protein Sequence Impact Analysis ===")
        
        # Prepare Inputs with optional Directionality Splitting
        inputs_arg = []
        lost_file = os.path.join(step4_annotated, 'UFM1_dependent', 'rmats_sig_annotated.tsv') 
        preserved_file = os.path.join(step4_annotated, 'UFM1_independent', 'rmats_sig_annotated.tsv')
        
        if args.direction:
            print("[INFO] splitting Step 5 outputs by direction (Inc/Exc) for Step 6...")
            
            def split_and_save(infile, out_dir_base):
                if not os.path.exists(infile):
                    print(f"[WARN] Input file {infile} for splitting not found.")
                    return None, None
                
                df = pd.read_csv(infile, sep='\t')
                # Split by IncLevelDifference
                inc = df[df['IncLevelDifference'] > 0]
                exc = df[df['IncLevelDifference'] < 0]
                
                inc_path = os.path.join(out_dir_base, "dPSI_positive.tsv")
                exc_path = os.path.join(out_dir_base, "dPSI_negative.tsv")
                
                os.makedirs(out_dir_base, exist_ok=True)
                inc.to_csv(inc_path, sep='\t', index=False)
                exc.to_csv(exc_path, sep='\t', index=False)
                
                return inc_path, exc_path

            # Split Lost
            l_inc, l_exc = split_and_save(lost_file, os.path.join(step4_annotated, 'lost', 'split'))
            # Split Preserved
            p_inc, p_exc = split_and_save(preserved_file, os.path.join(step4_annotated, 'preserved', 'split'))
            
            if l_inc and l_exc and p_inc and p_exc:
                inputs_arg = [
                    f"UFM1_dependent_dPSI_positive={l_inc}",
                    f"UFM1_dependent_dPSI_negative={l_exc}",
                    f"UFM1_independent_dPSI_positive={p_inc}",
                    f"UFM1_independent_dPSI_negative={p_exc}"
                ]
            else:
                print("[WARN] Splitting failed (files missing?). Reverting to standard UFM1_dependent/UFM1_independent.")
                inputs_arg = [f"UFM1_dependent={lost_file}", f"UFM1_independent={preserved_file}"]
        else:
            inputs_arg = [f"UFM1_dependent={lost_file}", f"UFM1_independent={preserved_file}"]

        # Always run Step 6 to generate per-event tables needed for Step 7.
        cmd6 = ["mamba", "run", "-n", "splicing-functional", "python3", 
                os.path.join(script_dir, "protein_primary_sequence_impact.py"), 
                "--inputs"] + inputs_arg + ["--outdir", step4_out]
        
        if prot_fasta_file and not args.no_fasta:
            print(f"[INFO] Using Protein FASTA for alignment analysis: {prot_fasta_file}")
            cmd6.extend(["--protein_fasta", prot_fasta_file])
        else:
            print("[INFO] Skipping protein alignment analysis (no FASTA or --no-fasta set). Generating other impact tables.")
        
        _run_and_log(cmd6, "step04_seq_impact")

        per_event_file = os.path.join(step4_out, "per_event_compact_for_plotting.tsv")
        if os.path.exists(per_event_file):
            plot_output_file = os.path.join(step4_out, "impact_fractions.png")
            cmd6b = ["mamba", "run", "-n", "splicing-functional", "python3", os.path.join(script_dir, "plot_impact_fractions.py"), "--per_event_file", per_event_file, "--output_file", plot_output_file]
            _run_and_log(cmd6b, "step04_plot_impact_fractions")
            
            # RI Distribution Plot
            if "RI" in args.event_types:
               print("[INFO] Plotting RI Position Distribution with Genomic Background...")
               ri_plot_file = os.path.join(step4_out, "RI_Position_Distribution.png")
               
               # Generate Background
               bg_file = os.path.join(run_outdir, "genome_RI_background_counts.tsv")
               if not os.path.exists(bg_file): # Or always regen?
                   print("[INFO] Generating Genome-wide RI Background...")
                   cmd_bg = ["mamba", "run", "-n", "splicing-functional", "python3", os.path.join(script_dir, "generate_genomic_RI_background.py"), "--gtf", gtf_file, "--out_tsv", bg_file]
                   _run_and_log(cmd_bg, "step04_generate_background")
               
               print("[INFO] Plotting RI Odds Ratios and Saving Stats Table...")
               stats_file = os.path.join(step4_out, "RI_Position_Stats.tsv")
               cmd6c = ["mamba", "run", "-n", "splicing-functional", "python3", os.path.join(script_dir, "plot_RI_distribution.py"), "--per_event_file", per_event_file, "--output_file", ri_plot_file, "--background_counts", bg_file, "--stats_table", stats_file]
               _run_and_log(cmd6c, "step04_plot_RI_distribution")

            # SE Distribution Plot
            if "SE" in args.event_types:
               print("[INFO] Plotting SE Position Distribution with Genomic Background...")
               se_plot_file = os.path.join(step4_out, "SE_Position_Distribution.png")
               
               # Generate Background
               bg_file_se = os.path.join(run_outdir, "genome_SE_background_counts.tsv")
               if not os.path.exists(bg_file_se): 
                   print("[INFO] Generating Genome-wide SE Background (Exons)...")
                   cmd_bg_se = ["mamba", "run", "-n", "splicing-functional", "python3", os.path.join(script_dir, "generate_genomic_SE_background.py"), "--gtf", gtf_file, "--out_tsv", bg_file_se]
                   _run_and_log(cmd_bg_se, "step04_generate_se_background")
               
               print("[INFO] Plotting SE Odds Ratios and Saving Stats Table...")
               stats_file_se = os.path.join(step4_out, "SE_Position_Stats.tsv")
               cmd6d = ["mamba", "run", "-n", "splicing-functional", "python3", os.path.join(script_dir, "plot_SE_distribution.py"), "--per_event_file", per_event_file, "--output_file", se_plot_file, "--background_counts", bg_file_se, "--stats_table", stats_file_se]
               _run_and_log(cmd6d, "step04_plot_SE_distribution")

    # Anchored Density Analysis moved to Step 14 as per user request

    # Step 5: Add Frame Shift Density (Formerly Step 9)
    if should_run(5):
        print("[INFO] === Step 5: Frame Shift Density Analysis ===")
        # Double check directory exists
        os.makedirs(step5_out, exist_ok=True)
        
        input_step7 = os.path.join(step4_out, "per_event_compact_for_plotting.tsv")
        if os.path.exists(input_step7):
            cmd_base = [
                "mamba", "run", "-n", "splicing-functional", "python3", os.path.join(script_dir, "add_frame_shify_density.py"),
                "--per_event", input_step7, "--gtf", gtf_file,
                "--normalize", args.normalize, "--nperm", args.nperm, "--nbins", args.nbins,
                "--bin_stats", "--annotate_binstats", "--overlay_datasets", "--one_plot_per_class",
                "--bin_fisher", "--write_bin_fisher_tsv", "--event_types"] + args.event_types
            if args.pool: cmd_base.append("--pool")

            if args.direction:
                # Run Inclusion Analysis -> dPSI_positive
                print("[INFO] Step 5: Running dPSI_positive (Inclusion) analysis...")
                cmd_inc = cmd_base + ["--out_prefix", f"{step5_out}/fig4_dPSI_positive", "--dataset_order", "UFM1_dependent_dPSI_positive,UFM1_independent_dPSI_positive"]
                _run_and_log(cmd_inc, "step05_frameshift_density_dPSI_positive")
                
                # Run Exclusion Analysis -> dPSI_negative
                print("[INFO] Step 5: Running dPSI_negative (Exclusion) analysis...")
                cmd_exc = cmd_base + ["--out_prefix", f"{step5_out}/fig4_dPSI_negative", "--dataset_order", "UFM1_dependent_dPSI_negative,UFM1_independent_dPSI_negative"]
                _run_and_log(cmd_exc, "step05_frameshift_density_dPSI_negative")
            else:
                # Standard Analysis
                cmd7 = cmd_base + ["--out_prefix", f"{step5_out}/fig4", "--dataset_order", "UFM1_dependent,UFM1_independent"] 
                _run_and_log(cmd7, "step05_frameshift_density")


        else:
            print(f"[WARN] Input file for Step 7 not found: {input_step7}. Skipping Step 7.")

    # Step 6: AA Features from Lost and Preserved (Formerly Step 5)
    if should_run(6):
        if dna_fasta_file:
            print("[INFO] === Step 6: Extracting AA Features ===")
            aa_features_cmd = ["mamba", "run", "-n", "splicing-functional", "bash", 
                               os.path.join(script_dir, "AAfeatures.sh"), 
                               "--lost", os.path.join(step1_out, "UFM1_dependent.tsv"), 
                               "--preserved", os.path.join(step1_out, "UFM1_independent.tsv"), 
                               "-g", dna_fasta_file, 
                               "-o", step6_out,
                               "--bg_gtf", gtf_file,
                               "--bg_cache_dir", os.path.join(args.cache_dir, "bg_exons")]
            if args.direction:
                aa_features_cmd.append("--direction")
            _run_and_log(aa_features_cmd, "step06_AA_features")
        else:
            print("[INFO] Skipping Step 6: AA Features (No DNA FASTA found)")

    # Step 7: Compute Protein Attributes (Biophysical Properties) (Formerly Step 6)
    if should_run(7):
        # Only run if Step 8 was likely possible or if input exists
        # Actually Step 9 takes input from Step 8 output.
        print("[INFO] === Step 7: Computing Protein Attributes (Biophysical Properties) ===")
        biophys_cmd = ["mamba", "run", "-n", "splicing-functional", "python3", 
                       os.path.join(script_dir, "biophysical_properties.py"),
                       "--indir", step6_out,
                       "--out_tsv", os.path.join(step7_out, "protein_attributes_properties.tsv")]
        _run_and_log(biophys_cmd, "step07_protein_attributes_calc")

    # Step 8: RI Motif & Deep Dive Analysis (New Feature)
    if should_run(8):
        if dna_fasta_file:
            print("[INFO] === Step 8: Motif & Deep Dive Analysis ===")
            step8_out = os.path.join(run_outdir, "step08_motif_analysis")
            os.makedirs(step8_out, exist_ok=True)
            
            # Generate SRSF3_SRp20 Seqlogo
            custom_meme = os.path.join(project_root, "data/motifs/custom_SRSF3_SRSF7.meme")
            if os.path.exists(custom_meme):
                print("[INFO] Generating SRSF3_SRp20 Sequence Logo...")
                cmd_logo = ["mamba", "run", "-n", "splicing-functional", "Rscript",
                            os.path.join(script_dir, "plot_custom_logo.R"),
                            custom_meme,
                            os.path.join(step8_out, "SRSF3_SRp20_logo.pdf")]
                _run_and_log(cmd_logo, "step08_srsf3_logo")
            
            # 10A/B: RI Analysis
            # 10A/B: RI Analysis
            if "RI" in args.event_types:
                print("[INFO] === Step 8: RI Motif & Deep Dive Analysis ===")
                
                # Only run directional analyses if --direct flag is set
                if args.direct:
                    print("[INFO] Running directional DeepDive analyses (dPSI_positive and dPSI_negative)...")
                    analyses = [
                         ("combined", "UFM1_dependent.tsv", "UFM1_independent.tsv"),
                         ("dPSI_positive", "UFM1_dependent_dPSI_positive.tsv", "UFM1_independent_dPSI_positive.tsv"),
                         ("dPSI_negative", "UFM1_dependent_dPSI_negative.tsv", "UFM1_independent_dPSI_negative.tsv")
                    ]
                else:
                    print("[INFO] Running combined DeepDive analysis only (use --direct for directional analysis)...")
                    analyses = [
                         ("combined", "UFM1_dependent.tsv", "UFM1_independent.tsv")
                    ]
                
                for name, l_name, p_name in analyses:
                     l_path = os.path.join(step1_out, l_name)
                     p_path = os.path.join(step1_out, p_name)
                     
                     if not os.path.exists(l_path) or not os.path.exists(p_path):
                         print(f"[WARN] Skipping RI sub-analysis '{name}': inputs missing ({l_name}, {p_name}).")
                         continue

                     print(f"[INFO] Running RI Deep Dive: {name}...")
                     step10_ri_out = os.path.join(step8_out, f"RI_DeepDive_{name}")
                     os.makedirs(step10_ri_out, exist_ok=True)
                     
                     # 10A: Extract Sequences
                     cmd_10a = ["mamba", "run", "-n", "splicing-functional", "python3",
                                os.path.join(script_dir, "extract_ri_motifs.py"),
                                "--lost", l_path,
                                "--preserved", p_path,
                                "--genome_fasta", dna_fasta_file,
                                "--outdir", step10_ri_out]
                     _run_and_log(cmd_10a, f"step08_a_extract_ri_{name}")
                     
                     # 10B: Deep Dive
                     lost_fa = os.path.join(step10_ri_out, "UFM1_dependent.intron.fa")
                     pres_fa = os.path.join(step10_ri_out, "UFM1_independent.intron.fa")
                     
                     if os.path.exists(lost_fa) and os.path.exists(pres_fa) and gtf_file:
                          cmd_10b = ["mamba", "run", "-n", "splicing-functional", "python3",
                                    os.path.join(script_dir, "analyze_ri_vs_constitutive.py"),
                                    "--gtf", gtf_file,
                                    "--genome_fasta", dna_fasta_file,
                                    "--lost_tsv", l_path,
                                    "--preserved_tsv", p_path,
                                    "--lost_fa", lost_fa,
                                    "--preserved_fa", pres_fa,
                                    "--outdir", step10_ri_out,
                                    "--script_dir", script_dir]
                          
                          if args.motif_db:
                              cmd_10b.append(f"--motif_db={args.motif_db}")
                              
                          _run_and_log(cmd_10b, f"step08_b_deep_dive_ri_{name}")
                     else:
                          print(f"[WARN] Skipping Step 10B RI {name} (Missing inputs).")

            # 10C: SE Analysis
            if "SE" in args.event_types:
                 print("[INFO] === Step 10C: SE Motif & Deep Dive Analysis ===")
                 
                 # Only run directional analyses if --direct flag is set
                 if args.direct:
                     print("[INFO] Running directional DeepDive analyses (dPSI_positive and dPSI_negative)...")
                     analyses = [
                         ("combined", "UFM1_dependent.tsv", "UFM1_independent.tsv"),
                         ("dPSI_positive", "UFM1_dependent_dPSI_positive.tsv", "UFM1_independent_dPSI_positive.tsv"),
                         ("dPSI_negative", "UFM1_dependent_dPSI_negative.tsv", "UFM1_independent_dPSI_negative.tsv")
                     ]
                 else:
                     print("[INFO] Running combined DeepDive analysis only (use --direct for directional analysis)...")
                     analyses = [
                         ("combined", "UFM1_dependent.tsv", "UFM1_independent.tsv")
                     ]
                 
                 for name, l_name, p_name in analyses:
                     l_path = os.path.join(step1_out, l_name)
                     p_path = os.path.join(step1_out, p_name)
                     
                     if not os.path.exists(l_path) or not os.path.exists(p_path):
                         print(f"[WARN] Skipping SE sub-analysis '{name}': inputs missing.")
                         continue
                         
                     print(f"[INFO] Running SE Deep Dive: {name}...")
                     step10_se_out = os.path.join(step8_out, f"SE_DeepDive_{name}")
                     os.makedirs(step10_se_out, exist_ok=True)
                     
                     # 1. Extract Sequences
                     cmd_se_a = ["mamba", "run", "-n", "splicing-functional", "python3",
                                os.path.join(script_dir, "extract_se_motifs.py"),
                                "--lost", l_path,
                                "--preserved", p_path,
                                "--genome_fasta", dna_fasta_file,
                                "--outdir", step10_se_out]
                     _run_and_log(cmd_se_a, f"step08_a_extract_se_{name}")
                     
                     # 2. Analyze vs Constitutive
                     lost_prefix = os.path.join(step10_se_out, "UFM1_dependent")
                     pres_prefix = os.path.join(step10_se_out, "UFM1_independent")
                     
                     if os.path.exists(lost_prefix + ".exon.fa") and os.path.exists(pres_prefix + ".exon.fa"):
                         cmd_se_b = ["mamba", "run", "-n", "splicing-functional", "python3",
                                    os.path.join(script_dir, "analyze_se_vs_constitutive.py"),
                                    "--gtf", gtf_file,
                                    "--genome_fasta", dna_fasta_file,
                                    "--lost_tsv", l_path,
                                    "--preserved_tsv", p_path,
                                    "--lost_prefix", lost_prefix,
                                    "--preserved_prefix", pres_prefix,
                                    "--outdir", step10_se_out,
                                    "--script_dir", script_dir]
                         
                         if args.motif_db:
                             cmd_se_b.append(f"--motif_db={args.motif_db}")
                             
                         _run_and_log(cmd_se_b, f"step08_b_deep_dive_se_{name}")
                     else:
                         print(f"[WARN] SE FASTA extraction failed or incomplete for {name}. Skipping analysis.")
                     

             
            if "RI" not in args.event_types and "SE" not in args.event_types:
                 print("[INFO] Skipping Step 10 actions (Requires 'RI' or 'SE' in --event_types).")
            
            # Step 10D: Comprehensive ESE/RBP/AME Analysis
            if "RI" in args.event_types:
                print("[INFO] === Step 10D: ESE/RBP Motif Analysis with AME ===")
                step10_ese_out = os.path.join(step8_out, "ESE_Analysis")
                os.makedirs(step10_ese_out, exist_ok=True)
                
                # Get motif database
                motif_db = getattr(args, 'motif_db', os.path.join(project_root, "data/motifs/CisBP_Human_All.meme"))
                if args.species.lower() == "mouse":
                    motif_db = os.path.join(project_root, "data/motifs/CisBP_Mouse_All.meme")
                elif args.species.lower() == "arabidopsis":
                    motif_db = os.path.join(project_root, "data/motifs/CisBP_Arabidopsis_All.meme")
                
                ese_cmd = ["mamba", "run", "-n", "splicing-functional", "python3",
                           os.path.join(script_dir, "run_ese_analysis.py"),
                           "--events_rds", os.path.join(step1_out, "UFM1_events_rich.rds"),
                           "--genome_fasta", dna_fasta_file,
                           "--gtf", gtf_file,
                           "--outdir", step10_ese_out,
                           "--motif_db", motif_db]
                _run_and_log(ese_cmd, "step10d_ese_analysis")
        else:
             print("[INFO] Skipping Step 10 (Requires DNA FASTA).")

             
    # Step 9: Mechanism Investigation (Stalling & Adjacency)
    if should_run(9):
        print("[INFO] === Step 9: Mechanism Investigation (Stalling & Adjacency) ===")
        step9_out = os.path.join(run_outdir, "step09_mechanism_investigation")
        os.makedirs(step9_out, exist_ok=True)
        step8_out = os.path.join(run_outdir, "step08_motif_analysis")
        
        # 11A & 11C: Stalling Motifs & Seq Properties
        # Loop through directions if available: neg (Lost) and pos (Preserved but maybe different context?)
        # Stalling is usually relevant for "Lost" events (did gaining the intron cause stalling? or losing it?)
        # User wants split for "psI" (pos?) and "psE" (neg?).
        
        directions_to_process = ["dPSI_negative", "dPSI_positive"] if args.direction else ["dPSI_negative"] # Default to negative?
        # Actually user ran with --direction.
        # Step 10 produced SE_DeepDive_neg and SE_DeepDive_pos.
        
        for direction in directions_to_process:
             # Identify input directory from Step 10
             # Note: For RI, it is RI_DeepDive_... For SE, SE_DeepDive_...
             # We need to detect which event type ran.
             # If both, prioritize SE for Stalling? Or loop both?
             # User's current failure is on RI run.
             
             target_dirs = []
             if "SE" in args.event_types: target_dirs.append(("SE", os.path.join(step8_out, f"SE_DeepDive_{direction}")))
             if "RI" in args.event_types: target_dirs.append(("RI", os.path.join(step8_out, f"RI_DeepDive_{direction}")))
             
             for etype, deep_dir in target_dirs:
                 if not os.path.exists(deep_dir): continue
                 
                 print(f"[INFO] Step 11A ({etype}_{direction}): Running Stalling Motif Analysis...")
                 
                 # Define inputs
                 # For SE: lost.exon.fa vs preserved.exon.fa
                 # For RI: lost.intron.fa vs preserved.intron.fa (maybe? Stalling in intron?)
                 # RI stalling usually checks Retained Intron sequence (Preserved Intron).
                 # So "Preserved" in RI means Intron is kept. "Lost" means spliced out.
                 # We want to check sequences of Preserved (Retained) vs Lost (Spliced/Constitutive-like?).
                 # extract_ri_motifs produces lost.intron.fa and preserved.intron.fa
                 
                 if etype == "SE":
                     input_cands = [
                         ("UFM1_dependent", os.path.join(deep_dir, "UFM1_dependent.exon.fa")), 
                         ("UFM1_independent", os.path.join(deep_dir, "UFM1_independent.exon.fa")),
                         ("Genome", os.path.join(deep_dir, "constitutive_exons.exon.fa"))
                     ]
                 else:
                     input_cands = [
                         ("UFM1_dependent", os.path.join(deep_dir, "UFM1_dependent.intron.fa")), 
                         ("UFM1_independent", os.path.join(deep_dir, "UFM1_independent.intron.fa")),
                         ("Genome", os.path.join(deep_dir, "constitutive_introns.intron.fa"))
                     ]
                 
                 # Verify files exist
                 valid_groups = []
                 for name, fpath in input_cands:
                     if os.path.exists(fpath): valid_groups.append(f"{name}:{fpath}")
                 
                 if valid_groups:
                     out_subdir = os.path.join(step9_out, f"{etype}_{direction}")
                     os.makedirs(out_subdir, exist_ok=True)
                     
                     # 1. Stalling
                     cmd_11a = ["mamba", "run", "-n", "splicing-functional", "python3", 
                                os.path.join(script_dir, "analyze_stalling_motifs.py"),
                                "--groups"] + valid_groups + ["--outdir", out_subdir]
                     _run_and_log(cmd_11a, f"step09_a_stalling_{etype}_{direction}")
                     
                     # 2. GC/Len Check
                     # check_gc_len needs explicit lost/pres inputs
                     if len(input_cands) >= 2 and os.path.exists(input_cands[0][1]) and os.path.exists(input_cands[1][1]):
                          print(f"[INFO] Step 11C ({etype}_{direction}): Running Sequence Property Check...")
                          cmd_11c = ["mamba", "run", "-n", "splicing-functional", "python3",
                                     os.path.join(script_dir, "check_gc_len.py"),
                                     "--lost", input_cands[0][1],
                                     "--preserved", input_cands[1][1]]
                          _run_and_log(cmd_11c, f"step09_c_seq_properties_{etype}_{direction}")
                 else:
                     print(f"[WARN] inputs missing for {etype}_{direction} in {deep_dir}")

        # 11B: Adjacency (SE-RI) - Unchanged as it depends on Step 1 tables which are already split if needed
        # But analyze_se_ri_adjacency.py processes the full table usually.
        # User output showed Step 11B ran.
        l_path = os.path.join(step1_out, "UFM1_dependent.tsv")
        p_path = os.path.join(step1_out, "UFM1_independent.tsv")
        if os.path.exists(l_path):
             print("[INFO] Step 11B: Running SE-RI Adjacency Analysis...")
             adj_out = os.path.join(step9_out, "adjacency")
             os.makedirs(adj_out, exist_ok=True)
             cmd_11b = ["mamba", "run", "-n", "splicing-functional", "python3", 
                        os.path.join(script_dir, "analyze_se_ri_adjacency.py"),
                        "--se_files", l_path, p_path,
                        "--ri_files", l_path, p_path,
                        "--outdir", adj_out]
             _run_and_log(cmd_11b, "step09_b_adjacency")

    # Step 10: Positional Motif Enrichment (RNA Map)
    if should_run(10):
        if dna_fasta_file and args.motif_db:
             print("[INFO] === Step 10: Positional Motif Enrichment (RNA Map) ===")
             step10_out = os.path.join(run_outdir, "step10_rna_maps")
             os.makedirs(step10_out, exist_ok=True)
             
             # Inputs:
             # We prefer dPSI_negative (Exclusion/Loss) events as they are "Lost".
             # If direction splits exist:
             lost_path = os.path.join(step1_out, "UFM1_dependent_dPSI_negative.tsv")
             preserved_path = os.path.join(step1_out, "UFM1_independent_dPSI_negative.tsv")
             
             if not os.path.exists(lost_path):
                 # Fallback to combined if split not found
                 lost_path = os.path.join(step1_out, "UFM1_dependent.tsv")
                 preserved_path = os.path.join(step1_out, "UFM1_independent.tsv")
                 print(f"[INFO] Using combined files for RNA Map (dPSI_negative not found).")
             
             # Constitutive Background:
             # From Step 10 - check RI or SE directories
             # We prioritize RI as per current focus.
             constitutive_bed = None
             if "RI" in args.event_types:
                 # Check Step 10 RI Negative dir
                 c_cand = os.path.join(run_outdir, "step08_motif_analysis", "RI_DeepDive_dPSI_negative", "constitutive_introns.intron.bed")
                 if os.path.exists(c_cand): constitutive_bed = c_cand
             
             if not constitutive_bed and "SE" in args.event_types:
                 # Check SE
                 c_cand = os.path.join(run_outdir, "step08_motif_analysis", "SE_DeepDive_dPSI_negative", "constitutive_exons.exon.bed") 
                 # Wait, RNA map usually specific to event type. Mixing RI/SE might be weird.
                 # User requested RNA Map for RI specifically in context of "Intron Retention motfis".
                 pass

             if os.path.exists(lost_path) and os.path.exists(preserved_path):
                 cmd_12 = ["mamba", "run", "-n", "splicing-functional", "python3",
                           os.path.join(script_dir, "plot_rna_map.py"),
                           "--lost", lost_path,
                           "--preserved", preserved_path,
                           "--genome", dna_fasta_file,
                           "--motif_db", args.motif_db,
                           "--motifs", *(args.motifs if args.motifs else ["RBMY1D", "MBNL1", "PCBP2", "PCBP3", "SRSF5", "SRSF3", "SRSF3_SRp20"]),
                           "--outdir", step10_out,
                           "--window", "50",
                           "--smooth", "40"]
                 
                 if constitutive_bed:
                     cmd_12.extend(["--constitutive", constitutive_bed])
                 else:
                     print("[WARN] Constitutive background BED not found. Plotting without background.")
                     
                 _run_and_log(cmd_12, "step10_rna_map")
             else:
                 print(f"[WARN] Skipping Step 10: Input files missing ({lost_path}, {preserved_path})")
        else:
             print("[INFO] Skipping Step 12 (Requires DNA FASTA and Motif DB).")
             
    # Step 11: Anchored Intron GC Content Analysis
    # Merged into Step 14
    if should_run(11):
        print("[INFO] Step 13 is merged into Step 14. Please run Step 14.")
    # Step 11: Genomic & Sequence Associations (miRNA, NMD) - Uses files from Step 1
    if should_run(11):
         print("[INFO] === Step 11: Genomic & Sequence Associations ===")
         
         step11_out = os.path.join(run_outdir, "step11_genomic_associations")
         os.makedirs(step11_out, exist_ok=True)
         
         # Check Step 1 for the RDS file
         events_rds = os.path.join(step1_out, "UFM1_events_rich.rds")
         
         if not os.path.exists(events_rds):
             print(f"[WARN] Step 14 Input {events_rds} missing! Did Step 1 run? Attempting to generate it now...")
             # Fallback: Run export if missing
             cmd_14 = ["mamba", "run", "-n", "splicing-functional", "Rscript",
                    os.path.join(script_dir, "export_events.R"),
                    "--dependent", os.path.join(step1_out, "UFM1_dependent.tsv"),
                    "--independent", os.path.join(step1_out, "UFM1_independent.tsv"),
                    "--outdir", step1_out]
             _run_and_log(cmd_14, "step11_export_events_fallback")
         
         # Genomic Associations (Merged from Step 15)
         print("[INFO] === Genomic Associations (miRNA, NMD) ===")
         
         # Note: Step 1 now outputs to step1_out, but we can organize them into step11_out if preferred.
         # For now, let's keep outputs in step11_out but use input from step1_out.
         if os.path.exists(events_rds) and gtf_file and dna_fasta_file:
             # A. miRNA Isoform Analysis (Optional - only run if --run-mirna is set)
             if args.run_mirna:
                 print("[INFO] Step 11A: miRNA Isoform Analysis...")
                 mirna_out = os.path.join(step11_out, "mirna_isoform_analysis")
                 os.makedirs(mirna_out, exist_ok=True)
                 
                 cmd_14a = ["mamba", "run", "-n", "splicing-functional", "Rscript",
                            os.path.join(script_dir, "analyze_mirna_isoforms.R"),
                            "--gtf", gtf_file,
                            "--fasta", dna_fasta_file,
                            "--events", events_rds,
                            "--species", args.species,
                            "--outdir", mirna_out]
                 _run_and_log(cmd_14a, "step11a_mirna_isoforms")
             else:
                 print("[INFO] Skipping Step 11A: miRNA Isoform Analysis (use --run-mirna to enable)")
             
             # B. NMD Analysis (Features & Metagene)
             print("[INFO] Step 14B: NMD Susceptibility Analysis...")
             nmd_out = os.path.join(step11_out, "nmd_analysis")
             os.makedirs(nmd_out, exist_ok=True)
             
             # Risk & Length
             cmd_14b1 = ["mamba", "run", "-n", "splicing-functional", "Rscript",
                         os.path.join(script_dir, "analyze_nmd_features.R"),
                         "--gtf", gtf_file,
                         "--events", events_rds,
                         "--outdir", nmd_out]
             _run_and_log(cmd_14b1, "step14b_nmd_features")
             
             # Metagene Plot
             cmd_14b2 = ["mamba", "run", "-n", "splicing-functional", "Rscript",
                         os.path.join(script_dir, "plot_nmd_metagene.R"),
                         "--gtf", gtf_file,
                         "--events", events_rds,
                         "--outdir", nmd_out]
             _run_and_log(cmd_14b2, "step14b_nmd_metagene")

             # Feature Length Analysis (Introns & 3'UTRs)
             print("[INFO] Step 14B-3: Feature Length Analysis...")
             feat_len_out = os.path.join(step11_out, "feature_lengths")
             os.makedirs(feat_len_out, exist_ok=True)
             
             cmd_14b3 = ["mamba", "run", "-n", "splicing-functional", "Rscript",
                         os.path.join(script_dir, "analyze_feature_lengths.R"),
                         "--gtf", gtf_file,
                         "--events", events_rds,
                         "--outdir", feat_len_out]
             _run_and_log(cmd_14b3, "step14b_feature_lengths")
             
             # EJC + PTC Analysis
             print("[INFO] Step 14B-4: EJC & PTC Analysis...")
             cmd_14b4 = ["mamba", "run", "-n", "splicing-functional", "Rscript",
                         os.path.join(script_dir, "EJC_PTC.R"),
                         "--gtf", gtf_file,
                         "--fasta", dna_fasta_file,
                         "--events", events_rds,
                         "--outdir", nmd_out]
             _run_and_log(cmd_14b4, "step14b_ejc_ptc")
             
             # Anchored Stop Codon Intron Density
             print("[INFO] Step 14B-5: Anchored Stop Codon Intron Density...")
             cmd_14b5 = ["mamba", "run", "-n", "splicing-functional", "Rscript",
                         os.path.join(script_dir, "anchored_stopcodon_intron_density.R"),
                         "--gtf", gtf_file,
                         "--events", events_rds,
                         "--outdir", nmd_out,
                         "--window", str(args.anchored_window)]
             _run_and_log(cmd_14b5, "step14b_stopcodon_density")
             
             # Step 14B-6: 3'SS Sequence Logo (+/- 20bp)
             print("[INFO] Step 14B-6: 3'SS Sequence Logo (+/- 20bp)...")
             step14_3ss_out = os.path.join(step11_out, "3ss_sequence_logos")
             os.makedirs(step14_3ss_out, exist_ok=True)
             
             cmd_14b6_extract = ["mamba", "run", "-n", "splicing-functional", "python3",
                                os.path.join(script_dir, "extract_3ss_sequences.py"),
                                "--lost", os.path.join(step1_out, "UFM1_dependent.tsv"),
                                "--preserved", os.path.join(step1_out, "UFM1_independent.tsv"),
                                "--gtf", gtf_file,
                                "--genome_fasta", dna_fasta_file,
                                "--outdir", step14_3ss_out,
                                "--script_dir", script_dir]
             _run_and_log(cmd_14b6_extract, "step14b_extract_3ss")

             cmd_14b6_plot = ["mamba", "run", "-n", "splicing-functional", "Rscript",
                             os.path.join(script_dir, "plot_3ss_logo.R"),
                             "--dep", os.path.join(step14_3ss_out, "UFM1_dependent.3ss_20bp.fa"),
                             "--indep", os.path.join(step14_3ss_out, "UFM1_independent.3ss_20bp.fa"),
                             "--const", os.path.join(step14_3ss_out, "Constitutive.3ss_20bp.fa"),
                             "--out", os.path.join(step14_3ss_out, "3ss_sequence_logo.pdf")]
             _run_and_log(cmd_14b6_plot, "step14b_plot_3ss_logo")
             
         else:
             print("[WARN] Skipping Genomic Associations: Inputs missing (RDS, GTF, or FASTA). Run Step 14 Export first.")

         # C. Anchored Intron GC (Merged from Step 13)
         if dna_fasta_file:
             print("[INFO] === Step 14C: Anchored Intron GC Content Analysis (Merged Step 13) ===")
             step14_gc_out = os.path.join(step11_out, "anchored_intron_gc")
             os.makedirs(step14_gc_out, exist_ok=True)

             dependent_tsv = os.path.join(step1_out, "UFM1_dependent.tsv")
             independent_tsv = os.path.join(step1_out, "UFM1_independent.tsv")

             if os.path.exists(dependent_tsv) and os.path.exists(independent_tsv):
                 cmd_gc = ["mamba", "run", "-n", "splicing-functional", "python3",
                           os.path.join(script_dir, "analyze_anchored_intron_gc.py"),
                           "--dependent", dependent_tsv,
                           "--independent", independent_tsv,
                           "--gtf", gtf_file,
                           "--genome", dna_fasta_file,
                           "--outdir", step14_gc_out,
                           "--window", str(args.anchored_window)]
                 _run_and_log(cmd_gc, "step14c_anchored_intron_gc")
             else:
                 print("[WARN] Skipping Anchored GC: Input TSVs missing.")

         # D. Anchored Density Plot Data (Merged from Step 6)
         print("[INFO] === Step 14D: Anchored Density Analysis (Merged from Step 6) ===")
         step6_out_dir = os.path.join(run_outdir, "step04_protein_sequence_impact")
         per_event_file_local = os.path.join(step6_out_dir, "per_event_compact_for_plotting.tsv")
         
         if os.path.exists(per_event_file_local):
             event_type_str = "_".join(args.event_types)
             out_filename = f"anchored_density_{event_type_str}_plot_data.tsv"
             out_filepath = os.path.join(step11_out, out_filename) 
             
             cmd_dens = ["mamba", "run", "-n", "splicing-functional", "python3",
                         os.path.join(script_dir, "analyze_anchored_density.py"),
                         "--input_table", per_event_file_local,
                         "--gtf", gtf_file,
                         "--outdir", step11_out,
                         "--outfile", out_filepath,
                         "--window", str(args.anchored_window)]
             _run_and_log(cmd_dens, "step14d_anchored_density")
         else:
             print(f"[WARN] Skipping Anchored Density: Step 6 Output ({per_event_file_local}) missing.")


    # ===========================================================================
    # STEP 12: PROTEIN PROPERTY METAGENE ANALYSIS (All Species)
    # ===========================================================================
    if should_run(12):
        print("\n[INFO] === Step 12: Protein Property Metagene Analysis ===")
        step12_out = os.path.join(run_outdir, "step12_protein_metagene")
        os.makedirs(step12_out, exist_ok=True)
        
        dep_tsv = os.path.join(step1_out, "UFM1_dependent.tsv")
        indep_tsv = os.path.join(step1_out, "UFM1_independent.tsv")
        
        if os.path.exists(dep_tsv) and os.path.exists(indep_tsv) and gtf_file and dna_fasta_file:
            print(f"[INFO] Running protein metagene for {args.species}...")
            import tempfile, subprocess
            
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_dep.txt') as f:
                dep_genes_file = f.name
                subprocess.run(f"awk -F'\\t' '{{print $2}}' {dep_tsv} | tail -n +2 | sort -u > {dep_genes_file}", shell=True, check=True)
            
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_indep.txt') as f:
                indep_genes_file = f.name
                subprocess.run(f"awk -F'\\t' '{{print $2}}' {indep_tsv} | tail -n +2 | sort -u > {indep_genes_file}", shell=True, check=True)
            
            cmd_12 = ["mamba", "run", "-n", "splicing-functional", "python",
                     os.path.join(script_dir, "protein_properties_metagene.py"),
                     "--gtf", gtf_file, "--genome", dna_fasta_file,
                     "--dep_genes", dep_genes_file, "--indep_genes", indep_genes_file,
                     "--outdir", step12_out, "--window", "30", "--step", "10", "--nbins", "100"]
            _run_and_log(cmd_12, "step12_protein_metagene")
            
            os.unlink(dep_genes_file); os.unlink(indep_genes_file)
        else:
            print("[WARN] Skipping Protein Metagene: Missing files")


    # ===========================================================================
    # STEP 13: SUBCELLULAR LOCALIZATION ENRICHMENT (All Species)
    # ===========================================================================
    if should_run(13):
        print("\n[INFO] === Step 13: Subcellular Localization Enrichment ===")
        step13_out = os.path.join(run_outdir, "step13_subcellular_localization")
        os.makedirs(step13_out, exist_ok=True)
        
        events_rds = os.path.join(step1_out, "UFM1_events_rich.rds")
        
        if os.path.exists(events_rds):
            print(f"[INFO] Running subcellular localization for {args.species}...")
            cmd_13 = ["mamba", "run", "-n", "splicing-functional", "Rscript",
                     os.path.join(script_dir, "analyze_subcellular_localization.R"),
                     "--species", args.species,
                     "--rds", events_rds,
                     "--outdir", step13_out]
            _run_and_log(cmd_13, "step13_subcellular_localization")
        else:
            print(f"[WARN] Skipping Step 13: Missing RDS file ({events_rds})")


    # ===========================================================================
    # STEP 14: GO ENRICHMENT ANALYSIS (BP & MF)
    # ===========================================================================
    if should_run(14):
        print("\n[INFO] === Step 14: GO Enrichment Analysis (BP & MF) ===")
        step14_out = os.path.join(run_outdir, "step14_go_enrichment")
        os.makedirs(step14_out, exist_ok=True)
        
        events_rds = os.path.join(step1_out, "UFM1_events_rich.rds")
        
        if os.path.exists(events_rds):
            print(f"[INFO] Running GO enrichment for {args.species}...")
            cmd_14 = ["mamba", "run", "-n", "splicing-functional", "Rscript",
                     os.path.join(script_dir, "analyze_go_enrichment.R"),
                     "--species", args.species,
                     "--rds", events_rds,
                     "--outdir", step14_out]
            _run_and_log(cmd_14, "step14_go_enrichment")
        else:
            print(f"[WARN] Skipping Step 14: Missing RDS file ({events_rds})")

    # ===========================================================================
    # STEP 15: CROSS-SPECIES CONSERVATION ANALYSIS (Eukaryotic GO Landscape)
    # ===========================================================================
    if should_run(15):
        print("\n[INFO] === Step 15: Cross-Species Conservation Analysis ===")
        viz_out = os.path.join(args.outdir, "Cross_species_conservation")
        
        # This script performs the refined Eukaryotic Landscape analysis
        # capturing conserved functional themes across Human, Mouse, and Arabidopsis.
        cmd_15 = ["mamba", "run", "-n", "splicing-functional", "Rscript",
                 os.path.join(script_dir, "plot_eukaryotic_go_landscape.R"),
                 "--outdir", viz_out,
                 "--results_dir", args.outdir]
        
        _run_and_log(cmd_15, "step15_cross_species_conservation")

    print("       PIPELINE COMPLETED SUCCESSFULLY")
    print("="*60 + "\n")
    print(f"Results are available in: {run_outdir}")

if __name__ == "__main__":
    main()
