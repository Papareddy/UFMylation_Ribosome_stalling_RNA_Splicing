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
    
    parser.add_argument("--dpsi", default="0.15", help="Legacy argument. NOTE: dPSI thresholds are now event-specific (SE=0.2, Others=0.1) and set internally in Step 1.")
    parser.add_argument("--min-reads", default="20", help="Minimum reads threshold.")
    parser.add_argument("--normalize", default="log2ratio", help="Normalization method for frame shift density.")
    parser.add_argument("--nperm", default="1000", help="Number of permutations for statistics.")
    parser.add_argument("--nbins", default="5", help="Number of bins for density plots.")
    parser.add_argument("--pool", action="store_true", help="Pool data before processing.")
    parser.add_argument("--event_types", nargs='+', default=["SE"], help="List of splicing event types to include. Options: SE, A3SS, A5SS, MXE, RI.")
    parser.add_argument("--fdr", type=float, default=0.05, help="FDR threshold for rMATS filtering.")
    parser.add_argument("--fdr_domain", type=float, default=0.01, help="FDR threshold for domain enrichment significance.")
    parser.add_argument("--background", choices=['genome', 'rmats'], default='genome', help="Background for enrichment analysis: 'genome' (all genes) or 'rmats' (only genes tested in rMATS).")
    parser.add_argument("--direction", action="store_true", help="Enable direction-based splitting (dPSI_positive vs dPSI_negative) in AAfeatures.sh.")
    parser.add_argument("--cache-dir", default="data/cache", help="Directory to store persistent cache (TxDb, Biomart results).")
    parser.add_argument("--motif-db", help="Path to MEME Motif DB for AME analysis in Step 10. (default: CisBP for the chosen species)", default=argparse.SUPPRESS)
    
    # Execution Control
    parser.add_argument("--start-step", type=int, default=1, help="Start pipeline from this step (1-9).")
    parser.add_argument("--steps", type=int, nargs='+', help="Run only specific steps (e.g. 3 4). Overrides --start-step.")
    
    args = parser.parse_args()
    if not hasattr(args, 'motif_db'): args.motif_db = None
    return args

def _run_and_log(cmd, log_name):
    """Helper to run shell commands and log status to file and console."""
    os.makedirs("logs", exist_ok=True)
    log_file = os.path.join("logs", f"{log_name}.log")
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
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # --- File Detection & Setup ---
    data_dir = os.path.join(script_dir, "data", args.species)
    # Input Directories
    if args.fraction:
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
            "human": os.path.join(script_dir, "data/motifs/CisBP_Human_All.meme"),
            "mouse": os.path.join(script_dir, "data/motifs/CisBP_Mouse_All.meme"),
            "arabidopsis": os.path.join(script_dir, "data/motifs/CisBP_Arabidopsis_All.meme")
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
    step2_out = os.path.join(run_outdir, "step02_splice_impact")
    step3_out = os.path.join(run_outdir, "step03_domain_enrichment")
    step4_out = os.path.join(run_outdir, "step04_protein_attributes")
    step5_out = os.path.join(run_outdir, "step05_functional_impact")       # Old Step 7
    step6_out = os.path.join(run_outdir, "step06_protein_sequence_impact") # Old Step 8
    step7_out = os.path.join(run_outdir, "step07_frameshift_density")      # Old Step 9
    step8_out = os.path.join(run_outdir, "step08_aa_features")             # Old Step 5
    step9_out = os.path.join(run_outdir, "step09_biophysical_properties")  # Old Step 6(renamed)
    
    for d in [step1_out, step2_out, step3_out, step4_out, step5_out, step6_out, step7_out, step8_out, step9_out]:
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
        _run_and_log(["mamba", "run", "-n", "splicing-functional", "Rscript", os.path.join(script_dir, "src/prepare_rmats_data.R"), 
                        f"--wt_dir={wt_dir}", f"--ufm_dir={ufm_dir}", f"--outdir={step1_out}", 
                        f"--fdr={args.fdr}", f"--dpsi={args.dpsi}", f"--min_reads={args.min_reads}",
                        f"--event_types={event_types_str}"], "step01_data_prep")

    # Step 2: SpliceImpactR Integration
    if should_run(2):
        print("[INFO] === Step 2: Extracting Domains & Protein Attributes ===")
        # Note: R script outputs to 'output_dir' (run_outdir here). We move files after generation.
        run_splice_impact(args.species, step1_out, run_outdir, gtf_file, args.cache_dir, args.direction)
        
        # Move files to respective step folders for organization
        domain_tsv_src = os.path.join(run_outdir, 'domain_enrichment.tsv')
        domain_tsv_dst = os.path.join(step3_out, 'domain_enrichment.tsv')
        if os.path.exists(domain_tsv_src):
            shutil.move(domain_tsv_src, domain_tsv_dst)
            print(f"[INFO] Moved domain_enrichment.tsv to {step3_out}")

        biophys_tsv_src = os.path.join(run_outdir, 'biophysical_enrichment.tsv')
        biophys_tsv_dst = os.path.join(step4_out, 'biophysical_enrichment.tsv')
        if os.path.exists(biophys_tsv_src):
            shutil.move(biophys_tsv_src, biophys_tsv_dst)
            print(f"[INFO] Moved biophysical_enrichment.tsv to {step4_out}")
    
    # Step 3: Domain Enrichment Plotting
    if should_run(3):
        print("[INFO] === Step 3: Domain Enrichment Plotting ===")
        domain_input = os.path.join(step3_out, 'domain_enrichment.tsv')
        if not os.path.exists(domain_input):
             print(f"[WARN] Step 3 Input {domain_input} missing! Did you run Step 2?")
        
        enrich_cmd = ["mamba", "run", "-n", "splicing-functional", "python3",
                      os.path.join(script_dir, "src/analyze_domain_enrichment.py"),
                      f"--enrichment={domain_input}",
                      f"--fdr={args.fdr_domain}",
                      f"--outfile={os.path.join(step3_out, 'Volcano_Domain_Enrichment_Comparison.png')}"]
        _run_and_log(enrich_cmd, "step03_domain_enrichment")

    # Step 4: Protein Attribute Enrichment Plotting
    if should_run(4):
        print("[INFO] === Step 4: Protein Attribute Enrichment Plotting ===")
        biophys_input = os.path.join(step4_out, 'biophysical_enrichment.tsv')
        
        bio_cmd = ["mamba", "run", "-n", "splicing-functional", "python3",
                   os.path.join(script_dir, "src/analyze_domain_enrichment.py"),
                   f"--enrichment={biophys_input}",
                   f"--fdr={args.fdr_domain}",
                   "--protein-attributes",
                   f"--outfile={os.path.join(step4_out, 'Protein_Attributes_Enrichment.png')}"]
        _run_and_log(bio_cmd, "step04_protein_attributes")


    # Step 5: Splicing Functional Impact (Formerly Step 7)
    if should_run(5):
        print("[INFO] === Step 5: Functional Impact Classification ===")
        _run_and_log(["mamba", "run", "-n", "splicing-functional", "python3", 
                        os.path.join(script_dir, "src/splicing_functional_impat.py"), 
                        "--lost_table", os.path.join(step1_out, "UFM1_dependent.tsv"), 
                        "--preserved_table", os.path.join(step1_out, "UFM1_independent.tsv"), 
                        "--gtf", gtf_file, "--outdir", step5_out, "--advanced"], "step05_functional_impact")

    # Step 6: Protein Primary Sequence Impact (Formerly Step 8)
    if should_run(6):
        print("[INFO] === Step 6: Protein Sequence Impact ===")
        
        # Prepare Inputs with optional Directionality Splitting
        inputs_arg = []
        lost_file = os.path.join(step5_out, 'UFM1_dependent', 'rmats_sig_annotated.tsv') 
        preserved_file = os.path.join(step5_out, 'UFM1_independent', 'rmats_sig_annotated.tsv')
        
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
            l_inc, l_exc = split_and_save(lost_file, os.path.join(step5_out, 'lost', 'split'))
            # Split Preserved
            p_inc, p_exc = split_and_save(preserved_file, os.path.join(step5_out, 'preserved', 'split'))
            
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
                os.path.join(script_dir, "src/protein_primary_sequence_impact.py"), 
                "--inputs"] + inputs_arg + ["--outdir", step6_out]
        
        if prot_fasta_file and not args.no_fasta:
            print(f"[INFO] Using Protein FASTA for alignment analysis: {prot_fasta_file}")
            cmd6.extend(["--protein_fasta", prot_fasta_file])
        else:
            print("[INFO] Skipping protein alignment analysis (no FASTA or --no-fasta set). Generating other impact tables.")
        
        _run_and_log(cmd6, "step06_seq_impact")

        per_event_file = os.path.join(step6_out, "per_event_compact_for_plotting.tsv")
        if os.path.exists(per_event_file):
            plot_output_file = os.path.join(step6_out, "impact_fractions.png")
            cmd6b = ["mamba", "run", "-n", "splicing-functional", "python3", os.path.join(script_dir, "src/plot_impact_fractions.py"), "--per_event_file", per_event_file, "--output_file", plot_output_file]
            _run_and_log(cmd6b, "step06_plot_impact_fractions")
            
            # RI Distribution Plot
            if "RI" in args.event_types:
               print("[INFO] Plotting RI Position Distribution with Genomic Background...")
               ri_plot_file = os.path.join(step6_out, "RI_Position_Distribution.png")
               
               # Generate Background
               bg_file = os.path.join(run_outdir, "genome_RI_background_counts.tsv")
               if not os.path.exists(bg_file): # Or always regen?
                   print("[INFO] Generating Genome-wide RI Background...")
                   cmd_bg = ["mamba", "run", "-n", "splicing-functional", "python3", os.path.join(script_dir, "src/generate_genomic_RI_background.py"), "--gtf", gtf_file, "--out_tsv", bg_file]
                   _run_and_log(cmd_bg, "step06_generate_background")
               
               print("[INFO] Plotting RI Odds Ratios and Saving Stats Table...")
               stats_file = os.path.join(step6_out, "RI_Position_Stats.tsv")
               cmd6c = ["mamba", "run", "-n", "splicing-functional", "python3", os.path.join(script_dir, "src/plot_RI_distribution.py"), "--per_event_file", per_event_file, "--output_file", ri_plot_file, "--background_counts", bg_file, "--stats_table", stats_file]
               _run_and_log(cmd6c, "step06_plot_RI_distribution")

            # SE Distribution Plot
            if "SE" in args.event_types:
               print("[INFO] Plotting SE Position Distribution with Genomic Background...")
               se_plot_file = os.path.join(step6_out, "SE_Position_Distribution.png")
               
               # Generate Background
               bg_file_se = os.path.join(run_outdir, "genome_SE_background_counts.tsv")
               if not os.path.exists(bg_file_se): 
                   print("[INFO] Generating Genome-wide SE Background (Exons)...")
                   cmd_bg_se = ["mamba", "run", "-n", "splicing-functional", "python3", os.path.join(script_dir, "src/generate_genomic_SE_background.py"), "--gtf", gtf_file, "--out_tsv", bg_file_se]
                   _run_and_log(cmd_bg_se, "step06_generate_se_background")
               
               print("[INFO] Plotting SE Odds Ratios and Saving Stats Table...")
               stats_file_se = os.path.join(step6_out, "SE_Position_Stats.tsv")
               cmd6d = ["mamba", "run", "-n", "splicing-functional", "python3", os.path.join(script_dir, "src/plot_SE_distribution.py"), "--per_event_file", per_event_file, "--output_file", se_plot_file, "--background_counts", bg_file_se, "--stats_table", stats_file_se]
               _run_and_log(cmd6d, "step06_plot_SE_distribution")

    # Step 7: Add Frame Shift Density (Formerly Step 9)
    if should_run(7):
        print("[INFO] === Step 7: Frame Shift Density Analysis ===")
        # Double check directory exists
        os.makedirs(step7_out, exist_ok=True)
        
        input_step7 = os.path.join(step6_out, "per_event_compact_for_plotting.tsv")
        if os.path.exists(input_step7):
            cmd_base = [
                "mamba", "run", "-n", "splicing-functional", "python3", os.path.join(script_dir, "src/add_frame_shify_density.py"),
                "--per_event", input_step7, "--gtf", gtf_file,
                "--normalize", args.normalize, "--nperm", args.nperm, "--nbins", args.nbins,
                "--bin_stats", "--annotate_binstats", "--overlay_datasets", "--one_plot_per_class",
                "--bin_fisher", "--write_bin_fisher_tsv", "--event_types"] + args.event_types
            if args.pool: cmd_base.append("--pool")

            if args.direction:
                # Run Inclusion Analysis -> dPSI_positive
                print("[INFO] Step 7: Running dPSI_positive (Inclusion) analysis...")
                cmd_inc = cmd_base + ["--out_prefix", f"{step7_out}/fig4_dPSI_positive", "--dataset_order", "UFM1_dependent_dPSI_positive,UFM1_independent_dPSI_positive"]
                _run_and_log(cmd_inc, "step07_frameshift_density_dPSI_positive")
                
                # Run Exclusion Analysis -> dPSI_negative
                print("[INFO] Step 7: Running dPSI_negative (Exclusion) analysis...")
                cmd_exc = cmd_base + ["--out_prefix", f"{step7_out}/fig4_dPSI_negative", "--dataset_order", "UFM1_dependent_dPSI_negative,UFM1_independent_dPSI_negative"]
                _run_and_log(cmd_exc, "step07_frameshift_density_dPSI_negative")
            else:
                # Standard Analysis
                cmd7 = cmd_base + ["--out_prefix", f"{step7_out}/fig4", "--dataset_order", "UFM1_dependent,UFM1_independent"] 
                _run_and_log(cmd7, "step07_frameshift_density")
        else:
            print(f"[WARN] Input file for Step 7 not found: {input_step7}. Skipping Step 7.")

    # Step 8: AA Features from Lost and Preserved (Formerly Step 5)
    if should_run(8):
        if dna_fasta_file:
            print("[INFO] === Step 8: Extracting AA Features ===")
            aa_features_cmd = ["mamba", "run", "-n", "splicing-functional", "bash", 
                               os.path.join(script_dir, "src/AAfeatures.sh"), 
                               "--lost", os.path.join(step1_out, "UFM1_dependent.tsv"), 
                               "--preserved", os.path.join(step1_out, "UFM1_independent.tsv"), 
                               "-g", dna_fasta_file, 
                               "-o", step8_out,
                               "--bg_gtf", gtf_file,
                               "--bg_cache_dir", os.path.join(args.cache_dir, "bg_exons")]
            if args.direction:
                aa_features_cmd.append("--direction")
            _run_and_log(aa_features_cmd, "step08_AA_features")
        else:
            print("[INFO] Skipping Step 8: AA Features (No DNA FASTA found)")

    # Step 9: Compute Protein Attributes (Biophysical Properties) (Formerly Step 6)
    if should_run(9):
        # Only run if Step 8 was likely possible or if input exists
        # Actually Step 9 takes input from Step 8 output.
        print("[INFO] === Step 9: Computing Protein Attributes (Biophysical Properties) ===")
        biophys_cmd = ["mamba", "run", "-n", "splicing-functional", "python3", 
                       os.path.join(script_dir, "src/biophysical_properties.py"),
                       "--indir", step8_out,
                       "--out_tsv", os.path.join(step9_out, "protein_attributes_properties.tsv")]
        _run_and_log(biophys_cmd, "step09_protein_attributes_calc")

    # Step 10: RI Motif & Deep Dive Analysis (New Feature)
    if should_run(10):
        if dna_fasta_file:
            print("[INFO] === Step 10: Motif & Deep Dive Analysis ===")
            step10_out = os.path.join(run_outdir, "step10_motif_analysis")
            os.makedirs(step10_out, exist_ok=True)
            
            # 10A/B: RI Analysis
            # 10A/B: RI Analysis
            if "RI" in args.event_types:
                print("[INFO] === Step 10: RI Motif & Deep Dive Analysis (Directional) ===")
                
                # Use same combined/pos/neg logic as SE
                analyses = [
                     ("combined", "UFM1_dependent.tsv", "UFM1_independent.tsv"),
                     ("dPSI_positive", "UFM1_dependent_dPSI_positive.tsv", "UFM1_independent_dPSI_positive.tsv"),
                     ("dPSI_negative", "UFM1_dependent_dPSI_negative.tsv", "UFM1_independent_dPSI_negative.tsv")
                ]
                
                for name, l_name, p_name in analyses:
                     l_path = os.path.join(step1_out, l_name)
                     p_path = os.path.join(step1_out, p_name)
                     
                     if not os.path.exists(l_path) or not os.path.exists(p_path):
                         print(f"[WARN] Skipping RI sub-analysis '{name}': inputs missing ({l_name}, {p_name}).")
                         continue

                     print(f"[INFO] Running RI Deep Dive: {name}...")
                     step10_ri_out = os.path.join(step10_out, f"RI_DeepDive_{name}")
                     os.makedirs(step10_ri_out, exist_ok=True)
                     
                     # 10A: Extract Sequences
                     cmd_10a = ["mamba", "run", "-n", "splicing-functional", "python3",
                                os.path.join(script_dir, "src/extract_ri_motifs.py"),
                                "--lost", l_path,
                                "--preserved", p_path,
                                "--genome_fasta", dna_fasta_file,
                                "--outdir", step10_ri_out]
                     _run_and_log(cmd_10a, f"step10_a_extract_ri_{name}")
                     
                     # 10B: Deep Dive
                     lost_fa = os.path.join(step10_ri_out, "lost.intron.fa")
                     pres_fa = os.path.join(step10_ri_out, "preserved.intron.fa")
                     
                     if os.path.exists(lost_fa) and os.path.exists(pres_fa) and gtf_file:
                          cmd_10b = ["mamba", "run", "-n", "splicing-functional", "python3",
                                    os.path.join(script_dir, "src/analyze_ri_vs_constitutive.py"),
                                    "--gtf", gtf_file,
                                    "--genome_fasta", dna_fasta_file,
                                    "--lost_tsv", l_path,
                                    "--preserved_tsv", p_path,
                                    "--lost_fa", lost_fa,
                                    "--preserved_fa", pres_fa,
                                    "--outdir", step10_ri_out,
                                    "--script_dir", os.path.join(script_dir, "src")]
                          
                          if args.motif_db:
                              cmd_10b.append(f"--motif_db={args.motif_db}")
                              
                          _run_and_log(cmd_10b, f"step10_b_deep_dive_ri_{name}")
                     else:
                          print(f"[WARN] Skipping Step 10B RI {name} (Missing inputs).")

            # 10C: SE Analysis
            if "SE" in args.event_types:
                 print("[INFO] === Step 10C: SE Motif & Deep Dive Analysis (Directional) ===")
                 
                 # Define sub-analyses: 'combined', 'pos' (dPSI>0), 'neg' (dPSI<0)
                 analyses = [
                     ("combined", "UFM1_dependent.tsv", "UFM1_independent.tsv"),
                     ("dPSI_positive", "UFM1_dependent_dPSI_positive.tsv", "UFM1_independent_dPSI_positive.tsv"),
                     ("dPSI_negative", "UFM1_dependent_dPSI_negative.tsv", "UFM1_independent_dPSI_negative.tsv")
                 ]
                 
                 for name, l_name, p_name in analyses:
                     l_path = os.path.join(step1_out, l_name)
                     p_path = os.path.join(step1_out, p_name)
                     
                     if not os.path.exists(l_path) or not os.path.exists(p_path):
                         print(f"[WARN] Skipping SE sub-analysis '{name}': inputs missing ({l_name}, {p_name}). (Maybe run Step 1 again?)")
                         continue
                     
                     print(f"[INFO] Running SE Deep Dive: {name}...")
                     step10_se_out = os.path.join(step10_out, f"SE_DeepDive_{name}")
                     os.makedirs(step10_se_out, exist_ok=True)
                     
                     cmd_10c = ["mamba", "run", "-n", "splicing-functional", "python3", 
                                os.path.join(script_dir, "src/analyze_se_vs_constitutive.py"),
                                "--lost", l_path,
                                "--preserved", p_path,
                                "--genome_fasta", dna_fasta_file,
                                "--outdir", step10_se_out]
                     
                     if args.motif_db:
                         cmd_10c.append(f"--motif_db={args.motif_db}")
                         
                     _run_and_log(cmd_10c, f"step10_c_se_{name}")
             
            if "RI" not in args.event_types and "SE" not in args.event_types:
                 print("[INFO] Skipping Step 10 actions (Requires 'RI' or 'SE' in --event_types).")
        else:
             print("[INFO] Skipping Step 10 (Requires DNA FASTA).")
             
    # Step 11: Mechanism Investigation (Stalling & Adjacency)
    if should_run(11):
        print("[INFO] === Step 11: Mechanism Investigation (Stalling & Adjacency) ===")
        step11_out = os.path.join(run_outdir, "step11_mechanism_investigation")
        os.makedirs(step11_out, exist_ok=True)
        step10_out = os.path.join(run_outdir, "step10_motif_analysis")
        
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
             if "SE" in args.event_types: target_dirs.append(("SE", os.path.join(step10_out, f"SE_DeepDive_{direction}")))
             if "RI" in args.event_types: target_dirs.append(("RI", os.path.join(step10_out, f"RI_DeepDive_{direction}")))
             
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
                     input_cands = [("UFM1_dependent", os.path.join(deep_dir, "lost.exon.fa")), ("UFM1_independent", os.path.join(deep_dir, "preserved.exon.fa"))]
                 else:
                     input_cands = [("UFM1_dependent", os.path.join(deep_dir, "lost.intron.fa")), ("UFM1_independent", os.path.join(deep_dir, "preserved.intron.fa"))]
                 
                 # Verify files exist
                 valid_groups = []
                 for name, fpath in input_cands:
                     if os.path.exists(fpath): valid_groups.append(f"{name}:{fpath}")
                 
                 if valid_groups:
                     out_subdir = os.path.join(step11_out, f"{etype}_{direction}")
                     os.makedirs(out_subdir, exist_ok=True)
                     
                     # 1. Stalling
                     cmd_11a = ["mamba", "run", "-n", "splicing-functional", "python3", 
                                os.path.join(script_dir, "src/analyze_stalling_motifs.py"),
                                "--groups"] + valid_groups + ["--outdir", out_subdir]
                     _run_and_log(cmd_11a, f"step11_a_stalling_{etype}_{direction}")
                     
                     # 2. GC/Len Check
                     # check_gc_len needs explicit lost/pres inputs
                     if len(input_cands) >= 2 and os.path.exists(input_cands[0][1]) and os.path.exists(input_cands[1][1]):
                          print(f"[INFO] Step 11C ({etype}_{direction}): Running Sequence Property Check...")
                          cmd_11c = ["mamba", "run", "-n", "splicing-functional", "python3",
                                     os.path.join(script_dir, "src/check_gc_len.py"),
                                     "--lost", input_cands[0][1],
                                     "--preserved", input_cands[1][1]]
                          _run_and_log(cmd_11c, f"step11_c_seq_properties_{etype}_{direction}")
                 else:
                     print(f"[WARN] inputs missing for {etype}_{direction} in {deep_dir}")

        # 11B: Adjacency (SE-RI) - Unchanged as it depends on Step 1 tables which are already split if needed
        # But analyze_se_ri_adjacency.py processes the full table usually.
        # User output showed Step 11B ran.
        l_path = os.path.join(step1_out, "UFM1_dependent.tsv")
        p_path = os.path.join(step1_out, "UFM1_independent.tsv")
        if os.path.exists(l_path):
             print("[INFO] Step 11B: Running SE-RI Adjacency Analysis...")
             adj_out = os.path.join(step11_out, "adjacency")
             os.makedirs(adj_out, exist_ok=True)
             cmd_11b = ["mamba", "run", "-n", "splicing-functional", "python3", 
                        os.path.join(script_dir, "src/analyze_se_ri_adjacency.py"),
                        "--se_files", l_path, p_path,
                        "--ri_files", l_path, p_path,
                        "--outdir", adj_out]
             _run_and_log(cmd_11b, "step11_b_adjacency")
        
    print("\n" + "="*60)
    print("       PIPELINE COMPLETED SUCCESSFULLY")
    print("="*60 + "\n")
    print(f"Results are available in: {run_outdir}")
        
if __name__ == "__main__":
    main()
