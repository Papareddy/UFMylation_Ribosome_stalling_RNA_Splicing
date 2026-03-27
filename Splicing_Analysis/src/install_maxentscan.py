#!/usr/bin/env python3
import os
import urllib.request
import tarfile
import shutil
import stat

def install_maxentscan(dest_dir="src/maxentscan"):
    if os.path.exists(dest_dir):
        print(f"[INFO] MaxEntScan directory exists: {dest_dir}. Skipping download.")
        return

    url = "http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz"
    tar_path = "fordownload.tar.gz"
    
    print(f"[INFO] Downloading MaxEntScan from {url}...")
    try:
        # User might not have internet access in restricted env, but we'll try.
        # If this fails, we might need a backup (e.g. mocking it).
        urllib.request.urlretrieve(url, tar_path)
    except Exception as e:
        print(f"[ERROR] Download failed: {e}")
        return

    print(f"[INFO] Extracting to {dest_dir}...")
    try:
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(path="temp_maxent")
            
        # The tar usually contains a 'fordownload' directory
        src_inner = os.path.join("temp_maxent", "fordownload")
        if os.path.exists(src_inner):
            shutil.move(src_inner, dest_dir)
        else:
            # Maybe flat?
            shutil.move("temp_maxent", dest_dir)
            
        # Cleanup
        os.remove(tar_path)
        if os.path.exists("temp_maxent"): shutil.rmtree("temp_maxent")
        
        # Verify and Fix Permissions
        for f in ["score5.pl", "score3.pl"]:
            p = os.path.join(dest_dir, f)
            if os.path.exists(p):
                st = os.stat(p)
                os.chmod(p, st.st_mode | stat.S_IEXEC)
                
        print("[SUCCESS] MaxEntScan installed.")
        
    except Exception as e:
        print(f"[ERROR] Installation failed: {e}")

if __name__ == "__main__":
    install_maxentscan()
