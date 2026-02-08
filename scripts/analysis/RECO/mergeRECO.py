import os
import subprocess

# Paths for MC
base_path = "/afs/cern.ch/user/s/scattola/FCCAnalyses/Ztautau/treemaker/"
data_dir = os.path.join(base_path, "p8_ee_Ztautau_ecm91")

def merge_and_clean(target_dir):
    output_file = f"{target_dir}.root"
    # hadd -f overwrites the target if it exists
    cmd = ["hadd", "-f", output_file, f"{target_dir}/chunk_*.root"]
    
    print(f">>> Merging chunks in {target_dir} into {output_file}...")
    # Use shell=True to allow the wildcard '*' to be expanded by the shell
    result = subprocess.run(" ".join(cmd), shell=True)
    
    if result.returncode == 0:
        print(f"Success! Deleting chunks in {target_dir}...")
        subprocess.run(f"rm {target_dir}/chunk_*.root", shell=True)
    else:
        print(f"Error merging {target_dir}. Chunks were NOT deleted.")

# Run for both data and templates
merge_and_clean(data_dir)
