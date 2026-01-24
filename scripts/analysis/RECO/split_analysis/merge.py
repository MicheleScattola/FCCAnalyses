import os
import subprocess
import glob

# Configuration
base_path = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/split_analysis/"
output_file_templates = os.path.join(base_path, "templates_merged.root")
output_file_data = os.path.join(base_path, "output_merged.root")

# Find all templates_*.root files
template_files = sorted(glob.glob(os.path.join(base_path, "templates_*.root")))

# Exclude templates_histograms.root and templates_merged.root
exclude_files = {"templates_histograms.root", "templates_merged.root"}
template_files = [f for f in template_files if os.path.basename(f) not in exclude_files]

# Find all output_*.root files
output_files = sorted(glob.glob(os.path.join(base_path, "output_*.root")))

# Exclude output_merged.root
output_files = [f for f in output_files if os.path.basename(f) != "output_merged.root"]

if not template_files and not output_files:
    print(f"No files found to merge in {base_path}")
    exit(1)

# Merge template files
if template_files:
    print(f"Found {len(template_files)} template files:")
    for f in template_files:
        print(f"  - {os.path.basename(f)}")

    cmd = ["hadd", "-f", output_file_templates] + template_files

    print(f"\n>>> Merging {len(template_files)} template files into {output_file_templates}...")
    result = subprocess.run(cmd)

    if result.returncode == 0:
        print(f"\nSuccess! Merged file created: {output_file_templates}")
        print(f"\nDeleting individual template files...")
        for f in template_files:
            try:
                os.remove(f)
                print(f"  - Deleted {os.path.basename(f)}")
            except OSError as e:
                print(f"  - Error deleting {os.path.basename(f)}: {e}")
    else:
        print(f"\nError merging template files.")

# Merge output files
if output_files:
    print(f"\nFound {len(output_files)} output files:")
    for f in output_files:
        print(f"  - {os.path.basename(f)}")

    cmd = ["hadd", "-f", output_file_data] + output_files

    print(f"\n>>> Merging {len(output_files)} output files into {output_file_data}...")
    result = subprocess.run(cmd)

    if result.returncode == 0:
        print(f"\nSuccess! Merged file created: {output_file_data}")
        print(f"\nDeleting individual output files...")
        for f in output_files:
            try:
                os.remove(f)
                print(f"  - Deleted {os.path.basename(f)}")
            except OSError as e:
                print(f"  - Error deleting {os.path.basename(f)}: {e}")
    else:
        print(f"\nError merging output files.")
