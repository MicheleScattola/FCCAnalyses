import os
import subprocess
import glob

# Configuration
base_path = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/segmented_analysis/"
output_file = os.path.join(base_path, "templates_merged.root")

# Find all templates_*.root files
template_files = sorted(glob.glob(os.path.join(base_path, "templates_*.root")))

if not template_files:
    print(f"No templates_*.root files found in {base_path}")
    exit(1)

print(f"Found {len(template_files)} template files:")
for f in template_files:
    print(f"  - {os.path.basename(f)}")

# Merge all template files into one
cmd = ["hadd", "-f", output_file] + template_files

print(f"\n>>> Merging {len(template_files)} template files into {output_file}...")
result = subprocess.run(cmd)

if result.returncode == 0:
    print(f"\nSuccess! Merged file created: {output_file}")
    print(f"\nDeleting individual template files...")
    for f in template_files:
        try:
            os.remove(f)
            print(f"  - Deleted {os.path.basename(f)}")
        except OSError as e:
            print(f"  - Error deleting {os.path.basename(f)}: {e}")
else:
    print(f"\nError merging template files.")
