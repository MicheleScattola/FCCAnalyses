import os

# Configuration
events_dir = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Ztautau_ecm91/"
output_dir = "/afs/cern.ch/user/s/scattola/FCCAnalyses/scripts/analysis/RECO/split_analysis/lists"
output_dir_root = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/split_analysis/"
os.makedirs(output_dir, exist_ok=True)

# Discover all ROOT files in the events directory
all_files = []
for entry in sorted(os.listdir(events_dir)):
    candidate = os.path.join(events_dir, entry)
    if entry.endswith(".root") and os.path.isfile(candidate):
        all_files.append(candidate)

# 1. Create Template Lists (First 50M events - 100k events per file -> 500 files)
n_template_files = 200
template_files = all_files[:n_template_files]

# Split templates into 20 jobs of 1M events each (10 files/job at 100k events/file)
template_files_per_job = 10
n_template_jobs = 20

job_manifest = []

for i in range(n_template_jobs):
    start_idx = i * template_files_per_job
    end_idx = start_idx + template_files_per_job
    chunk_files = template_files[start_idx:end_idx]

    list_name = f"{output_dir}/templates_{i}.txt"
    output_root = f"templates_{i}.root"

    with open(list_name, "w") as f:
        for line in chunk_files:
            f.write(line + "\n")

    abs_list_path = os.path.abspath(list_name)
    job_manifest.append(f"{abs_list_path} {output_root} templates.py")

# 2. Create Lists for Fits (Next 20M events each - 10 files per job)
# Start after the template files
start_idx = n_template_files
data_files_per_job = 10 # 1M events
n_data_jobs = 20

for i in range(n_data_jobs):
    end_idx = start_idx + data_files_per_job
    chunk_files = all_files[start_idx:end_idx]
    
    list_name = f"{output_dir}/data_{i}.txt"
    output_root = f"output_{i}.root"
    
    # Save the chunk list
    with open(list_name, "w") as f:
        for line in chunk_files:
            f.write(line + "\n")
            
    # Add to manifest: "path/to/input_list.txt  output_filename.root"
    # Use absolute paths for safety
    abs_list_path = os.path.abspath(list_name)
    job_manifest.append(f"{abs_list_path} {output_root} analysis.py")
    
    start_idx = end_idx

# Save the Manifest file for Condor
with open("job_manifest.txt", "w") as f:
    for line in job_manifest:
        f.write(line + "\n")

print(f"Created {n_template_jobs} template lists and {n_data_jobs} fit lists.")
print("Generated 'job_manifest.txt' for Condor.")
