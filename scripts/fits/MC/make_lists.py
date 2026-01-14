import glob
import os

def write_list(filename, file_list):
    with open(filename, 'w') as f:
        for path in file_list:
            # write absolute path just to be 100% safe
            f.write(f"{os.path.abspath(path)}\n") 
    print(f"Created {filename} with {len(file_list)} files.")

# --- Configuration ---
# NOTE: Ensure this path ends with a slash or os.path.join logic handles it
input_dir = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Ztautau_ecm91/"
pattern   = "events_*.root"

# --- Main Logic ---
full_search_path = os.path.join(input_dir, pattern)
all_files = sorted(glob.glob(full_search_path))

if len(all_files) == 0:
    print(f"ERROR: No files found in {input_dir}")
    exit(1)

# Split: 10 for templates, 100 for data
templates_files = all_files[0:10]
analysis_files  = all_files[10:40]

write_list("files_templates.txt", templates_files)
write_list("files_data.txt",      analysis_files)