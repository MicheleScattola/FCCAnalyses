import os
import subprocess
import glob
import fnmatch
import time

# Configuration
eos_base_path = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/angular_analysis/"
xrootd_base = "root://eoshome-s.cern.ch//eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/angular_analysis/"
local_tmp_dir = "/tmp/scattola/"

# Create local temp directory
os.makedirs(local_tmp_dir, exist_ok=True)

# Local temp files
local_templates_merged = os.path.join(local_tmp_dir, "templates_merged.root")
local_output_merged = os.path.join(local_tmp_dir, "output_merged.root")

# Final EOS destinations
eos_templates_merged = os.path.join(eos_base_path, "templates_merged.root")
eos_output_merged = os.path.join(eos_base_path, "output_merged.root")

# Build explicit lists using xrootd protocol for reading from EOS
template_files = [xrootd_base + f"templates_{i}.root" for i in range(25)]
output_files = [xrootd_base + f"output_{i}.root" for i in range(25)]

if not template_files and not output_files:
    print(f"No files found to merge in {base_path}")
    exit(1)

script_start = time.time()
print(f"=== Merge Script Started: {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")

# Merge template files locally, then copy to EOS
if template_files:
    print(f"Found {len(template_files)} template files:")
    for i, f in enumerate(template_files, 1):
        print(f"  [{i:2d}] {os.path.basename(f)}")

    # Step 1: Merge to local disk
    merge_start = time.time()
    cmd = ["hadd", "-f", local_templates_merged] + template_files

    print(f"\n>>> Starting template merge ({len(template_files)} files)...")
    print(f">>> Reading from EOS, writing to local: {local_templates_merged}\n")
    result = subprocess.run(cmd)

    merge_elapsed = time.time() - merge_start

    if result.returncode != 0:
        print(f"\n✗ Error merging template files (exit code: {result.returncode})")
        exit(1)
    
    print(f"\n✓ Template merge completed in {merge_elapsed:.2f}s")
    print(f"✓ Local merged file: {local_templates_merged}")
    
    # Step 2: Copy to EOS
    print(f"\n>>> Copying merged file to EOS...")
    copy_start = time.time()
    copy_cmd = ["xrdcp", "-f", local_templates_merged, xrootd_base + "templates_merged.root"]
    result = subprocess.run(copy_cmd)
    copy_elapsed = time.time() - copy_start
    
    if result.returncode != 0:
        print(f"\n✗ Error copying to EOS (exit code: {result.returncode})")
        exit(1)
    
    print(f"✓ Copy to EOS completed in {copy_elapsed:.2f}s")
    print(f"✓ EOS file: {eos_templates_merged}")
    
    # Step 3: Delete individual files from EOS
    print(f"\nDeleting {len(template_files)} individual template files from EOS...")
    delete_count = 0
    for i in range(25):
        eos_file = os.path.join(eos_base_path, f"templates_{i}.root")
        try:
            os.remove(eos_file)
            delete_count += 1
            print(f"  [{i+1:2d}/{len(template_files)}] ✓ Deleted templates_{i}.root")
        except OSError as e:
            print(f"  [{i+1:2d}/{len(template_files)}] ✗ Error deleting templates_{i}.root: {e}")
    print(f"✓ Deleted {delete_count}/{len(template_files)} template files")
    
    # Step 4: Clean up local temp file
    try:
        os.remove(local_templates_merged)
        print(f"✓ Cleaned up local temp file")
    except OSError as e:
        print(f"⚠ Warning: Could not delete local temp file: {e}")

# Merge output files locally, then copy to EOS
if output_files:
    print(f"\nFound {len(output_files)} output files:")
    for i, f in enumerate(output_files, 1):
        print(f"  [{i:2d}] {os.path.basename(f)}")

    # Step 1: Merge to local disk
    merge_start = time.time()
    cmd = ["hadd", "-f", local_output_merged] + output_files

    print(f"\n>>> Starting output merge ({len(output_files)} files)...")
    print(f">>> Reading from EOS, writing to local: {local_output_merged}\n")
    result = subprocess.run(cmd)

    merge_elapsed = time.time() - merge_start

    if result.returncode != 0:
        print(f"\n✗ Error merging output files (exit code: {result.returncode})")
        exit(1)
    
    print(f"\n✓ Output merge completed in {merge_elapsed:.2f}s")
    print(f"✓ Local merged file: {local_output_merged}")
    
    # Step 2: Copy to EOS
    print(f"\n>>> Copying merged file to EOS...")
    copy_start = time.time()
    copy_cmd = ["xrdcp", "-f", local_output_merged, xrootd_base + "output_merged.root"]
    result = subprocess.run(copy_cmd)
    copy_elapsed = time.time() - copy_start
    
    if result.returncode != 0:
        print(f"\n✗ Error copying to EOS (exit code: {result.returncode})")
        exit(1)
    
    print(f"✓ Copy to EOS completed in {copy_elapsed:.2f}s")
    print(f"✓ EOS file: {eos_output_merged}")
    
    # Step 3: Delete individual files from EOS
    print(f"\nDeleting {len(output_files)} individual output files from EOS...")
    delete_count = 0
    for i in range(25):
        eos_file = os.path.join(eos_base_path, f"output_{i}.root")
        try:
            os.remove(eos_file)
            delete_count += 1
            print(f"  [{i+1:2d}/{len(output_files)}] ✓ Deleted output_{i}.root")
        except OSError as e:
            print(f"  [{i+1:2d}/{len(output_files)}] ✗ Error deleting output_{i}.root: {e}")
    print(f"✓ Deleted {delete_count}/{len(output_files)} output files")
    
    # Step 4: Clean up local temp file
    try:
        os.remove(local_output_merged)
        print(f"✓ Cleaned up local temp file")
    except OSError as e:
        print(f"⚠ Warning: Could not delete local temp file: {e}")

script_elapsed = time.time() - script_start
print(f"\n{'='*60}")
print(f"✓ All merges completed successfully!")
print(f"✓ Total script time: {script_elapsed:.2f}s ({int(script_elapsed//60)}m {int(script_elapsed%60)}s)")
print(f"✓ Script ended: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"{'='*60}")
