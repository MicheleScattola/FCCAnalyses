import os
import subprocess
import glob
import fnmatch
import time

# Configuration
base_path = "/eos/user/s/scattola/FCCAnalyses/Ztautau/treemaker/RECO/angular_analysis/"
output_file_templates = os.path.join(base_path, "templates_merged.root")
output_file_data = os.path.join(base_path, "output_merged.root")

# Build explicit lists for templates_i.root and output_i.root (i = 0..19)
template_files = [os.path.join(base_path, f"templates_{i}.root") for i in range(20)]
output_files = [os.path.join(base_path, f"output_{i}.root") for i in range(20)]

if not template_files and not output_files:
    print(f"No files found to merge in {base_path}")
    exit(1)

script_start = time.time()
print(f"=== Merge Script Started: {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")

# Merge template files directly to EOS
if template_files:
    print(f"Found {len(template_files)} template files:")
    for i, f in enumerate(template_files, 1):
        print(f"  [{i:2d}] {os.path.basename(f)}")

    merge_start = time.time()
    cmd = ["hadd", "-f", output_file_templates] + template_files

    print(f"\n>>> Starting template merge ({len(template_files)} files)...")
    print(f">>> Target: {output_file_templates}\n")
    result = subprocess.run(cmd)

    merge_elapsed = time.time() - merge_start

    if result.returncode == 0:
        print(f"\n✓ Template merge completed in {merge_elapsed:.2f}s")
        print(f"✓ Merged file: {output_file_templates}")
        
        print(f"\nDeleting {len(template_files)} individual template files...")
        delete_count = 0
        for i, f in enumerate(template_files, 1):
            try:
                os.remove(f)
                delete_count += 1
                print(f"  [{i:2d}/{len(template_files)}] ✓ Deleted {os.path.basename(f)}")
            except OSError as e:
                print(f"  [{i:2d}/{len(template_files)}] ✗ Error deleting {os.path.basename(f)}: {e}")
        print(f"✓ Deleted {delete_count}/{len(template_files)} template files")
    else:
        print(f"\n✗ Error merging template files (exit code: {result.returncode})")
        exit(1)

# Merge output files directly to EOS
if output_files:
    print(f"\nFound {len(output_files)} output files:")
    for i, f in enumerate(output_files, 1):
        print(f"  [{i:2d}] {os.path.basename(f)}")

    merge_start = time.time()
    cmd = ["hadd", "-f", output_file_data] + output_files

    print(f"\n>>> Starting output merge ({len(output_files)} files)...")
    print(f">>> Target: {output_file_data}\n")
    result = subprocess.run(cmd)

    merge_elapsed = time.time() - merge_start

    if result.returncode == 0:
        print(f"\n✓ Output merge completed in {merge_elapsed:.2f}s")
        print(f"✓ Merged file: {output_file_data}")
        
        print(f"\nDeleting {len(output_files)} individual output files...")
        delete_count = 0
        for i, f in enumerate(output_files, 1):
            try:
                os.remove(f)
                delete_count += 1
                print(f"  [{i:2d}/{len(output_files)}] ✓ Deleted {os.path.basename(f)}")
            except OSError as e:
                print(f"  [{i:2d}/{len(output_files)}] ✗ Error deleting {os.path.basename(f)}: {e}")
        print(f"✓ Deleted {delete_count}/{len(output_files)} output files")
    else:
        print(f"\n✗ Error merging output files (exit code: {result.returncode})")
        exit(1)

script_elapsed = time.time() - script_start
print(f"\n{'='*60}")
print(f"✓ All merges completed successfully!")
print(f"✓ Total script time: {script_elapsed:.2f}s ({int(script_elapsed//60)}m {int(script_elapsed%60)}s)")
print(f"✓ Script ended: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"{'='*60}")
