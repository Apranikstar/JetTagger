import sys
import os
import argparse
import subprocess
import shutil
from datetime import date
import time
import concurrent.futures
import ROOT

# ________________________________________________________________________________
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--indir",
        help="Path to input directory",
        default="/eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/",
    )
    parser.add_argument(
        "--outdir",
        help="Path to output directory",
        default=os.path.join(os.getcwd(), "output"),
    )
    parser.add_argument(
        "--sample",
        help="Sample name template (use {process} as placeholder)",
        default="mgp8_pp_{process}_HT_2000_100000_5f_84TeV",
    )
    parser.add_argument("--ncpus", help="Number of CPUs", type=int, default=64)
    parser.add_argument(
        "--opt",
        help="Option: 1=run stage 1, 2=run stage 2, 3=all, 4=clean",
        default="3",
    )

    args = parser.parse_args()
    indir = args.indir
    outdir = args.outdir
    ncpus = args.ncpus
    sample = args.sample
    opt = args.opt

    # Hardcoded flavors list
    flavors = ["bb", "cc", "ss", "gg", "qq", "tautau", "tt"]

    flavor_to_process = {
    "qq": "jj",  # THIS is the fix!
    "bb": "bb",
    "cc": "cc",
    "ss": "ss",
    "gg": "gg",
    "tautau": "tautau",
    "tt": "tt",
}

    outtmpdir = os.path.join(os.getcwd(), "tmp")

    # Clean and create temporary directories
    if os.path.exists(outtmpdir):
        shutil.rmtree(outtmpdir)
    os.makedirs(outtmpdir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)

    # Prepare stage1 filenames
    #stage1_files = {f: f"{outtmpdir}/stage1_H{f}.root" for f in flavors}
    stage1_files = {
    f: f"{outtmpdir}/stage1_{flavor_to_process[f]}.root" for f in flavors
    }

    for f in flavors:
        #sample_f = sample.format(process=f)
        sample_f = sample.format(process=flavor_to_process[f])
        edm_files = f"{indir}/{sample_f}/*.root"

        # Run stage 1
        if opt in ["1", "3"]:
            cmd_stage1 = [
                "fccanalysis", "run", "stage1.py",
                "--output", stage1_files[f],
                "--files-list", edm_files,
                "--ncpus", str(ncpus)
            ]
            print("\nRunning stage 1:\n", " ".join(cmd_stage1), "\n")
            subprocess.run(" ".join(cmd_stage1), shell=True, check=True)

        # Run stage 2
        if opt in ["2", "3"]:
            nevents = count_events(stage1_files[f])
            nevents_per_thread = int(nevents / ncpus)

            stage2_files = {}
            commands_stage2 = []
            stage2_final_file = f"{outtmpdir}/stage2_{f}.root"
            stage2_wild_files = f"{outtmpdir}/stage2_{f}_*.root"
            hadd_cmd = f"hadd -f {stage2_final_file} {stage2_wild_files}"

            for i in range(ncpus):
                stage2_file = f"{outtmpdir}/stage2_{f}_{i}.root"
                stage2_files[i] = stage2_file
                nstart = i * nevents_per_thread
                nend = nstart + nevents_per_thread

                cmd_stage2 = [
                    "python", "stage2.py",
                    stage1_files[f],
                    stage2_file,
                    str(nstart),
                    str(nend),
                ]
                commands_stage2.append(cmd_stage2)

            # Run stage 2 in parallel using ProcessPoolExecutor for better CPU usage
            with concurrent.futures.ProcessPoolExecutor(max_workers=ncpus) as executor:
                futures = [executor.submit(run_command, cmd) for cmd in commands_stage2]
                concurrent.futures.wait(futures)

            # Merge stage 2 files
            print(f"\nMerging stage 2 files into: {stage2_final_file}\n{hadd_cmd}\n")
            subprocess.run(hadd_cmd, shell=True, check=True)

            # Copy to output directory
            shutil.copy(stage2_final_file, outdir)
            print(f"File copied to: {outdir}")

            # Cleanup
            print("Cleaning up temporary files...")
            for file_path in [stage2_final_file, stage1_files[f]]:
                if os.path.exists(file_path):
                    os.remove(file_path)
            subprocess.run(f"rm -f {stage2_wild_files}", shell=True)
            print("Done.")

    # Option 4: clean
    if opt == "4":
        print("Cleaning temporary directory...")
        shutil.rmtree(outtmpdir, ignore_errors=True)
        print("Done.")

# ________________________________________________________________________________
def run_command(command):
    print(f"Running command: {' '.join(command)}")
    subprocess.run(command, check=True)

# ________________________________________________________________________________
def count_events(file, tree_name="events"):
    root_file = ROOT.TFile.Open(file)
    tree = root_file.Get(tree_name)
    return tree.GetEntries()

# ________________________________________________________________________________
if __name__ == "__main__":
    main()
