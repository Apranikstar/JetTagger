import sys
import os
import argparse
import subprocess
import shutil
from datetime import date
import concurrent.futures
import ROOT


def get_process_from_sample(sample):
    import re
    match = re.search(r"mgp8_pp_(.+?)_(?:HT|Q)", sample)
    if match:
        return match.group(1)
    else:
        raise ValueError(f"Could not extract process from sample name: {sample}")


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
    parser.add_argument("--sample", help="sample name", default="mgp8_pp_tt_HT_2000_100000_5f_84TeV")

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

    flavor_to_process = {
        "jj": "jj",  # THIS is the fix!
        "qq": "qq",
        "bb": "bb",
        "cc": "cc",
        "ss": "ss",
        "gg": "gg",
        "tautau": "tautau",
        "tt": "tt",
        "uuddssgg": "qq",
    }
    process = get_process_from_sample(sample)
    if process not in flavor_to_process:
        raise ValueError(f"Process {process} not recognized in flavor_to_process mapping.")
    proc_name = flavor_to_process[process]

    outtmpdir = os.path.join(os.getcwd(), "tmp")

    # Clean and create temporary directories
    if os.path.exists(outtmpdir):
        shutil.rmtree(outtmpdir)
    os.makedirs(outtmpdir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)

    # Stage1
    stage1_files = f"{outtmpdir}/stage1_{proc_name}.root"
    sample_f = sample.replace("XX", process)
    edm_files = f"{indir}/{sample_f}/*.root"

    # -------------------------------------------------------------------------
    # Run stage 1
    if opt in ["1", "3"]:
        cmd_stage1 = [
            "fccanalysis", "run", "stage1.py",
            "--output", stage1_files,
            "--files-list", edm_files,
            "--ncpus", str(ncpus)
        ]
        print("\nRunning stage 1:\n", " ".join(cmd_stage1), "\n")
        subprocess.run(" ".join(cmd_stage1), shell=True, check=True)

    # -------------------------------------------------------------------------
    # Run stage 2
    if opt in ["2", "3"]:
        nevents = count_events(stage1_files)
        nworkers = min(ncpus, nevents)  # never more workers than events
        nevents_per_worker = nevents // nworkers
        remainder = nevents % nworkers

        stage2_files = {}
        commands_stage2 = []
        stage2_final_file = f"{outtmpdir}/stage2_{proc_name}.root"
        stage2_wild_files = f"{outtmpdir}/stage2_{proc_name}_*.root"
        hadd_cmd = f"hadd -f {stage2_final_file} {stage2_wild_files}"

        for i in range(nworkers):
            nstart = i * nevents_per_worker + min(i, remainder)
            nend = nstart + nevents_per_worker + (1 if i < remainder else 0)
            if i == nworkers - 1:
                nend = nevents

            stage2_file = f"{outtmpdir}/stage2_{proc_name}_{i}.root"
            stage2_files[i] = stage2_file
            cmd_stage2 = [
                "python", "stage2.py",
                stage1_files,
                stage2_file,
                str(nstart),
                str(nend),
            ]
            print(f"Worker {i}: processing events [{nstart}, {nend})")
            commands_stage2.append(cmd_stage2)

        # Run workers in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=nworkers) as executor:
            futures = {executor.submit(run_command, cmd): cmd for cmd in commands_stage2}
            for fut in concurrent.futures.as_completed(futures):
                cmd = futures[fut]
                try:
                    fut.result()
                except Exception as e:
                    print(f"Worker for command {' '.join(cmd)} failed with: {e}")
                    raise

        # Merge stage 2 files
        print(f"\nMerging stage 2 files into: {stage2_final_file}\n{hadd_cmd}\n")
        subprocess.run(hadd_cmd, shell=True, check=True)

        # Copy to output directory
        shutil.copy(stage2_final_file, outdir)
        print(f"File copied to: {outdir}")

        # Cleanup
        print("Cleaning up temporary files...")
        for file_path in [stage2_final_file, stage1_files]:
            if os.path.exists(file_path):
                os.remove(file_path)
        subprocess.run(f"rm -f {stage2_wild_files}", shell=True)
        print("Done.")

    # -------------------------------------------------------------------------
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
