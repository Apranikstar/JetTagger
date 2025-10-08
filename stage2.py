import os
import uproot
import torch
import math

def process_root_file(
    filepath: str,
    file_name: str,
    tree: str,
    output_dir: str,
    chunk_size: int = 100_000
):
    # --- Ensure output directory exists ---
    os.makedirs(output_dir, exist_ok=True)

    # List of per-constituent keys
    pfcand_keys = [
        'pfcand_isMu', 'pfcand_isEl', 'pfcand_isChargedHad', 'pfcand_isGamma',
        'pfcand_isNeutralHad', 'pfcand_e', 'pfcand_p', 'pfcand_theta', 'pfcand_phi',
        'pfcand_charge', 'pfcand_type', 'pfcand_erel', 'pfcand_erel_log',
        'pfcand_thetarel', 'pfcand_phirel', 'pfcand_dxy', 'pfcand_dz', 'pfcand_phi0',
        'pfcand_C', 'pfcand_ct', 'pfcand_dptdpt', 'pfcand_dxydxy', 'pfcand_dzdz',
        'pfcand_dphidphi', 'pfcand_detadeta', 'pfcand_dxydz', 'pfcand_dphidxy',
        'pfcand_phidz', 'pfcand_phictgtheta', 'pfcand_dxyctgtheta', 'pfcand_dlambdadz',
        'pfcand_cctgtheta', 'pfcand_phic', 'pfcand_dxyc', 'pfcand_cdz',
        'pfcand_btagSip2dVal', 'pfcand_btagSip2dSig', 'pfcand_btagSip3dVal',
        'pfcand_btagSip3dSig', 'pfcand_btagJetDistVal', 'pfcand_btagJetDistSig'
    ]

    # List of jet-level branches to include directly
    jet_branches = [
        'jet_mass', 'jet_p', 'jet_e', 'jet_phi', 'jet_theta', 'jet_pT',
        'jet_nnhad', 'jet_ngamma', 'jet_nchad', 'jet_nel', 'jet_nmu', 'jet_nconst'
    ]

    # Open ROOT file
    file_path = os.path.join(filepath, file_name)
    file = uproot.open(file_path)[tree]
    num_events = file.num_entries

    leading_jet_index = 0
    subleading_jet_index = 1

    num_chunks = math.ceil(num_events / chunk_size)

    for chunk_idx in range(num_chunks):
        start = chunk_idx * chunk_size
        stop = min((chunk_idx + 1) * chunk_size, num_events)
        print(f"Processing {file_name}: events {start} → {stop-1}")

        save_dict = {}

        # --- Process per-constituent keys ---
        for key in pfcand_keys:
            arr = file[key].array(library="np", entry_start=start, entry_stop=stop)
            jets_list = []
            for i in range(len(arr)):  # only within this chunk
                jets_list.append(arr[i][leading_jet_index])
                jets_list.append(arr[i][subleading_jet_index])
            max_len = max(len(jet) for jet in jets_list) if jets_list else 0
            num_jets = len(jets_list)
            jets_tensor = torch.zeros((num_jets, max_len), dtype=torch.float64)
            for i, jet in enumerate(jets_list):
                jets_tensor[i, :len(jet)] = torch.tensor(jet, dtype=torch.float64)
            save_dict[key] = jets_tensor

        # --- Initialize per-jet labels ---
        num_jets_total = (stop - start) * 2
        for label in [
            "recojet_isG","recojet_isB","recojet_isTHAD","recojet_isTLEP",
            "recojet_isWHAD","recojet_isZHAD","recojet_isTAU","recojet_isUUDDSS"
        ]:
            save_dict[label] = torch.zeros(num_jets_total, dtype=torch.int)

        # --- Assign label based on file name ---
        if "_gg_" in file_name:
            save_dict["recojet_isG"][:] = 1
        elif "_bb_" in file_name:
            save_dict["recojet_isB"][:] = 1
        elif "_thadthad_" in file_name:
            save_dict["recojet_isTHAD"][:] = 1
        elif "_tleptlep_" in file_name:
            save_dict["recojet_isTLEP"][:] = 1
        elif "_whadwhad_" in file_name:
            save_dict["recojet_isWHAD"][:] = 1
        elif "_zhadzhad_" in file_name:
            save_dict["recojet_isZHAD"][:] = 1
        elif "_tautau_" in file_name:
            save_dict["recojet_isTAU"][:] = 1
        elif "_uudds_" in file_name:
            save_dict["recojet_isUUDDSS"][:] = 1

        # --- Add jet-level branches directly ---
        for branch in jet_branches:
            arr = file[branch].array(library="np", entry_start=start, entry_stop=stop)
            branch_list = []
            for i in range(len(arr)):
                branch_list.append(arr[i][leading_jet_index])
                branch_list.append(arr[i][subleading_jet_index])
            save_dict[branch] = torch.tensor(branch_list, dtype=torch.float64)

        # --- Save chunk to output directory ---
        output_path = os.path.join(output_dir, f"{os.path.splitext(file_name)[0]}_chunk{chunk_idx}.pt")
        torch.save(save_dict, output_path)
        print(f"Saved {output_path}")

# --- Example usage ---
if __name__ == "__main__":
    tree = "events"
    input_dir = "/eos/user/h/hfatehi/JetTagger/TESTPIPELINE/"
    output_dir = "/eos/user/h/hfatehi/JetTagger/TESTPIPELINE/processedJets"

    for fname in os.listdir(input_dir):
        if fname.endswith(".root"):
            process_root_file(
                filepath=input_dir,
                file_name=fname,
                tree=tree,
                output_dir=output_dir,
                chunk_size=100_000
            )
