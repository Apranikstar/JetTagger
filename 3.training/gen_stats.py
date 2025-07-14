import glob
import json
import awkward as ak
import uproot
import numpy as np

# List of variables from your YAML config
variables = [
    "pfcand_erel_log", "pfcand_thetarel", "pfcand_phirel", "pfcand_dptdpt", "pfcand_detadeta",
    "pfcand_dphidphi", "pfcand_dxydxy", "pfcand_dzdz", "pfcand_dxydz", "pfcand_dphidxy",
    "pfcand_dlambdadz", "pfcand_dxyc", "pfcand_dxyctgtheta", "pfcand_phic", "pfcand_phidz",
    "pfcand_phictgtheta", "pfcand_cdz", "pfcand_cctgtheta", "pfcand_charge", "pfcand_isMu",
    "pfcand_isEl", "pfcand_isChargedHad", "pfcand_isGamma", "pfcand_isNeutralHad", "pfcand_type",
    "pfcand_dxy", "pfcand_dz", "pfcand_btagSip2dVal", "pfcand_btagSip2dSig", "pfcand_btagSip3dVal",
    "pfcand_btagSip3dSig", "pfcand_btagJetDistVal", "pfcand_btagJetDistSig", "pfcand_theta", "pfcand_phi"
]

# Your ROOT files
root_files = glob.glob("/afs/cern.ch/work/m/mmalekho/JetTagger/3.training/JetTagger-pTMin50/stage2_*.root")

# Open the files
arrays = uproot.concatenate(root_files, expressions=variables)

# Compute stats
stats = {}
for var in variables:
    flat_array = ak.flatten(arrays[var], axis=None)
    flat_array = flat_array[~ak.is_none(flat_array) & ~np.isnan(flat_array)]
    mean = float(np.mean(flat_array))
    std = float(np.std(flat_array))
    stats[var] = {"mean": mean, "std": std}

# Save to JSON
with open("data_stats.json", "w") as f:
    json.dump(stats, f, indent=2)

print("✅ data_stats.json generated.")
