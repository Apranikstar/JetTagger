import os, copy
import urllib.request

from config import (
    variables_pfcand,
    variables_jet,
    variables_event,
)

processList = {
    "mgp8_pp_tt_HT_2000_100000_5f_84TeV": {"fraction": 0.0001,},
    "mgp8_pp_jj_HT_2000_100000_5f_84TeV": {"fraction": 0.0001,},
}


includePaths = ["pxpypz.h"]

outputDir = "/eos/user/h/hfatehi/test"

inputDir = "/eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II"

nCPUS = -1

model_name = "model_Ep40_R39_1TeV_NM"

url_model_dir = "https://raw.githubusercontent.com/Apranikstar/JetTagger/main/3.training/inference"
url_preproc = f"{url_model_dir}/{model_name}.json"
url_model   = f"{url_model_dir}/{model_name}.onnx"

model_dir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
local_preproc = f"{model_dir}/{model_name}.json"
local_model   = f"{model_dir}/{model_name}.onnx"

def get_file_path(url, filename):
    if os.path.exists(filename):
        return os.path.abspath(filename)
    else:
        local_name = os.path.basename(filename)
        print(f"Downloading {url} → {local_name}")
        urllib.request.urlretrieve(url, local_name)
        return os.path.abspath(local_name)

weaver_preproc = get_file_path(url_preproc, local_preproc)
weaver_model   = get_file_path(url_model, local_model)

from jetFlavourHelper import JetFlavourHelper
from jetClusteringHelper import InclusiveJetClusteringHelper

jetFlavourHelper = None
jetClusteringHelper = None

# Mandatory: RDFanalysis class where the user defines the operations on the TTree
class RDFanalysis:
    # __________________________________________________________
    # Mandatory: analysers function to define the analysers to process,
    # please make sure you return the last dataframe
    def analysers(df):
        global jetClusteringHelper
        global jetFlavourHelper

        from config import collections, njets, ptcut, coneRadius

        ## define jet clustering parameters
        jetClusteringHelper = InclusiveJetClusteringHelper(
            collections["PFParticles"], coneRadius, ptcut
        )

        ## run jet clustering
        df = jetClusteringHelper.define(df)

        ## define jet flavour tagging parameters
        jetFlavourHelper = JetFlavourHelper(
            collections,
            jetClusteringHelper.jets,
            jetClusteringHelper.constituents,
            njets,
        )

        ## define observables for tagger
        df = jetFlavourHelper.define(df)
        df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df)

        ## compute invariant mass of two leading jets
        # df = df.Define(
        #     "jet_p4",
        #     "JetConstituentsUtils::compute_tlv_jets({})".format(jetClusteringHelper.jets),
        # )
        # df = df.Define(
        #     "event_invariant_mass",
        #     "JetConstituentsUtils::InvariantMass(jet_p4[0], jet_p4[1])",
        # )
        df = df.Filter("event_njet > 1")
        #df = df.Define("recoT", "scores_recojet_isT")
        #df = df.Filter("recoT > 0.90")

        return df

    # __________________________________________________________
    # Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = jetFlavourHelper.outputBranches()
        return branchList
