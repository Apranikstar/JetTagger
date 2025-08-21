processList = {

    # xsecs need to be scaled by 280/989 ...for xsec of ee -> H ...

    # Semileptonic processes
    "mgp8_pp_jj_HT_2000_100000_5f_84TeV": {"fraction": 1,},
    "mgp8_pp_tt_HT_2000_100000_5f_84TeV": {"fraction": 1,},
}

outputDir = "./output"
inputDir = "/eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/"
nCPUS = -1

from config import (
    variables_pfcand,
    variables_jet,
    variables_event,
)

from ONNXRuntime.python.jetFlavourHelper import JetFlavourHelper
from FastJet.python.jetClusteringHelper import InclusiveJetClusteringHelper

jetFlavourHelper = None
jetClusteringHelper = None


# Mandatory: RDFanalysis class where the user defines the operations on the TTree
class RDFanalysis:
    # __________________________________________________________
    # Mandatory: analysers function to define the analysers to process
    # Please make sure you return the last dataframe, in this example it is df
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
        # df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df)

        ## compute invariant mass of two leading jets
        df = df.Define(
            "jet_p4",
            "JetConstituentsUtils::compute_tlv_jets({})".format(jetClusteringHelper.jets),
        )
        df = df.Define(
            "event_invariant_mass",
            "JetConstituentsUtils::InvariantMass(jet_p4[0], jet_p4[1])",
        )

        ## compute sum of four-momenta of constituents
        df = df.Define(
            "sumTLVs1",
            "JetConstituentsUtils::sum_tlv_constituents({})".format(
                jetClusteringHelper.constituents
            ),
        )
        
        
        df = df.Define("sumTLVs1M", "sumTLVs1[0].M()")
        df = df.Define("sumTLVs2M", "sumTLVs1[1].M()")

        return df

    # __________________________________________________________
    # Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = ["sumTLVs1M", "sumTLVs2M", "event_invariant_mass"]
        return branchList


