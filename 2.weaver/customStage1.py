from config import (
    variables_pfcand,
    variables_jet,
    variables_event,
)

from ONNXRuntime.python.jetFlavourHelper import JetFlavourHelper
from FastJet.python.jetClusteringHelper import InclusiveJetClusteringHelper

jetFlavourHelper = None
jetClusteringHelper = None

# Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis:
    # __________________________________________________________
    # Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        global jetClusteringHelper
        global jetFlavourHelper

        from config import collections, njets, ptcut, coneRadius

        ## define jet clustering parameters
        jetClusteringHelper = InclusiveJetClusteringHelper(collections["PFParticles"], coneRadius, ptcut)

        ## run jet clustering
        df = jetClusteringHelper.define(df)

        ## define jet flavour tagging parameters
        jetFlavourHelper = JetFlavourHelper(
            collections,
            jetClusteringHelper.jets,
            jetClusteringHelper.constituents,
            njets
        )

        ## define observables for tagger
        df = jetFlavourHelper.define(df)
        #df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df)

        ## compute invariant mass of two leading jets
        df = df.Define("jet_p4", "JetConstituentsUtils::compute_tlv_jets({})".format(jetClusteringHelper.jets))
        df = df.Define("event_invariant_mass", "JetConstituentsUtils::InvariantMass(jet_p4[0], jet_p4[1])")
                ## compute sum of four-momenta of constituents
        df = df.Define(
            "sumTLVs",
            "JetConstituentsUtils::sum_tlv_constituents({})".format(
                jetClusteringHelper.constituents
            ),
        )
        
        
       # df = df.Define("sumTLVs1M", "sumTLVs[0].M()")
       # df = df.Define("sumTLVs2M", "sumTLVs[1].M()")
       # df = df.Filter(" 125 < sumTLVs1M && sumTLVs1M < 225")
       # df = df.Filter(" 125 < sumTLVs2M && sumTLVs2M < 225")
        df = df.Redefine("jet_mass", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].M(), sumTLVs[1].M()})")

# and filter on the vector elements
        df = df.Filter("125 < jet_mass[0] && jet_mass[0] < 225")
        df = df.Filter("125 < jet_mass[1] && jet_mass[1] < 225")
        return df

    # __________________________________________________________
    # Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branches_pfcand = list(variables_pfcand.keys())
        branches_jet = list(variables_jet.keys())
        branches_event = list(variables_event.keys())

        branchList = branches_event + branches_jet + branches_pfcand #+ ["jet_mass"]

        return branchList

