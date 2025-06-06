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

        ## compute invariant mass of two leading jets
        df = df.Define("jet_p4", "JetConstituentsUtils::compute_tlv_jets({})".format(jetClusteringHelper.jets))
        df = df.Define("event_invariant_mass", "JetConstituentsUtils::InvariantMass(jet_p4[0], jet_p4[1])")

        return df

    # __________________________________________________________
    # Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branches_pfcand = list(variables_pfcand.keys())
        branches_jet = list(variables_jet.keys())
        branches_event = list(variables_event.keys())

        branchList = branches_event + branches_jet + branches_pfcand

        return branchList
