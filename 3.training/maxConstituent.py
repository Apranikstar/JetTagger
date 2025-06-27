import uproot
import awkward as ak

def get_max_constituents_info(filename, treename="tree", branch="pfcand_e"):
    with uproot.open(filename) as file:
        tree = file[treename]
        pfcand_e = tree[branch].arrays(library="ak")[branch]
        n_constituents = ak.num(pfcand_e)
        max_index = ak.argmax(n_constituents)
        max_count = n_constituents[max_index]
        return max_index, max_count, pfcand_e[max_index]

# Analyze stage2_tt.root
tt_index, tt_count, tt_energies = get_max_constituents_info("stage2_tt.root")
print("== stage2_tt.root ==")
print(f"Event index with max constituents: {tt_index}")
print(f"Number of constituents: {tt_count}")
print(f"Energies: {tt_energies}\n")

# Analyze stage2_jj.root
jj_index, jj_count, jj_energies = get_max_constituents_info("stage2_jj.root")
print("== stage2_jj.root ==")
print(f"Event index with max constituents: {jj_index}")
print(f"Number of constituents: {jj_count}")
print(f"Energies: {jj_energies}")
