from Bio.PDB import PDBParser, DSSP
import os
import sys

def parse_dssp(pdb_file):
    # Constants
    dssp_path = "/cs/labs/dina/noabirman/NES_hackathon/3d-biology-hackathon/dssp/mkdssp"
    pdb_file = "your_model.pdb"
    helix_codes = {"H", "G", "I"}  # Alpha, 3-10, pi helix
    exposed_threshold = 0.25  # RSA threshold

    # Parse structure
    parser = PDBParser()
    structure = parser.get_structure("model", pdb_file)
    model = structure[0]

    # Run DSSP
    dssp = DSSP(model, pdb_file, dssp=dssp_path)

    # Create boolean lists
    is_helix = []
    is_exposed = []

    # DSSP returns a dict keyed by (chain_id, res_id)
    for key in dssp:
        ss = dssp[key][2]  # Secondary structure code
        rsa = dssp[key][3]  # Relative Solvent Accessibility (RSA)

        is_helix.append(ss in helix_codes)
        is_exposed.append(rsa >= exposed_threshold)

    return is_helix, is_exposed

if __name__ == "__main__":
    pdb_path = sys.argv[1]
    is_helix, is_exposed = parse_dssp(pdb_path)

    print("Helix:", is_helix[:10])
    print("Exposed:", is_exposed[:10])
