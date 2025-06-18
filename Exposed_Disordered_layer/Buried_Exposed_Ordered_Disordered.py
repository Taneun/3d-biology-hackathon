from Bio.PDB import PDBParser, DSSP
import os
import sys

def parse_dssp(pdb_path, dssp_exe="/cs/labs/dina/noabirman/NES_hackathon/3d-biology-hackathon/dssp/mkdssp", rsa_thresh=0.25):
    from Bio.PDB import PDBParser, DSSP

    helix_codes = {"H", "G", "I"}

    # Parse structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_path)
    model = structure[0]

    # Run DSSP
    dssp = DSSP(model, pdb_path, dssp=dssp_exe)

    # Initialize outputs
    is_helix = []
    is_exposed = []

    for key in dssp.keys():
        dssp_data = dssp[key]

        ss = dssp_data[1]      # secondary structure (1-letter code)
        #rsa = dssp_data[2]     # relative solvent accessibility (0-1)
        try:
            rsa = float(dssp_data[2])
        except ValueError:
            rsa = 0.0  # fallback: treat as buried if invalid

        is_helix.append(ss in helix_codes)
        is_exposed.append(rsa >= rsa_thresh)

    return is_helix, is_exposed



if __name__ == "__main__":
    pdb_path = sys.argv[1]
    is_helix, is_exposed = parse_dssp(pdb_path)

    print("Helix:", is_helix[:10])
    print("Exposed:", is_exposed[:10])
