from Bio.PDB import PDBParser, DSSP
import os
import sys

def parse_dssp(pdb_file):
    dssp_path = "/cs/labs/dina/noabirman/NES_hackathon/3d-biology-hackathon/dssp/mkdssp"

    # Parse structure
    structure = PDBParser(QUIET=True).get_structure("model", pdb_file)
    model = structure[0]

    # Run DSSP (requires mkdssp installed)
    #dssp = DSSP(model, pdb_file)
    dssp = DSSP(model, pdb_file, dssp=dssp_path)

    helix_list = []
    exposure_list = []

    for key in dssp.keys():
        aa_info = dssp[key]

        # Secondary structure code (H, B, E, G, I, T, S, or ' ')
        ss = aa_info[2]
        helix = ss in {"H", "G", "I"}  # Alpha, 3-10, Pi helix

        # Relative solvent accessibility (RSA)
        rsa = aa_info[3]
        exposed = rsa > 0.25  # Change threshold as needed

        helix_list.append(helix)
        exposure_list.append(exposed)

    return helix_list, exposure_list

if __name__ == "__main__":
    pdb_path = sys.argv[1]
    helix_list, exposure_list = parse_dssp(pdb_path)

    # Example output
    for i, (h, e) in enumerate(zip(helix_list, exposure_list)):
        print(f"Residue {i+1:4}: Helix={h} | Exposed={e}")
