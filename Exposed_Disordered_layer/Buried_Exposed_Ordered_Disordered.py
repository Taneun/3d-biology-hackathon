import os
from Bio.PDB import PDBParser, DSSP
from Bio.PDB.Polypeptide import three_to_one

def parse_af_model(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_path)
    model = next(structure.get_models())
    dssp = DSSP(model, pdb_path)

    data = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != ' ':
                continue
            res_id = residue.id[1]
            try:
                aa = three_to_one(residue.resname)
            except:
                aa = 'X'
            try:
                dssp_key = (chain.id, (' ', res_id, ' '))
                rsa = dssp[dssp_key][3]  # relative solvent accessibility
            except Exception:
                rsa = None
            bfactor = residue['CA'].bfactor if 'CA' in residue else None
            data.append({
                'chain': chain.id,
                'res_id': res_id,
                'aa': aa,
                'rsa': rsa,
                'plddt': bfactor,
            })
    return data

def parse_iupred(iupred_path):
    scores = []
    with open(iupred_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    score = float(parts[2])
                    scores.append(score)
                except:
                    continue
    return scores

def label_residues(res_data, iupred_scores=None):
    for i, r in enumerate(res_data):
        r['buried'] = r['rsa'] is not None and r['rsa'] < 0.2
        r['disordered'] = False
        if r['plddt'] is not None and r['plddt'] < 70:
            r['disordered'] = True
        if iupred_scores and i < len(iupred_scores):
            r['iupred'] = iupred_scores[i]
            if r['iupred'] > 0.5:
                r['disordered'] = True
    return res_data

def filter_nes_candidates(res_data, window_size=12, threshold=0.5):
    candidates = []
    for i in range(len(res_data) - window_size + 1):
        window = res_data[i:i + window_size]
        buried_frac = sum(r['buried'] for r in window) / window_size
        disordered_frac = sum(r['disordered'] for r in window) / window_size
        if buried_frac < threshold and disordered_frac < threshold:
            candidates.append({
                'start_res_id': window[0]['res_id'],
                'end_res_id': window[-1]['res_id'],
                'sequence': ''.join(r['aa'] for r in window),
                'buried_fraction': buried_frac,
                'disordered_fraction': disordered_frac
            })
    return candidates

def main(pdb_path, iupred_path=None, window_size=12, threshold=0.5):
    print(f"Parsing AlphaFold model: {pdb_path}")
    res_data = parse_af_model(pdb_path)

    if iupred_path and os.path.exists(iupred_path):
        print(f"Loading IUPred3 disorder predictions: {iupred_path}")
        iupred_scores = parse_iupred(iupred_path)
    else:
        iupred_scores = None

    print("Labeling residues...")
    res_data = label_residues(res_data, iupred_scores)

    print(f"Scanning for NES candidate regions (window={window_size}, threshold={threshold})...")
    candidates = filter_nes_candidates(res_data, window_size, threshold)

    print(f"Found {len(candidates)} candidate regions:")
    for c in candidates:
        print(f" - {c['start_res_id']}–{c['end_res_id']}: {c['sequence']} "
              f"(buried={c['buried_fraction']:.2f}, disordered={c['disordered_fraction']:.2f})")

    return candidates

if __name__ == "__main__":
    # Replace with your paths
    pdb_path = "AF-P12345-F1-model_v4.pdb"
    iupred_path = "P12345.iupred"
    main(pdb_path, iupred_path)
