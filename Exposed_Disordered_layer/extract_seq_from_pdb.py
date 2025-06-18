import sys
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

pdb_file = sys.argv[1]
out_fasta = sys.argv[2]

parser = PDBParser(QUIET=True)
structure = parser.get_structure("model", pdb_file)

# Extract sequence from first chain with residues
for model in structure:
    for chain in model:
        residues = [res for res in chain if res.id[0] == ' ']
        seq = ''.join(res.resname for res in residues)
        # Convert 3-letter to 1-letter
        from Bio.Data import SCOPData
        aa_seq = ''.join(SCOPData.protein_letters_3to1.get(res.resname.capitalize(), 'X') for res in residues)

        record = SeqRecord(Seq(aa_seq), id=chain.id, description="")
        SeqIO.write(record, out_fasta, "fasta")
        break  # Only first chain
    break
