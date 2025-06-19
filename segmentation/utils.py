import numpy as np
import ruptures as rpt
import argparse
import csv

def read_csv_sequences(csv_path):
    '''
    Reads in CSV file containing multiple sequences.
    CSV must have columns "Fasta Header" and "Sequence".
    Returns dictionary holding multiple sequences with UniProt IDs as keys.
    '''

    sequences = dict()

    with open(csv_path, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file)

        for row in csv_reader:
            if row['Sequence'] == '':
                continue
            uniprot_id = row['uniprotID']
            sequence = row['Sequence']

            uniprot_id = uniprot_id.strip()

            # Store the sequence in uppercase with gaps removed
            sequences[uniprot_id] = sequence.upper().replace("-", "")

    return sequences
    # return dict(list(sequences.items())[:2])

def get_protein_segments(emb_dict, max_bkps_per100aa=3):
    """
    Define the boundaries of protein segments using change point analysis.

    Parameters:
    emb_dict: a dictionary of UniProt protein IDs as keys and ProtT5 per-residue embeddings as values
    max_bkps_per100aa: an integer that represents the maximum number of boundaries per 100 amino acids in the protein

    Returns:
    protein_segments: a dictionary of Uniprot protein IDs as keys and a list of boundaries as values
    """

    protein_segments = {}

    # for each protein and embedding
    for protein, emb in emb_dict.items():
        n_boundaries = max(int(emb.shape[0]*max_bkps_per100aa/100),1)

        # try to get the breakpoints from change point analysis
        try:
            alg = rpt.Window(width=30, model='rbf', jump=1).fit(emb)
            # n_bkps is the maximum number of breakpoints allowed in this protein
            boundaries = alg.predict(n_bkps=n_boundaries)
            # convert boundaries to segments
            my_segments = [[0,boundaries[0]]]
            my_segments.extend([[boundaries[ii], boundaries[ii+1]] for ii in range(len(boundaries)-1)])
        # this algorithm can fail if the maximum number of breakpoints is unreasonable
        except:
            print("Failed", protein, "number of boundaries", n_boundaries)
            my_segments = "Failed"

        protein_segments[protein] = my_segments

    return protein_segments

def get_protein_segment_embeddings(emb_dict, protein_segment_boundaries):
    """
    Calculates the segment embedding for each segment defined by the change point analysis

    Parameters:
    emb_dict: a dictionary of UniProt protein IDs as keys and ProtT5 per-residue embeddings as values
    protein_segment_boundaries: a dictionary of Uniprot protein IDs as keys and a list of boundaries as values

    Returns:
    protein_segment_embeddings: a dictionary of "ID start-stop" as keys and segment embeddings as values
        IDs are UniProt protein IDs, start and stop corresponds to the positions of the segment in the 
        original protein sequence, and segment embeddings are vectors of size 1x1024
    """

    protein_segment_embeddings = {}

    # for each protein with an embedding
    for protein, emb in emb_dict.items():

        # if the segmentation did not fail
        if protein_segment_boundaries[protein] != "Failed":

            # generate a "segment embedding" for each segment of the protein
            for segment in protein_segment_boundaries[protein]:
                seg_emb = np.mean(emb[segment[0]:segment[1],:], axis=0)
                key = protein+" "+str(segment[0])+"-"+str(segment[1])

                protein_segment_embeddings[key] = seg_emb

    return protein_segment_embeddings


def print_info(avg_length, max_seq_len, n_long, seq_list):
    print('########################################')
    print('Example sequence: {}\n{}'.format(seq_list[0][0], seq_list[0][1]))
    print('########################################')
    print('Total number of sequences: {}'.format(len(seq_list)))
    print("Average sequence length: {}".format(avg_length))
    print("Number of sequences >{}: {}".format(max_seq_len, n_long))


def create_arg_parser():
    """"Creates and returns the ArgumentParser object."""

    # Instantiate the parser
    parser = argparse.ArgumentParser(description=(
            't5_embedder.py creates T5 embeddings for a given text ' +
            ' file containing sequence(s) in FASTA-format.'))

    # Required positional argument
    parser.add_argument('-i', '--input', required=True, type=str,
                        help='A path to a fasta-formatted text file containing protein sequence(s).')

    # Optional positional argument
    parser.add_argument('-o', '--output', required=True, type=str,
                        help='A path for saving the created embeddings as NumPy npz file.')

    # Required positional argument
    parser.add_argument('--model', required=False, type=str,
                        default=None,
                        help='A path to a directory holding the checkpoint for a pre-trained model')

    # Optional argument
    parser.add_argument('--per_protein', type=int,
                        default=0,
                        help="Whether to return per-residue embeddings (0: default) or the mean-pooled per-protein representation (1).")
    return parser

