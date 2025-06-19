import h5py
from protT5_embedder import get_embeddings_from_csv
from esm_embedder import get_esm_embeddings_from_csv
from utils import *
import warnings
warnings.filterwarnings("ignore")

def process_protein_embeddings(
        csv_path,
        seg_bounds_path,
        save_whole_emb_to_hdf5=False,
        whole_emb_path=None,
        save_seg_emb_to_hdf5=False,
        seg_emb_path=None,
        model_dir="",
        per_protein=False,
        max_residues=4000,
        max_seq_len=4000,
        max_batch=100,
        max_bkps_per100aa=3,
        model_type='t5'
):
    """
    Process protein sequences to generate embeddings and segment boundaries.

    Args:
        csv_path (str): Path to the csv file containing protein sequences
        seg_bounds_path (str): Path to save the protein segment boundaries TSV file
        save_whole_emb_to_hdf5 (bool): Whether to save whole protein embeddings to HDF5
        whole_emb_path (str): Path for whole protein embeddings HDF5 file (required if save_whole_emb_to_hdf5=True)
        save_seg_emb_to_hdf5 (bool): Whether to save segment embeddings to HDF5
        seg_emb_path (str): Path for segment embeddings HDF5 file (required if save_seg_emb_to_hdf5=True)
        model_dir (str): Directory for the ProtT5 model
        per_protein (bool): Whether to process per protein
        max_residues (int): Maximum number of residues
        max_seq_len (int): Maximum sequence length
        max_batch (int): Maximum batch size
        max_bkps_per100aa (int): Maximum breakpoints per 100 amino acids for segmentation

    Returns:
        tuple: (emb_dict, protein_segments) - embeddings dictionary and protein segments

    Raises:
        ValueError: If required paths are not provided when saving options are enabled
    """

    # Validate inputs
    if save_whole_emb_to_hdf5 and not whole_emb_path:
        raise ValueError("whole_emb_path is required when save_whole_emb_to_hdf5 is True")

    if save_seg_emb_to_hdf5 and not seg_emb_path:
        raise ValueError("seg_emb_path is required when save_seg_emb_to_hdf5 is True")

    if model_type == 't5':
        # Generate ProtT5 embeddings
        emb_dict = get_embeddings_from_csv(
            csv_path=csv_path,
            model_dir=model_dir,
            per_protein=per_protein,
            max_residues=max_residues,
            max_seq_len=max_seq_len,
            max_batch=max_batch
        )
    elif model_type == 'esm':
        # Generate ESM embeddings
        emb_dict = get_esm_embeddings_from_csv(
            csv_path=csv_path,
            embedding_size=1280,  # ESM-2 embedding size
            per_protein=per_protein,
            max_residues=max_residues,
            max_seq_len=max_seq_len,
            max_batch=max_batch,
            print_seq_info=True,
            embedding_layer=33  # ESM-2 embedding layer
        )

    # Save per-residue embeddings of the proteins (optional)
    if save_whole_emb_to_hdf5:
        with h5py.File(str(whole_emb_path), "a") as hf:
            for sequence_id, embedding in emb_dict.items():
                hf.create_dataset(sequence_id, data=embedding)

    # Identify segment boundaries using change point analysis
    protein_segments = get_protein_segments(emb_dict, max_bkps_per100aa=max_bkps_per100aa)

    # Save segment boundaries
    with open(seg_bounds_path, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        for protein_id, protein_seg in protein_segments.items():
            writer.writerow([protein_id, protein_seg])

    # Make and save segment embeddings (optional)
    if save_seg_emb_to_hdf5:
        protein_segment_embeddings = get_protein_segment_embeddings(emb_dict, protein_segments)

        with h5py.File(str(seg_emb_path), "a") as hf:
            for sequence_key, embedding in protein_segment_embeddings.items():
                hf.create_dataset(sequence_key, data=embedding)

    return emb_dict, protein_segments


# Example usage:
if __name__ == "__main__":
    # Basic usage with original parameters
    # emb_dict, protein_segments = process_protein_embeddings(
    #     csv_path="NESDB_combined_database.csv",
    #     seg_bounds_path="NESDB_combined_segments.tsv",
    #     save_whole_emb_to_hdf5=False,
    #     whole_emb_path="NESDB_combined_whole_emb.hdf5",
    #     save_seg_emb_to_hdf5=False,
    #     seg_emb_path="NESDB_combined_seg_emb.hdf5"
    # )

    esm_emb_dict, esm_protein_segments = process_protein_embeddings(
        csv_path="NESDB_combined_database.csv",
        seg_bounds_path="ESM_NESDB_combined_segments.tsv",
        save_whole_emb_to_hdf5=False,
        whole_emb_path="NESDB_combined_whole_emb.hdf5",
        save_seg_emb_to_hdf5=False,
        seg_emb_path="NESDB_combined_seg_emb.hdf5",
        model_type='esm'
    )