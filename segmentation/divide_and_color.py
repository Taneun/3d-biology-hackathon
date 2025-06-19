import h5py
from protT5_embedder import get_embeddings_from_csv
from esm_embedder import get_esm_embeddings_from_csv
from utils import *
from plotein import *
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

def main():
    """
    Main function to run the protein annotation plotting pipeline.
    """
    import os
    if not os.path.exists("data"):
        os.makedirs("data")
    # Process ProtT5 embeddings
    print("############# ProT5 #############")
    emb_dict, protein_segments = process_protein_embeddings(
        csv_path="data/NESDB_combined_database.csv",
        seg_bounds_path="data/T5_NESDB_combined_segments.tsv",
        save_whole_emb_to_hdf5=False,
        whole_emb_path="data/T5_NESDB_combined_whole_emb.hdf5",
        save_seg_emb_to_hdf5=False,
        seg_emb_path="data/T5_NESDB_combined_seg_emb.hdf5"
    )
    # Process ESM embeddings
    print("############# ESM-2 #############")
    esm_emb_dict, esm_protein_segments = process_protein_embeddings(
        csv_path="data/NESDB_combined_database.csv",
        seg_bounds_path="data/ESM_NESDB_combined_segments.tsv",
        save_whole_emb_to_hdf5=False,
        whole_emb_path="data/ESM_NESDB_combined_whole_emb.hdf5",
        save_seg_emb_to_hdf5=False,
        seg_emb_path="data/ESM_NESDB_combined_seg_emb.hdf5",
        model_type='esm'
    )
    # File paths - modify these according to your file locations
    t5_tsv_file = "data/T5_NESDB_combined_segments.tsv"
    esm_tsv_file = "data/ESM_NESDB_combined_segments.tsv"
    t5_model_name = "ProtT5"
    esm_model_name = "ESM-2"
    csv_file = "data/NESDB_combined_database.csv"
    save_directory = "protein_plots"
    import os
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
    t5_segmentation_data = load_segmentation_data(t5_tsv_file)
    esm_segmentation_data = load_segmentation_data(esm_tsv_file)
    nes_data = load_nes_data(csv_file)

    t5_best_nes = calculate_precent_of_match(t5_segmentation_data, nes_data)
    esm_best_nes = calculate_precent_of_match(esm_segmentation_data, nes_data)
    # if there are proteins that are in both lists
    common_proteins = set(t5_best_nes) & set(esm_best_nes)
    chosen = ["Q06219", "O00255"]

    print("############# ProT5 data analysis #############")
    calculate_best_segment_length_distribution(t5_segmentation_data, nes_data, save_dir=save_directory, model_name=t5_model_name)
    # take top 10 proteins with highest NES motif matching percentage
    create_all_protein_plots(t5_segmentation_data, nes_data, from_list=chosen, save_dir=save_directory, model_name=t5_model_name)

    print("############# ESM-2 data analysis #############")
    calculate_best_segment_length_distribution(esm_segmentation_data, nes_data, save_dir=save_directory, model_name=esm_model_name)
    # take top 10 proteins with highest NES motif matching percentage
    create_all_protein_plots(esm_segmentation_data, nes_data, from_list=chosen, save_dir=save_directory, model_name=esm_model_name)

if __name__ == "__main__":
    main()
