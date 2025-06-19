import torch
import esm
import time
from tqdm import tqdm
from utils import read_csv_sequences, print_info

# All of ESM-2 pre-trained models by embedding size
ESM_MODELS_DICT = {320: esm.pretrained.esm2_t6_8M_UR50D,
                   480: esm.pretrained.esm2_t12_35M_UR50D,
                   640: esm.pretrained.esm2_t30_150M_UR50D,
                   1280: esm.pretrained.esm2_t33_650M_UR50D,
                   2560: esm.pretrained.esm2_t36_3B_UR50D,
                   5120: esm.pretrained.esm2_t48_15B_UR50D}


def get_esm_model(embedding_size=1280):
    """
    Retrieves a pre-trained ESM-2 model
    :param embedding_size: The ESM-2 model embedding size
    :return: esm_model, alphabet, batch_converter, device
    """

    if embedding_size not in ESM_MODELS_DICT:
        raise ValueError(f"ERROR: ESM does not have a trained model with embedding size of {embedding_size}.\n "
                         f"Please use one of the following embedding sized: {ESM_MODELS_DICT.keys()}")

    model, alphabet = ESM_MODELS_DICT[embedding_size]()
    batch_converter = alphabet.get_batch_converter()
    model.eval()  # disables dropout for deterministic results
    # check if GPU is available
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model.to(device)
    print(f"ESM model loaded to {device}")
    return model, alphabet, batch_converter, device


def get_esm_embeddings(pep_tuple_list, esm_model, alphabet, batch_converter, device, embedding_layer=33, sequence_embedding=True):
    """
    This function convert peptide sequence data into ESM sequence embeddings
    :param pep_tuple_list: peptide tuple list of format : [(name_1, seq_1), (name_2, seq_2), ...]
    :param esm_model: Pre-trained ESM-2 model
    :param alphabet: ESM-2 alphabet object
    :param batch_converter: ESM-2 batch_converter object
    :param device: GPU/CPU device
    :param embedding_layer: The desired embedding layer to get
    :param sequence_embedding: Whether to use a sequence embedding (default=True) or amino acid embedding
    :return: List of ESM-2 sequence/amino acids embeddings
    """
    batch_labels, batch_strs, batch_tokens = batch_converter(pep_tuple_list)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

    # Extract per-residue representations
    with torch.no_grad():
        results = esm_model(batch_tokens.to(device), repr_layers=[embedding_layer])
    token_representations = results["representations"][embedding_layer]

    # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
    representations = []
    for i, tokens_len in enumerate(batch_lens):
        embedding = token_representations[i, 1: tokens_len - 1]
        # Generate per-sequence representations via averaging
        if sequence_embedding:
            embedding = embedding.mean(dim=0, keepdim=False)
        representations.append(embedding.cpu().numpy())

    return representations

def get_esm_embeddings_from_csv(csv_path,
                   embedding_size=1280,  # ESM model embedding size instead of model_dir
                   per_protein=False,  # whether to derive per-protein embeddings
                   max_residues=4000,  # number of cumulative residues per batch
                   max_seq_len=4000,  # max length after which we switch to single-sequence processing to avoid OOM
                   max_batch=100,  # max number of sequences per single batch
                   print_seq_info=True,
                   embedding_layer=33):  # ESM embedding layer
    seq_dict = dict()
    emb_dict = dict()

    # Read in fasta
    seq_dict = read_csv_sequences(csv_path)
    model, alphabet, batch_converter, device = get_esm_model(embedding_size)

    avg_length = sum([len(seq) for _, seq in seq_dict.items()]) / len(seq_dict)
    n_long = sum([1 for _, seq in seq_dict.items() if len(seq) > max_seq_len])
    seq_list = sorted(seq_dict.items(), key=lambda kv: len(seq_dict[kv[0]]), reverse=True)

    if print_seq_info:
        print_info(avg_length, max_seq_len, n_long, seq_list)

    start = time.time()
    batch = list()
    for seq_idx, (pdb_id, seq) in enumerate(tqdm(seq_list, desc="Processing sequences", unit="seq"), 1):
        seq = seq.replace('U', 'X').replace('Z', 'X').replace('O', 'X')
        seq_len = len(seq)
        batch.append((pdb_id, seq, seq_len))

        # count residues in current batch and add the last sequence length to
        # avoid that batches with (n_res_batch > max_residues) get processed
        n_res_batch = sum([s_len for _, _, s_len in batch]) + seq_len
        if len(batch) >= max_batch or n_res_batch >= max_residues or seq_idx == len(seq_list) or seq_len > max_seq_len:
            # Prepare batch for ESM
            pdb_ids, seqs, seq_lens = zip(*batch)
            pep_tuple_list = [(pdb_id, seq) for pdb_id, seq in zip(pdb_ids, seqs)]
            batch = list()

            try:
                # Get ESM embeddings
                embeddings = get_esm_embeddings(pep_tuple_list, model, alphabet, batch_converter,
                                              device, embedding_layer, sequence_embedding=per_protein)
            except RuntimeError:
                print("RuntimeError during embedding for {} (L={}). Try lowering batch size. ".format(pdb_id, seq_len) +
                      "If single sequence processing does not work, you need more vRAM to process your protein.")
                continue

            # Store embeddings
            for batch_idx, identifier in enumerate(pdb_ids):
                s_len = seq_lens[batch_idx]
                emb = embeddings[batch_idx]

                if len(emb_dict) == 0:
                    print("Embedded protein {} with length {} to emb. of shape: {}".format(
                        identifier, s_len, emb.shape))

                emb_dict[identifier] = emb.squeeze()

    end = time.time()

    print('\n############# STATS #############')
    print('Total time: {:.2f}[s]; time/prot: {:.4f}[s]; avg. len= {:.2f}'.format(
        end - start, (end - start) / len(emb_dict), avg_length))
    return emb_dict




