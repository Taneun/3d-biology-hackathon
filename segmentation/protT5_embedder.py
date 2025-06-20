#https://github.com/agemagician/ProtTrans/blob/master/Embedding/prott5_embedder.py
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 18:33:22 2020

@author: mheinzinger (edits by Ami G. Sangster)


"""
import time
from pathlib import Path
import torch
from transformers import T5EncoderModel, T5Tokenizer
from tqdm import tqdm
from utils import read_csv_sequences, print_info

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print("Using device: {}".format(device))

def get_T5_model(model_dir, transformer_link = "Rostlab/prot_t5_xl_half_uniref50-enc"):
    print("Loading: {}".format(transformer_link))
    if model_dir is not None:
        print("##########################")
        print("Loading cached model from: {}".format(model_dir))
        print("##########################")
    model = T5EncoderModel.from_pretrained(transformer_link)#, cache_dir=model_dir, from_tf=True)
    model.full() if device=='cpu' else model.half() # only cast to full-precision if no GPU is available

    model = model.to(device)
    model = model.eval()
    vocab = T5Tokenizer.from_pretrained(transformer_link, do_lower_case=False )
    return model, vocab


def read_fasta( fasta_path ):
    '''
        Reads in fasta file containing multiple sequences.
        Returns dictionary of holding multiple sequences or only single 
        sequence, depending on input file.
    '''
    
    sequences = dict()
    with open( fasta_path, 'r' ) as fasta_f:
        for line in fasta_f:
            # get uniprot ID from header and create new entry
            if line.startswith('>'):
                uniprot_id = line.replace('>', '').strip()
                # replace tokens that are mis-interpreted when loading h5
                uniprot_id = uniprot_id.replace("/","_").replace(".","_")
                sequences[ uniprot_id ] = ''
            else:
                # repl. all whie-space chars and join seqs spanning multiple lines
                sequences[ uniprot_id ] += ''.join( line.split() ).upper().replace("-","") # drop gaps and cast to upper-case
                
    return sequences



def read_fasta_chunks( fasta_path, chunk_size=1000, chunk_iteration=0 ):
    '''
        Reads in fasta file containing multiple sequences.
        Returns dictionary of holding multiple sequences or only single 
        sequence, depending on input file.

        edited above function read in 'chunk_size' number of proteins from 
        the fasta file. If 'chunk_interation' is 0 then this reads the first 
        'chunk_size' sequences, if its 1 then it reads the next 'chunk_size'
        sequences (or less if its the end of the file)
    '''
    
    sequences = dict()
    protein_ii = 0
    protein_ii_start = chunk_size*chunk_iteration
    protein_ii_stop = chunk_size*(chunk_iteration+1)
    with open( fasta_path, 'r' ) as fasta_f:
        for line in fasta_f:
            # get uniprot ID from header and create new entry
            if protein_ii>protein_ii_start and protein_ii<protein_ii_stop:
                if line.startswith('>'):
                    uniprot_id = line.replace('>', '').strip()
                    # replace tokens that are mis-interpreted when loading h5
                    uniprot_id = uniprot_id.replace("/","_").replace(".","_")
                    sequences[ uniprot_id ] = ''
                    protein_ii+=1
                else:
                    # repl. all whie-space chars and join seqs spanning multiple lines
                    sequences[ uniprot_id ] += ''.join( line.split() ).upper().replace("-","") # drop gaps and cast to upper-case
            elif line.startswith('>'):
                protein_ii+=1
    return sequences


def get_embeddings_from_fasta(seq_path,
                   model_dir,
                   per_protein, # whether to derive per-protein (mean-pooled) embeddings
                   max_residues=4000, # number of cumulative residues per batch
                   max_seq_len=4000, # max length after which we switch to single-sequence processing to avoid OOM
                   max_batch=100 # max number of sequences per single batch
                ):

    seq_dict = dict()
    emb_dict = dict()

    # Read in fasta
    seq_dict = read_fasta( seq_path )
    model, vocab = get_T5_model(model_dir)

    print('########################################')
    print('Example sequence: {}\n{}'.format( next(iter(
            seq_dict.keys())), next(iter(seq_dict.values()))) )
    print('########################################')
    print('Total number of sequences: {}'.format(len(seq_dict)))

    avg_length = sum([ len(seq) for _, seq in seq_dict.items()]) / len(seq_dict)
    n_long     = sum([ 1 for _, seq in seq_dict.items() if len(seq)>max_seq_len])
    seq_dict   = sorted( seq_dict.items(), key=lambda kv: len( seq_dict[kv[0]] ), reverse=True )

    print("Average sequence length: {}".format(avg_length))
    print("Number of sequences >{}: {}".format(max_seq_len, n_long))

    start = time.time()
    batch = list()
    for seq_idx, (pdb_id, seq) in enumerate(seq_dict,1):
        seq = seq.replace('U','X').replace('Z','X').replace('O','X')
        seq_len = len(seq)
        seq = ' '.join(list(seq))
        batch.append((pdb_id,seq,seq_len))

        # count residues in current batch and add the last sequence length to
        # avoid that batches with (n_res_batch > max_residues) get processed
        n_res_batch = sum([ s_len for  _, _, s_len in batch ]) + seq_len
        if len(batch) >= max_batch or n_res_batch>=max_residues or seq_idx==len(seq_dict) or seq_len>max_seq_len:
            pdb_ids, seqs, seq_lens = zip(*batch)
            batch = list()

            token_encoding = vocab.batch_encode_plus(seqs, add_special_tokens=True, padding="longest")
            input_ids      = torch.tensor(token_encoding['input_ids']).to(device)
            attention_mask = torch.tensor(token_encoding['attention_mask']).to(device)

            try:
                with torch.no_grad():
                    embedding_repr = model(input_ids, attention_mask=attention_mask)
            except RuntimeError:
                print("RuntimeError during embedding for {} (L={}). Try lowering batch size. ".format(pdb_id, seq_len) +
                      "If single sequence processing does not work, you need more vRAM to process your protein.")
                continue

            # batch-size x seq_len x embedding_dim
            # extra token is added at the end of the seq
            for batch_idx, identifier in enumerate(pdb_ids):
                s_len = seq_lens[batch_idx]
                # slice-off padded/special tokens
                emb = embedding_repr.last_hidden_state[batch_idx,:s_len]

                if per_protein:
                    emb = emb.mean(dim=0)

                if len(emb_dict) == 0:
                    print("Embedded protein {} with length {} to emb. of shape: {}".format(
                        identifier, s_len, emb.shape))

                emb_dict[ identifier ] = emb.detach().cpu().numpy().squeeze()

    end = time.time()

    print('\n############# STATS #############')
    print('Total time: {:.2f}[s]; time/prot: {:.4f}[s]; avg. len= {:.2f}'.format(
            end-start, (end-start)/len(emb_dict), avg_length))
    return emb_dict


def get_embeddings_from_csv(csv_path,
                   model_dir,
                   per_protein,  # whether to derive per-protein (mean-pooled) embeddings
                   max_residues=4000,  # number of cumulative residues per batch
                   max_seq_len=4000,  # max length after which we switch to single-sequence processing to avoid OOM
                   max_batch=100,  # max number of sequences per single batch
                   print_seq_info=True):
    seq_dict = dict()
    emb_dict = dict()

    # Read in fasta
    seq_dict = read_csv_sequences(csv_path)
    model, vocab = get_T5_model(model_dir)

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
        seq = ' '.join(list(seq))
        batch.append((pdb_id, seq, seq_len))

        # count residues in current batch and add the last sequence length to
        # avoid that batches with (n_res_batch > max_residues) get processed
        n_res_batch = sum([s_len for _, _, s_len in batch]) + seq_len
        if len(batch) >= max_batch or n_res_batch >= max_residues or seq_idx == len(seq_list) or seq_len > max_seq_len:
            pdb_ids, seqs, seq_lens = zip(*batch)
            batch = list()

            token_encoding = vocab.batch_encode_plus(seqs, add_special_tokens=True, padding="longest")
            input_ids = torch.tensor(token_encoding['input_ids']).to(device)
            attention_mask = torch.tensor(token_encoding['attention_mask']).to(device)

            try:
                with torch.no_grad():
                    embedding_repr = model(input_ids, attention_mask=attention_mask)
            except RuntimeError:
                print("RuntimeError during embedding for {} (L={}). Try lowering batch size. ".format(pdb_id, seq_len) +
                      "If single sequence processing does not work, you need more vRAM to process your protein.")
                continue

            # batch-size x seq_len x embedding_dim
            # extra token is added at the end of the seq
            for batch_idx, identifier in enumerate(pdb_ids):
                s_len = seq_lens[batch_idx]
                # slice-off padded/special tokens
                emb = embedding_repr.last_hidden_state[batch_idx, :s_len]

                if per_protein:
                    emb = emb.mean(dim=0)

                if len(emb_dict) == 0:
                    print("Embedded protein {} with length {} to emb. of shape: {}".format(
                        identifier, s_len, emb.shape))

                emb_dict[identifier] = emb.detach().cpu().numpy().squeeze()

    end = time.time()

    print('\n############# STATS #############')
    print('Total time: {:.2f}[s]; time/prot: {:.4f}[s]; avg. len= {:.2f}'.format(
        end - start, (end - start) / len(emb_dict), avg_length))
    return emb_dict






