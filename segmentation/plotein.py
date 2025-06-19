import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import ast
from collections import defaultdict
import matplotlib.patches as patches
import seaborn as sns

def load_segmentation_data(tsv_file):
    """
    Load protein segmentation data from TSV file.

    Args:
        tsv_file (str): Path to TSV file with uniprotID and segmentation columns

    Returns:
        dict: Dictionary with uniprotID as key and list of segments as value
    """
    # Read TSV file without header
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['uniprotID', 'segmentation'])

    segmentation_data = {}
    for _, row in df.iterrows():
        uniprot_id = row['uniprotID']
        # Parse the segmentation string as a Python list
        segments = ast.literal_eval(row['segmentation'])
        segmentation_data[uniprot_id] = segments

    return segmentation_data


def load_nes_data(csv_file):
    """
    Load NES motif data from CSV file.

    Args:
        csv_file (str): Path to CSV file with NES motif annotations

    Returns:
        pandas.DataFrame: DataFrame containing NES motif data
    """
    return pd.read_csv(csv_file)


def calculate_overlap_percentage(segment_start, segment_end, nes_start, nes_end):
    """
    Calculate the percentage of NES motif that overlaps with a segment.

    Args:
        segment_start (int): Start position of segment
        segment_end (int): End position of segment
        nes_start (int): Start position of NES motif
        nes_end (int): End position of NES motif

    Returns:
        float: Percentage of NES motif overlapping with segment (0-100)
    """
    overlap_start = max(segment_start, nes_start)
    overlap_end = min(segment_end, nes_end)

    if overlap_start >= overlap_end:
        return 0.0

    overlap_length = overlap_end - overlap_start
    nes_length = nes_end - nes_start

    return (overlap_length / nes_length) * 100 if nes_length > 0 else 0.0


def find_best_nes_segments(segmentation_data, nes_data):
    """
    Find the segment with the best NES overlap for each protein.

    Args:
        segmentation_data (dict): Protein segmentation data
        nes_data (pandas.DataFrame): NES motif data

    Returns:
        dict: Dictionary with uniprotID as key and list of best NES segment info as value
    """
    best_nes_segments = defaultdict(list)

    for _, nes_row in nes_data.iterrows():
        uniprot_id = nes_row['uniprotID']
        nes_start = nes_row['start']
        nes_end = nes_row['end']
        nes_id = nes_row['NESdb_ID']
        protein_name = nes_row['name']

        if uniprot_id not in segmentation_data:
            continue

        segments = segmentation_data[uniprot_id]
        best_overlap = 0
        best_segment_length = 0
        best_segment_idx = -1

        for i, segment in enumerate(segments):
            segment_start, segment_end = segment
            overlap_pct = calculate_overlap_percentage(segment_start, segment_end, nes_start, nes_end)

            if overlap_pct > best_overlap:
                best_overlap = overlap_pct
                best_segment_length = segment_end - segment_start
                best_segment_idx = i

        if best_segment_idx >= 0 and best_overlap > 0:
            best_nes_segments[uniprot_id].append({
                'segment_idx': best_segment_idx,
                'overlap_pct': best_overlap,
                'nes_start': nes_start,
                'nes_end': nes_end,
                'nes_id': nes_id,
                'best_segment_length': best_segment_length,
                'protein_name': protein_name
            })

    return best_nes_segments


def plot_protein_annotation(uniprot_id, protein_name, segments, nes_info, save_path=None, model_name=None):
    """
    Create a linear annotation plot for a single protein.

    Args:
        uniprot_id (str): UniProt ID of the protein
        segments (list): List of [start, end] segments
        nes_info (list): List of dictionaries containing NES information
        save_path (str, optional): Path to save the plot
    """
    fig, ax = plt.subplots(figsize=(12, 4))

    # Calculate protein length
    protein_length = max([seg[1] for seg in segments])

    # Draw main protein line
    ax.plot([0, protein_length], [0, 0], 'k-', linewidth=3, label='Protein sequence')

    # Get NES segment indices for highlighting
    nes_segment_indices = set()
    if nes_info:
        nes_segment_indices = {info['segment_idx'] for info in nes_info}

    # Draw segments
    bar_height = 0.3
    label = True
    for i, (start, end) in enumerate(segments):
        # Highlight NES-containing segments
        if i in nes_segment_indices:
            # Draw highlighted segment with border
            nes_rect = patches.Rectangle((start, -bar_height / 2), end - start, bar_height,
                                         linewidth=1, edgecolor='black', facecolor='#f39b7f',
                                         alpha=0.8,
                                         label='Best matching segment' if i == min(nes_segment_indices) else "")
            ax.add_patch(nes_rect)
        else:
            # Draw regular segment
            rect = patches.Rectangle((start, -bar_height / 2), end - start, bar_height,
                                 linewidth=1, edgecolor='black', facecolor='#4dbbd5', alpha=0.6,
                                 label='Embedding segment' if label else "")
            label = False
            ax.add_patch(rect)

    # Add NES motif annotations
    labled = True
    if nes_info:
        for i, info in enumerate(reversed(nes_info)):
            nes_start = info['nes_start']
            nes_end = info['nes_end']

            # Draw NES motif as a thick line
            height = -0.2 - (i * 0.05)  # Adjust height for multiple NES motifs
            ax.plot([nes_start, nes_end], [height, height], linestyle='-', color='#e64b35', linewidth=3,
                    label='Annotated NES motif' if labled else "")
            labled = False

    # Formatting
    ax.set_xlim(-protein_length * 0.05, protein_length * 1.05)
    ax.set_ylim(-0.5, 0.8)
    ax.set_xlabel('Amino Acid Position')
    ax.set_ylabel('')
    protein_name = protein_name.split(' ')[0]  # Use only the first part of the protein name
    matching_precent = sum(info['overlap_pct'] for info in nes_info) / len(nes_info) if nes_info else 0
    if matching_precent > 0:
        ax.text(0.5, 0.9, f'NES Motif Matching: {matching_precent:.1f}%',
                transform=ax.transAxes, ha='center', fontsize=10, color='#e64b35')
    ax.set_title(f'{model_name} Embedding Based Protein Annotation\n{protein_name} (Length: {protein_length} aa)')
    ax.grid(True, alpha=0.3)
    ax.set_yticks([])
    ax.legend(loc='upper right')
    plt.tight_layout()

    if save_path:
        plt.savefig(f"{save_path}/{model_name}_{uniprot_id}.png", dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def create_all_protein_plots(segmentation_data, nes_data, from_list=None, save_dir=None, model_name=None):
    """
    Create annotation plots for all proteins.

    Args:
        segmentation_data (dict): Protein segmentation data
        nes_data (pandas.DataFrame): NES motif data
        save_dir (str, optional): Directory to save plots
    """
    # Find best NES segments
    best_nes_segments = find_best_nes_segments(segmentation_data, nes_data)

    for uniprot_id, segments in segmentation_data.items():
        if (from_list is not None) and (uniprot_id not in from_list):
            continue
        nes_info = best_nes_segments.get(uniprot_id, [])
        if len(nes_info) < 1:
            continue

        print(f"Processing {uniprot_id}...")
        if nes_info:
            print(f"  Found {len(nes_info)} NES motif(s)")
            for info in nes_info:
                print(
                    f"    NES_{info['nes_id']}: {info['overlap_pct']:.1f}% overlap with segment {info['segment_idx'] + 1}")
        else:
            print(f"  No NES motif(s) found")
        save_path = save_dir if save_dir else None
        prot_name = nes_info[0]['protein_name']
        plot_protein_annotation(uniprot_id, prot_name, segments, nes_info, save_path, model_name)


def calculate_precent_of_match(segmentation_data, nes_data, top=10):
    """
    Calculate the percentage of NES motif matching for each protein.

    Args:
        segmentation_data (dict): Protein segmentation data
        nes_data (pandas.DataFrame): NES motif data
        top (int): Number of top matches to return
    """
    best_matches = []
    match_data = {}
    # Find best NES segments
    best_nes_segments = find_best_nes_segments(segmentation_data, nes_data)

    for uniprot_id, segments in segmentation_data.items():
        nes_info = best_nes_segments.get(uniprot_id, [])
        if len(nes_info) < 1:
            continue
        total_overlap = sum(info['overlap_pct'] for info in nes_info) / len(nes_info)
        match_data[uniprot_id] = {
            'percent': total_overlap,
            'nes_count': len(nes_info)
        }

    # Sort by percentage first, then by NES count
    sorted_proteins = sorted(
        match_data.items(),
        key=lambda item: (item[1]['nes_count']),
        reverse=True
    )

    # Print top proteins with highest NES motif matching percentage
    print(f"Top {top} proteins with highest NES motif matching count:")
    for i, (uniprot_id, data) in enumerate(sorted_proteins):
        if i >= top:
            break
        print(f"{i + 1}. {uniprot_id}: {data['percent']:.2f}% (NES count: {data['nes_count']})")
        best_matches.append(uniprot_id)

    # Calculate and print average percentage of NES motif matching
    average_percent = sum(data['percent'] for _, data in match_data.items()) / len(match_data)
    print(f"Average percentage of NES motif matching: {average_percent:.2f}%")

    return best_matches

def calculate_best_segment_length_distribution(segmentation_data, nes_data, save_dir=None, model_name=None):
    """
    Calculate the best segment length distribution for all proteins.

    Args:
        segmentation_data (dict): Protein segmentation data
        nes_data (pandas.DataFrame): NES motif data
        save_dir (str, optional): Directory to save plots
    """
    all_chosen_segments = []
    all_segments_lengths = []

    # Calculate the lengths of all segments
    for segments in segmentation_data.values():
        for start, end in segments:
            all_segments_lengths.append(end - start)

    # Find best NES segments
    best_nes_segments = find_best_nes_segments(segmentation_data, nes_data)
    for uniprot_id, segments in segmentation_data.items():
        nes_info = best_nes_segments.get(uniprot_id, [])
        # add all the lengths of segments
        for info in nes_info:
            all_chosen_segments.append(info['best_segment_length'])
    # plot the distribution of segment lengths kde style, using seaborn
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.histplot(all_chosen_segments, element="step", color='#f39b7f', alpha=0.3, label=f'Best Segment Lengths (n={len(all_chosen_segments)})',
                 kde=True, stat="density",)
    sns.histplot(all_segments_lengths, element="step", color='#4dbbd5', alpha=0.3, label=f'All Segments Lengths (n={len(all_segments_lengths)})',
                 kde=True, stat="density",)
    plt.legend()
    plt.title(f'Distribution of Best Segment Lengths for NES Motifs ({model_name})')
    plt.xlabel('Segment Length (amino acids)')
    plt.ylabel('Density')
    plt.xlim(0, max(all_chosen_segments) * 1.05)
    plt.xticks(np.arange(0, max(all_chosen_segments) + 1, step=10))
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    if save_dir:
        plt.savefig(f"{save_dir}/{model_name}_segment_length_distribution.png", dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()
