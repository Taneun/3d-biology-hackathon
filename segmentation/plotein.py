import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import ast
from collections import defaultdict
import matplotlib.patches as patches


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


def calculate_percentage_of_segment(segment_start, segment_end, nes_start, nes_end):
    """
    Calculate the percentage of segment which is NES motif

    Args:
        segment_start (int): Start position of segment
        segment_end (int): End position of segment
        nes_start (int): Start position of NES motif
        nes_end (int): End position of NES motif

    Returns:
        float: Percentage of the segment which is NES motif (0-100)
    """
    overlap_start = max(segment_start, nes_start)
    overlap_end = min(segment_end, nes_end)

    if overlap_start >= overlap_end:
        return 0.0

    overlap_length = overlap_end - overlap_start
    segment_length = segment_end - segment_start

    return (overlap_length / segment_length) * 100 if segment_length > 0 else 0.0


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
        best_overlap_precent_of_segment = 0
        best_segment_idx = -1

        for i, segment in enumerate(segments):
            segment_start, segment_end = segment
            overlap_pct = calculate_overlap_percentage(segment_start, segment_end, nes_start, nes_end)

            if overlap_pct > best_overlap:
                best_overlap = overlap_pct
                best_overlap_precent_of_segment = calculate_percentage_of_segment(segment_start, segment_end, nes_start,
                                                                                  nes_end)
                best_segment_idx = i

        if best_segment_idx >= 0 and best_overlap > 0:
            best_nes_segments[uniprot_id].append({
                'segment_idx': best_segment_idx,
                'overlap_pct': best_overlap,
                'overlap_pct_of_segment': best_overlap_precent_of_segment,
                'nes_start': nes_start,
                'nes_end': nes_end,
                'nes_id': nes_id,
                'protein_name': protein_name
            })

    return best_nes_segments


def plot_protein_annotation(uniprot_id, protein_name, segments, nes_info, save_path=None):
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
    for i, (start, end) in enumerate(segments):
        # Highlight NES-containing segments
        if i in nes_segment_indices:
            # Draw highlighted segment with border
            nes_rect = patches.Rectangle((start, -bar_height / 2), end - start, bar_height,
                                         linewidth=2, edgecolor='red', facecolor='bisque',
                                         alpha=0.8,
                                         label='NES-containing segment' if i == min(nes_segment_indices) else "")
            ax.add_patch(nes_rect)
        else:
            # Draw regular segment
            rect = patches.Rectangle((start, -bar_height / 2), end - start, bar_height,
                                 linewidth=1, edgecolor='black', facecolor='lightsteelblue', alpha=0.6,
                                 label='Embedding segment' if i == 0 else "")

            ax.add_patch(rect)

        # Add segment labels
        # segment_center = (start + end) / 2
        # ax.text(segment_center, bar_height / 2 + 0.1, f'Seg_{i + 1}',
        #         ha='center', va='bottom', fontsize=8, rotation=0)

    # Add NES motif annotations
    labled = True
    if nes_info:
        for info in nes_info:
            nes_start = info['nes_start']
            nes_end = info['nes_end']
            nes_id = info['nes_id']

            # Draw NES motif as a thick line
            ax.plot([nes_start, nes_end], [-0.6, -0.6], 'r-', linewidth=4, alpha=0.8,
                    label='NES motif' if labled else "")
            labled = False

            # Add NES label
            # nes_center = (nes_start + nes_end) / 2
            # ax.text(nes_center, -0.8, f'NES_{nes_id}', ha='center', va='top',
            #         fontsize=8, color='red', weight='bold')

    # Formatting
    ax.set_xlim(-protein_length * 0.05, protein_length * 1.05)
    ax.set_ylim(-1.2, 0.8)
    ax.set_xlabel('Amino Acid Position')
    ax.set_ylabel('')
    protein_name = protein_name.split(' ')[0]  # Use only the first part of the protein name
    matching_precent = sum(info['overlap_pct_of_segment'] for info in nes_info) / len(nes_info) if nes_info else 0
    if matching_precent > 0:
        ax.text(0.5, 0.9, f'NES Motif Matching: {matching_precent:.1f}%',
                transform=ax.transAxes, ha='center', fontsize=10, color='darkred')
    ax.set_title(f'Protein Annotation: {protein_name} (Length: {protein_length} aa)')
    ax.grid(True, alpha=0.3)

    # Remove y-axis ticks
    ax.set_yticks([])

    # Add legend
    ax.legend(loc='upper right')

    plt.tight_layout()

    if save_path:
        plt.savefig(f"{save_path}/{uniprot_id}.png", dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def create_all_protein_plots(segmentation_data, nes_data, from_list=None, save_dir=None):
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
        plot_protein_annotation(uniprot_id, prot_name, segments, nes_info, save_path)

def calculate_precent_of_match(segmentation_data, nes_data, top=10):
    """
    Calculate the percentage of NES motif matching for each protein.

    Args:
        segmentation_data (dict): Protein segmentation data
        nes_data (pandas.DataFrame): NES motif data
        save_dir (str, optional): Directory to save plots
    """
    best_matches = []
    precent_of_match = dict()
    # Find best NES segments
    best_nes_segments = find_best_nes_segments(segmentation_data, nes_data)

    for uniprot_id, segments in segmentation_data.items():
        nes_info = best_nes_segments.get(uniprot_id, [])
        if len(nes_info) <= 1:
            continue
        total_overlap = sum(info['overlap_pct_of_segment'] for info in nes_info) / len(nes_info)
        precent_of_match[uniprot_id] = total_overlap
    # sort the dictionary by values in descending order
    sorted_precent_of_match = dict(sorted(precent_of_match.items(), key=lambda item: item[1], reverse=True))
    # print top 10 proteins with highest NES motif matching percentage
    print(f"Top {top} proteins with highest NES motif matching percentage:")
    for i, (uniprot_id, percent) in enumerate(sorted_precent_of_match.items()):
        if i >= top:
            break
        print(f"{i + 1}. {uniprot_id}: {percent:.2f}%")
        best_matches.append(uniprot_id)
    # print average percentage of NES motif matching
    average_percent = sum(sorted_precent_of_match.values()) / len(sorted_precent_of_match)
    print(f"Average percentage of NES motif matching: {average_percent:.2f}%")
    return best_matches

def main():
    """
    Main function to run the protein annotation plotting pipeline.
    """
    # File paths - modify these according to your file locations
    tsv_file = "NESDB_combined_segments.tsv"  # Replace with your TSV file path
    csv_file = "NESDB_combined_database.csv"  # Replace with your CSV file path
    save_directory = "protein_plots"  # Optional: directory to save plots
    import os
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)

    print("Loading segmentation data...")
    segmentation_data = load_segmentation_data(tsv_file)
    print(f"Loaded data for {len(segmentation_data)} proteins")

    print("Loading NES motif data...")
    nes_data = load_nes_data(csv_file)
    print(f"Loaded {len(nes_data)} NES motifs")

    best_nes = calculate_precent_of_match(segmentation_data, nes_data)
    # take top 10 proteins with highest NES motif matching percentage
    create_all_protein_plots(segmentation_data, nes_data, from_list=best_nes, save_dir=save_directory)


if __name__ == "__main__":
    main()
