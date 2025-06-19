# ProtSeg: Protein Embedding & NES Motif Analysis Pipeline
---

## Summary

ProtSeg is a Python script for extracting protein sequences, generating per-residue embeddings with ProtT5 or ESM-2 models, segmenting proteins via change-point analysis, and analyzing overlap with Nuclear Export Signal (NES) motifs. Designed for computational biologists and bioinformaticians, ProtSeg automates the end-to-end workflow—from PDB structure to interactive plots—so you can rapidly identify regions most likely to contain NES motifs based on embedding features and visualize results.

---

## Table of Contents

1. [Quick Start](#quick-start)  
2. [Installation](#installation)  
3. [Usage](#usage)  
4. [Feature Highlights](#feature-highlights)   
7. [Advanced Usage & Integrations](#advanced-usage--integrations)   
11. [License](#license)  
12. [Authors](#authors)  
13. [Acknowledgments & References](#acknowledgments--references)

---

## Quick Start

Getting ProtSeg up and running in 60 seconds:

```bash
git clone https://github.com/username/protseg.git
cd protseg
pip install -r requirements.txt
```
---

## Installation

### Prerequisites

* **Python 3.8+**
* **Git**
* **CUDA-enabled GPU** (optional, for faster embedding)
* **DSSP** executable for secondary-structure analysis (see `Buried_Exposed_Ordered_Disordered.py` for details) .

### Dependencies

Install via pip:

```bash
pip install -r requirements.txt
```

Key packages:

* `biopython`
* `transformers`, `torch`, `tqdm`
* `pandas`, `numpy`
* `matplotlib`, `seaborn`
* `ruptures`, `h5py`
---

## Usage

### 1. Extract Sequence from PDB

```bash
python extract_seq_from_pdb.py path/to/structure.pdb output.fasta
```

This script parses the first chain in a PDB file and writes a FASTA sequence .

### 2. Generate Embeddings & Segmentation

Run the full pipeline on NESDB data:
```bash
python3 divide_and_color.py
```

### 3. Visualize Results

The pipeline saves plots of embedding segments and NES overlaps under `protein_plots/`. For interactive plotting, import functions from `plotein.py`:

```python
from plotein import plot_protein_annotation, create_all_protein_plots
# load data and call plotting functions
```

---

## Feature Highlights

- **PDB → FASTA Extraction** (`extract_seq_from_pdb.py`)
- **ProtT5 Embeddings** via `protT5_embedder.py`   
- **ESM-2 Embeddings** via `esm_embedder.py` 
- **Change-Point Segmentation** with Ruptures (`utils.get_protein_segments`) 
- **NES Overlap Analysis** (`plotein.find_best_nes_segments`) 
- **Publication-Quality Plots** (`plotein.plot_protein_annotation`) 
---

## Advanced Usage & Integrations

- **Plugin New Models**: Implement `get_<model>_embeddings_from_csv` in `divide_and_color.py`.  
- **Jupyter Notebooks**: Use functions from `utils.py` and `plotein.py` for interactive analysis.  

## License

This project is licensed under the [MIT License](LICENSE).

---

## Authors

* **Noam Rosenmann** – [GitHub](https://github.com/Nowam) 
* **Shay Guttman** – [GitHub](https://github.com/shayGuttmann)
* **Noa Birman** – [GitHub](https://github.com/noabirman)
* **Dana Siton** – [GitHub](https://github.com/danasiton)
* **Tal Neumann** – [GitHub](https://github.com/Taneun) 
---

## Acknowledgments & References

* **ProtTrans**, **ESM-2**, **Ruptures**, **BioPython**
* Inspired by the NESDB database and Nuclear Export Signal analysis methods
* [MosesLab/ZeroShotProteinSegmentation](https://github.com/moses-lab/zero-shot-protein-segmentation.git)
* [ProtTrans/protT5](https://github.com/agemagician/ProtTrans)
* [FacebookResearch/esm](https://github.com/facebookresearch/esm)

