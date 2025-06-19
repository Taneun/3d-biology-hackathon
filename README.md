````markdown
# ProtSeg: Protein Embedding & NES Motif Analysis Pipeline

[![Build Status](https://img.shields.io/github/actions/workflow/status/username/protseg/ci.yml?branch=main)](https://github.com/username/protseg/actions)  
[![PyPI Version](https://img.shields.io/pypi/v/protseg)](https://pypi.org/project/protseg)  
[![License](https://img.shields.io/github/license/username/protseg)](./LICENSE)  
[![Python Versions](https://img.shields.io/pypi/pyversions/protseg)](https://pypi.org/project/protseg)

---

## Executive Summary

ProtSeg is a command-line toolkit and Python library for extracting protein sequences, generating per-residue embeddings with ProtT5 or ESM-2 models, segmenting proteins via change-point analysis, and analyzing overlap with Nuclear Export Signal (NES) motifs. Designed for computational biologists and bioinformaticians, ProtSeg automates the end-to-end workflow—from PDB structure to interactive plots—so you can rapidly identify regions most likely to contain NES motifs based on embedding features and visualize results :contentReference[oaicite:0]{index=0}.

---

## Visual Introduction

![ProtSeg Demo](docs/assets/demo.gif)

---

## Table of Contents

1. [Quick Start](#quick-start)  
2. [Installation](#installation)  
3. [Usage](#usage)  
4. [Feature Highlights](#feature-highlights)  
5. [Configuration & Customization](#configuration--customization)  
6. [Troubleshooting & FAQ](#troubleshooting--faq)  
7. [Advanced Usage & Integrations](#advanced-usage--integrations)  
8. [Contributing](#contributing)  
9. [Testing](#testing)  
10. [Roadmap & Known Issues](#roadmap--known-issues)  
11. [License](#license)  
12. [Authors & Contact](#authors--contact)  
13. [Acknowledgments & References](#acknowledgments--references)

---

## Quick Start

Getting ProtSeg up and running in 60 seconds:

```bash
git clone https://github.com/username/protseg.git
cd protseg
pip install -r requirements.txt

# Example: full pipeline on NESDB data
python divide_and_color.py
````

---

## Installation

### Prerequisites

* **Python 3.8+**
* **Git**
* **CUDA-enabled GPU** (optional, for faster embedding)
* **DSSP** executable for secondary-structure analysis (see \[Buried\_Exposed\_Ordered\_Disordered.py] for details) .

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

### Platform-Specific Notes

* **Linux/macOS**: Ensure `mkdssp` is on `PATH` or set `dssp_exe` in scripts.
* **Windows**: Use WSL for DSSP support or conda package: `conda install -c salilab dssp`.

---

## Usage

### 1. Extract Sequence from PDB

```bash
python extract_seq_from_pdb.py path/to/structure.pdb output.fasta
```

This script parses the first chain in a PDB file and writes a FASTA sequence .

### 2. Generate Embeddings & Segmentation

Run the main pipeline:

```bash
python divide_and_color.py \
  --csv data/NESDB_combined_database.csv \
  --model t5 \
  --seg-bounds data/T5_segments.tsv \
  --per-protein 0
```

Options:

* `--model`: `t5` (ProtT5) or `esm` (ESM-2)
* `--per-protein`: `0` (residue embeddings) or `1` (mean-pooled)
* Change CSV and TSV paths as needed .

### 3. Visualize Results

The pipeline saves plots of embedding segments and NES overlaps under `protein_plots/`. For interactive plotting, import functions from `plotein.py`:

````python
from plotein import plot_protein_annotation, create_all_protein_plots
# load data and call plotting functions
``` :contentReference[oaicite:4]{index=4}.

---

## Feature Highlights

- **PDB → FASTA Extraction** (`extract_seq_from_pdb.py`) :contentReference[oaicite:5]{index=5}  
- **ProtT5 Embeddings** via `protT5_embedder.py`   
- **ESM-2 Embeddings** via `esm_embedder.py` :contentReference[oaicite:6]{index=6}  
- **Change-Point Segmentation** with Ruptures (`utils.get_protein_segments`) :contentReference[oaicite:7]{index=7}  
- **NES Overlap Analysis** (`plotein.find_best_nes_segments`) :contentReference[oaicite:8]{index=8}  
- **Publication-Quality Plots** (`plotein.plot_protein_annotation`) :contentReference[oaicite:9]{index=9}  

---

## Configuration & Customization

All parameters can be overridden at runtime:

- **Segmentation**: `--max-bkps-per100aa`  
- **Batching**: `--max-residues`, `--max-seq-len`, `--max-batch`  
- **Model Cache**: `--model-dir` (for ProtT5)  
- **Output Paths**: `--whole-emb-path`, `--seg-emb-path`

Consult the docstring of `process_protein_embeddings` in [divide_and_color.py] for full details :contentReference[oaicite:10]{index=10}.

---

## Troubleshooting & FAQ

**Q:** *I get a `RuntimeError` OOM during embedding.*  
**A:** Lower `--max-batch` or `--max-seq-len`, or move to GPU.

**Q:** *DSSP not found?*  
**A:** Install via `conda install -c salilab dssp` or build from source, then set `dssp_exe` path in `Buried_Exposed_Ordered_Disordered.py` :contentReference[oaicite:11]{index=11}.

**Q:** *Plots not saving?*  
**A:** Ensure `protein_plots/` directory exists or specify `--save-dir`.

---

## Advanced Usage & Integrations

- **Plugin New Models**: Implement `get_<model>_embeddings_from_csv` in `divide_and_color.py`.  
- **Jupyter Notebooks**: Use functions from `utils.py` and `plotein.py` for interactive analysis.  
- **Docker**: Community-provided Dockerfile available in `docker/` (coming soon).

---

## Contributing

1. Fork the repository  
2. Create a feature branch: `git checkout -b feature/<your-feature>`  
3. Write tests and adhere to PEP8  
4. Commit with [Conventional Commits]  
5. Open a Pull Request against `main`

See [CONTRIBUTING.md] for full guidelines.

---

## Testing

```bash
# (Assuming a tests/ directory exists)
pytest tests/  # run all unit tests
````

Expected output: all tests pass with 100% coverage.

---

## Roadmap & Known Issues

* **Roadmap**

  * Support for additional language models (e.g., ProtBERT)
  * Interactive web dashboard
  * GPU-accelerated segmentation

* **Known Issues**

  * Very short proteins may fail segmentation (ruptures exception)
  * DSSP mismatch on non-standard residues
  * High memory use for very long sequences

---

## License

This project is licensed under the [MIT License](LICENSE).

---

## Authors & Contact

* **Ami Sangster** – [GitHub](https://github.com/username) – [ami@example.com](mailto:ami@example.com)
* **Contributor Name** – [GitHub](https://github.com/contributor) – [contributor@example.com](mailto:contributor@example.com)

---

## Acknowledgments & References

* **ProtTrans**, **ESM-2**, **Ruptures**, **BioPython**
* Inspired by the NESDB database and Nuclear Export Signal analysis methods
* [ProtTrans/protT5](https://github.com/agemagician/ProtTrans)
* [FacebookResearch/esm](https://github.com/facebookresearch/esm)

```
```
