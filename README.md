# 3D Biology Hackathon 2025
In this project we present a pipeline for predicting NESs in proteins, it includes:
- Segmenting protein sequences using Zero-Shot Protein Segmentation (see Acknowledgements) in `segmentation/`
- Finding Exposed and Helical Regions in `protein_prediction/exposed_helical_extraction/`
- Training a model to predict NESs using ESM embeddings in `preperation/`
- Predicting NESs across entire protein sequences and cross referencing with annotations found above in `protein_prediction/`

## Installation

### Prerequisites

* **Python 3.8+**
* **Git**
* **CUDA-enabled GPU** (optional, for faster embedding)
* **DSSP** executable for secondary-structure analysis (see `protein_prediction/exposed_helical_extraction` for details) .

### Dependencies

Install via pip:

```bash
pip install -r requirements.txt
```

## Usage

Note: as all data and model files exist, any step can be run independently

### 1. Generate Embeddings & Segmentation

Run the full pipeline on NESDB data:
```bash
python3 divide_and_color.py
```
Output: `segmentation/data/ESM_NESDB_combined_segments.tsv`

#### 1.1 Visualize Results

The pipeline saves plots of embedding segments and NES overlaps under `protein_plots/`. For interactive plotting, import functions from `plotein.py`:

```python
from plotein import plot_protein_annotation, create_all_protein_plots
# load data and call plotting functions
```
### 2. Train Model

```bash
cd preperation
python ex4.py 
```
Output: `preperation/models/trained_model_layer_33.pth`

### 3. Extract Exposed and Helical Regions

Run the `protein_prediction/exposed_helical_extraction/exposed_helical_analysis_and_extraction.ipynb` notebook (further documentation inside)

Output: `protein_prediction/exposed_helical_extraction/data/residue_annotation.csv`


### 4. Predict NESs across proteins and cross-reference with annotations
Run the `protein_prediction/protein_pred.ipynb` notebook (further documentation inside)

Output: `protein_prediction/data/false_positive_seqs.csv`
## Authors

* **Noam Rosenmann** 
* **Shay Guttman** 
* **Noa Birman**
* **Dana Siton** 
* **Tal Neumann**
---

## Acknowledgments & References

* **ProtTrans**, **ESM-2**, **Ruptures**, **BioPython**
* Inspired by the NESDB database and Nuclear Export Signal analysis methods
* [MosesLab/ZeroShotProteinSegmentation](https://github.com/moses-lab/zero-shot-protein-segmentation.git)
* [ProtTrans/protT5](https://github.com/agemagician/ProtTrans)
* [FacebookResearch/esm](https://github.com/facebookresearch/esm)

