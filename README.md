# collapse-motifs

Code accompanying the paper ["Unsupervised learning reveals landscape of local structural motifs across protein classes"](https://www.biorxiv.org/content/10.1101/2023.12.04.569990v1) by Alexander Derry and Russ B. Altman.

## Install dependencies

All dependencies can be installed by installing COLLAPSE following the instructions [here](https://github.com/awfderry/COLLAPSE/tree/main).

## Download fitted clustering model

Fitted clustering models for various number of clusters `k` are available for download [here](https://zenodo.org/records/10699466). Download into the `data` folder of this repo to avoid issues when importing in the scripts below.

The number of clusters controls the specificity of the resulting structural motifs; we recommend using `k=50000` for general use, and this is the value used in all analysis in our paper. 

## Embed a directory of PDB files using COLLAPSE

Given a set of protein structures located at `DATA_DIR`, this script will embed all residues with COLLAPSE and stored at `OUT_DIR` in [LMDB format](https://github.com/awfderry/COLLAPSE/tree/main?tab=readme-ov-file#embed-entire-dataset-of-pdb-files).

```python embed_pdb_dataset.py DATA_DIR OUT_DIR [--filetype] [--split_id] [--num_splits]```

Additional arguments are 
- `--filetype` (default `pdb`): file type of input structures, e.g. `pdb`,`pdb.gz`,`cif`

For large datasets, we recommend processing in parallel using the num_splits argument.

```python embed_pdb_dataset.py DATA_DIR OUT_DIR --split_id=$i --num_splits=NUM_SPLITS --filetype=pdb```

This produces NUM_SPLITS `tmp_` files in OUT_DIR. To combine all into the full dataset, run the following:

```python -m atom3d.datasets.scripts.combine_lmdb OUT_DIR/tmp_* OUT_DIR/full```

## Featurize proteins using pre-trained clustering model

The following script takes an LMDB dataset of protein structures embedded using COLLAPSE (see above) and computes the clusters associated with each embedded residue. The output of this script can be in one of two formats: (1) LMDB, which simply adds a new key `clusters` to each element and stores a new LMDB dataset, or (2) pickle, which saves a dictionary mapping from the protein id to the list of clusters in residue sequence order.

```python cluster_featurize.py DATA_PATH OUT_PATH [--out_fmt] [--k] [--split_id] [--num_splits]```

Optional arguments are
- `--out_fmt` (default `pkl`): Output format; valid options are either `lmdb` or `pkl`
- `--k` (default `50000`): Number of clusters from pre-fitted clustering model, assumed to be in `./data/pdb100_cluster_fit_{k}.pkl`
- `--split_id` and `--num_splits` are as above. Only used when `--out_fmt==lmdb`


## Compute TF-IDF fingerprints for clusterized dataset

This script computes TF-IDF fingerprints for a cluster-featurized dataset in `pkl` format, as produced by `cluster_featurize.py`.

```python cluster_tfidf.py cluster_sequences out_path [--id_to_label] [--k]```

Arguments are:
- `cluster_sequences`: dict mapping `id` to cluster sequence, saved in `pkl` format
- `out_path`: path where output representations will be saved (in `pkl` format with keys {`data`: TF-IDF fingerprints, `labels`: labels, `ids`: ids})
- `id_to_label` (optional): `pkl` file containing dict mapping from `id` (same as in `cluster_sequences`) to labels. If not give, all labels will be `None`.

## Scripts for reproducing analysis in paper

Analysis performed in the paper can be found in the following notebooks: `eval_clustering.ipynb`, `eval_fold_search.ipynb`, `klifs_analysis.ipynb`, and `mutations.ipynb`. Code for fitting clustering models is found in `run_clustering.py`.