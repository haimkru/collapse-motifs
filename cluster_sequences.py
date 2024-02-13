import sys
import numpy as np
import pandas as pd
import pickle
import collections as col
from tqdm import tqdm
from atom3d.datasets import load_dataset, make_lmdb_dataset
from collapse.utils import serialize, deserialize
from collapse.atom_info import aa_abbr

alphabet = dict(zip(range(len(aa_abbr)), aa_abbr))

k = sys.argv[1]

data_path = '/scratch/users/aderry/scop40/scop40_208_lmdb/full'

clusterer = deserialize(f'data/pdb100_cluster_fit_{k}.pkl')

dataset = load_dataset(data_path, 'lmdb')

sid_to_seq = {}
sid_to_label = {}
skip = 0
for item in tqdm(dataset):
    if item is None:
        skip += 1
        continue

    sid = item['id'].split('.')[0]
    emb = item['embeddings']
    clusters = clusterer.predict(emb)
    sid_to_seq[sid] = clusters
    # sid_to_label[sid] = item['label']
        
print(f'{skip} examples skipped (Nonetype)')

serialize(sid_to_seq, f'./data/scop40_208_cluster_sequences_{k}.pkl')
