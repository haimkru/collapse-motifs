import sys
import numpy as np
import pandas as pd
import pickle
import collections as col
from tqdm import tqdm
from atom3d.datasets import load_dataset, make_lmdb_dataset
from collapse.utils import serialize, deserialize
from collapse.atom_info import aa_abbr

split = sys.argv[1]
k = 50000

data_path = f'/scratch/users/aderry/psr_embedded/{split}'

clusterer = deserialize(f'data/pdb100_cluster_fit_{k}.pkl')

dataset = load_dataset(data_path, 'lmdb')

sid_to_seq = {}
sid_to_label = {}
skip = 0
for item in tqdm(dataset):
    if item is None:
        skip += 1
        continue

    sid = item['id']
    emb = item['embeddings']
    clusters = clusterer.predict(emb)
    sid_to_seq[sid] = clusters
    sid_to_label[sid] = item['scores']['gdt_ts']
        
print(f'{skip} examples skipped (Nonetype)')

serialize(sid_to_seq, f'./data/psr_cluster_sequences_{split}.pkl')
serialize(sid_to_label, f'./data/psr_cluster_labels_{split}.pkl')