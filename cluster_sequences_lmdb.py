import sys
import os
import numpy as np
import pandas as pd
import torch
import pickle
import collections as col
from tqdm import tqdm
import argparse
from atom3d.datasets import load_dataset, make_lmdb_dataset
from collapse.utils import serialize, deserialize


parser = argparse.ArgumentParser()
parser.add_argument('k', type=int)
parser.add_argument('--split_id', type=int, default=0)
parser.add_argument('--num_splits', type=int, default=1)
args = parser.parse_args()

data_path = '/scratch/users/aderry/pdb_lmdb_full/full'

clusterer = deserialize(f'data/pdb100_cluster_fit_{args.k}.pkl')

def cluster_transform(item):
    if item is None:
        return None
    clusters = clusterer.predict(item['embeddings'])
    item['id'] = item['id'].split('.')[0][3:]
    item['clusters'] = clusters
    return item

dataset = load_dataset(data_path, 'lmdb', transform=cluster_transform)
out_dir = f'/scratch/users/aderry/pdb_clusters_{args.k}'

if args.num_splits > 1:
    out_path = os.path.join(out_dir, f'tmp_{args.split_id}')
    os.makedirs(out_path, exist_ok=True)

    split_idx = np.array_split(np.arange(len(dataset)), args.num_splits)[args.split_id - 1]
    print(f'Processing split {args.split_id} with {len(split_idx)} examples...')
    dataset = torch.utils.data.Subset(dataset, split_idx)
else:
    os.makedirs(out_dir, exist_ok=True)
    out_path = out_dir

make_lmdb_dataset(dataset, out_path, serialization_format='pkl', filter_fn=lambda x: (x is None))