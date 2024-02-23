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

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('data_path', type=str)
    parser.add_argument('out_path', type=str)
    parser.add_argument('--k', type=int, default=50000)
    parser.add_argument('--split_id', type=int, default=0)
    parser.add_argument('--num_splits', type=int, default=1)
    args = parser.parse_args()
    
    clusterer = deserialize(f'./data/pdb100_cluster_fit_{k}.pkl')
    
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
    
    serialize(sid_to_seq, args.out_path)
