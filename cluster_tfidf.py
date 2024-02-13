import sys
import numpy as np
import pandas as pd
import collections as col
import pickle
from tqdm import tqdm
from atom3d.datasets import load_dataset
from collapse.utils import serialize, deserialize

k=int(sys.argv[1])

pdb_to_scop = {}
with open('../mappings/dir.des.scope.2.08-stable.txt') as f:
    for line in f:
        if line[0] == '#':
            continue
        else:
            l = line.strip().split('\t')
            _, cat, cls, pdbcd, desc = l
            if cat == 'px':
                pdb_to_scop[pdbcd] = cls

dataset = load_dataset('/scratch/users/aderry/scop40/scop40_208_lmdb/full', 'lmdb')
with open(f'data/scop40_208_cluster_sequences_{k}.pkl', 'rb') as f:
    sid_to_seq = pickle.load(f)

ids = []
labels = []
idf = np.ones(k)
tfidf_mat = np.zeros((len(sid_to_seq), k))
for i, (sid, seq) in enumerate(tqdm(sid_to_seq.items())):
    np.add.at(tfidf_mat[i,:], seq, 1)
    idf[seq] += 1
    labels.append(pdb_to_scop.get(sid))
    ids.append(sid)
    
tf = tfidf_mat / tfidf_mat.sum(1)[:,np.newaxis]

idf = np.log(float(len(sid_to_seq)) / idf)

tf_idf = tf * idf[np.newaxis, :]


with open(f'data/scop40_208_tfidf_{k}.pkl', 'wb') as f:
    pickle.dump({'data': tf_idf, 'labels': labels, 'ids': ids}, f, protocol=4)