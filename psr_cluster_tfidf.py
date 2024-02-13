import sys
import numpy as np
import pandas as pd
import collections as col
import pickle
from tqdm import tqdm
from atom3d.datasets import load_dataset
from collapse.utils import serialize, deserialize

k = 50000
split = sys.argv[1]

pdb_to_scop = {}
with open('../data/dir.des.scope.2.08-stable.txt') as f:
    for line in f:
        if line[0] == '#':
            continue
        else:
            l = line.strip().split('\t')
            _, cat, cls, pdbcd, desc = l
            if cat == 'px':
                pdb_to_scop[pdbcd] = cls


sid_to_seq = deserialize(f'./data/psr_cluster_sequences_{split}.pkl')
sid_to_label = deserialize(f'./data/psr_cluster_labels_{split}.pkl')

ids = []
labels = []
idf = np.ones(k)
tfidf_mat = np.zeros((len(sid_to_seq), k))
for i, (sid, seq) in enumerate(tqdm(sid_to_seq.items())):
    np.add.at(tfidf_mat[i,:], seq, 1)
    idf[seq] += 1
    labels.append(sid_to_label[sid])
    ids.append(sid)
    
tf = tfidf_mat / tfidf_mat.sum(1)[:,np.newaxis]

idf = np.log(float(len(sid_to_seq)) / idf)

tf_idf = tf * idf[np.newaxis, :]

#remove zero-valued features
# tf_idf = tf_idf[:,np.any(tf_idf != 0, axis=0)]

serialize({'data': tf_idf, 'labels': labels, 'ids': ids}, f'data/psr_tfidf_{split}.pkl')
    