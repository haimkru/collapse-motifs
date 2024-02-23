import sys
import numpy as np
import pandas as pd
import collections as col
import pickle
from tqdm import tqdm
from atom3d.datasets import load_dataset
from collapse.utils import serialize, deserialize

## For SCOP labeling
# pdb_to_scop = {}
# with open('../mappings/dir.des.scope.2.08-stable.txt') as f:
#     for line in f:
#         if line[0] == '#':
#             continue
#         else:
#             l = line.strip().split('\t')
#             _, cat, cls, pdbcd, desc = l
#             if cat == 'px':
#                 pdb_to_scop[pdbcd] = cls

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('cluster_sequences', type=str, help='cluster-featurized dataset produced by cluster_featurize.py (pkl format)')
    parser.add_argument('id_to_label', type=str, help='mapping from id to labels (pkl format)')
    parser.add_argument('out_path', type=str, help='output file path (.pkl format)')
    parser.add_argument('--k', type=int, default=50000)
    args = parser.parse_args()

    sid_to_seq =  deserialize(args.cluster_sequences)
    label_dict = deserialize(args.id_to_label)
    
    ids = []
    labels = []
    idf = np.ones(k)
    tfidf_mat = np.zeros((len(sid_to_seq), k))
    for i, (sid, seq) in enumerate(tqdm(sid_to_seq.items())):
        np.add.at(tfidf_mat[i,:], seq, 1)
        idf[seq] += 1
        labels.append(label_dict.get(sid))
        ids.append(sid)
        
    tf = tfidf_mat / tfidf_mat.sum(1)[:,np.newaxis]
    
    idf = np.log(float(len(sid_to_seq)) / idf)
    
    tf_idf = tf * idf[np.newaxis, :]
    
    serialize({'data': tf_idf, 'labels': labels, 'ids': ids}, args.out_path)