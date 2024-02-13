import os
import pickle
import argparse
import numpy as np
import pandas as pd
import collections as col
from fastdist import fastdist
import time
from tqdm import tqdm

def compute_sensitivity(sorted_idx, labels, ids, query_labels, query_ids):
    ct_dict = col.Counter(labels)
    
    sorted_labels = np.apply_along_axis(lambda x: labels[x], 0, sorted_idx)
    sorted_ids = np.apply_along_axis(lambda x: ids[x], 0, sorted_idx)
    
    correct_counts = np.array([ct_dict[l] - 1 for l in query_labels])
    self_hits = (sorted_ids == query_ids[:, np.newaxis])
    sorted_labels = sorted_labels[~self_hits].reshape(len(self_hits), -1)
    sorted_ids = sorted_ids[~self_hits].reshape(len(self_hits), -1)
    # print(correct_counts)
    # print(sorted_labels[:, :20])
    correct = sorted_labels == query_labels[:, np.newaxis]
    first_fp = np.argmin(correct, 1)
    sensitivity = first_fp / correct_counts
    return sensitivity

def main(args):
    with open(args.database, 'rb') as f:
        data = pickle.load(f)
    
    base_fname = args.database.split('/')[-1].split('.')[0]
    
    pdb_to_scop = {}
    with open('../mappings/dir.des.scope.2.01-stable.txt') as f:
        for line in f:
            if line[0] == '#':
                continue
            else:
                l = line.strip().split('\t')
                _, cat, cls, pdbcd, desc = l
                if cat == 'px':
                    pdb_to_scop[pdbcd] = cls
    
    embeddings = data['data']#[1:, :]
    labels = np.array(data['labels'])
    exist_idx = labels != None
    embeddings = embeddings[exist_idx]
    labels = labels[exist_idx]
    ids = np.array(data['ids'])[exist_idx]
    # labels = np.array([pdb_to_scop[i] for i in ids])
    ct_dict = col.Counter(labels)
    label_counts = np.array([ct_dict[l] for l in labels])
    labels_over_ct1 = label_counts > 1
    print(f'Keeping {np.sum(labels_over_ct1)} of {len(labels)} proteins with at least two SCOP family members')
    
    # rand = np.random.choice(np.where(labels_over_ct1)[0], size=args.num)
    emb_query = embeddings[labels_over_ct1]
    lab_query = labels[labels_over_ct1]
    id_query = ids[labels_over_ct1]
    
    start = time.time()
    
    cosines = fastdist.cosine_matrix_to_matrix(emb_query, embeddings)
    print('computed all-vs-all similarity in {:.2f} sec'.format(time.time() - start))
    sorted_idx = np.argsort(-1*cosines, 1)
    
    all_results = []
    for i in range(len(sorted_idx)):
        idx = sorted_idx[i]
        tgts = ids[idx]
        # tgt_scops = labels[idx]
        scores = cosines[i, idx]
        res = pd.DataFrame({'query': id_query[i], 'target': list(tgts), 'score': scores.tolist()})
        all_results.append(res)
    
    all_results = pd.concat(all_results)
    all_results = all_results[all_results.score > 0]
    all_results.to_csv(f'data/{base_fname}_search_results.csv', index=False)
    
    elapsed = time.time() - start
    print(f'Elapsed time: {elapsed} sec')
    
    sens_fam = compute_sensitivity(sorted_idx, labels, ids, lab_query, id_query)
    print(f'Family-level sensitivity: {np.mean(sens_fam)}')
    
    truncate_sf = np.vectorize(lambda x: '.'.join(x.split('.')[:3]))
    labels_sf = truncate_sf(labels)
    query_labels_sf = truncate_sf(lab_query)
    sens_sf = compute_sensitivity(sorted_idx, labels_sf, ids, query_labels_sf, id_query)
    print(f'Superfamily-level sensitivity: {np.mean(sens_sf)}')
    
    truncate_fold = np.vectorize(lambda x: '.'.join(x.split('.')[:2]))
    labels_fold = truncate_fold(labels)
    query_labels_fold = truncate_fold(lab_query)
    sens_fold = compute_sensitivity(sorted_idx, labels_fold, ids, query_labels_fold, id_query)
    print(f'Fold-level sensitivity: {np.mean(sens_fold)}')
    
    results = pd.DataFrame({'query': id_query, 'scop': lab_query, 'family': sens_fam, 'superfamily': sens_sf, 'fold': sens_fold})
    results.to_csv(f'data/{base_fname}_search_sensitivity.csv', index=False)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--database', type=str, default='data/scop40_201_tfidf_50000.pkl')
    parser.add_argument('--num', type=int, default=10)
    args = parser.parse_args()

    main(args)