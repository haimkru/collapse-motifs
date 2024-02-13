import os
import numpy as np
import pandas as pd
import pickle
from sklearn.cluster import MiniBatchKMeans, BisectingKMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
from collapse.utils import serialize, deserialize
import matplotlib.pyplot as plt
import seaborn as sns
import collections as col
import time

SEED = 77

def calc_entropy(df):
    counts = df.value_counts(['cluster', 'scop']).reset_index()
    counts['prob'] = counts['count'] / counts.groupby('cluster')['count'].transform('sum')
    counts['entropy'] = -counts.prob*np.log(counts.prob)
    cluster_entropy = counts.groupby('cluster')['entropy'].sum()
    avg_scop_count = counts.groupby('cluster')['scop'].count().mean()
    avg_samples = df.groupby('cluster').size().mean()
    return avg_scop_count, avg_samples, (counts.groupby('cluster')['count'].sum() / counts['count'].sum() * cluster_entropy).sum()

def parse_scop_file(fname):
    scop_names = {}
    pdb_to_scop = col.defaultdict(dict)
    with open(fname) as f:
        for line in f:
            if line[0] == '#':
                continue
            else:
                l = line.strip().split('\t')
                _, cat, cls, pdbcd, desc = l
                if cat in ['cl', 'cf', 'sf', 'fa']:
                    scop_names[cls] = f'{cls} ({desc})'
                elif cat == 'px':
                    pdbc = pdbcd[1:6]
                    loc_range = desc.split()[-1].split(':')[-1]
                    if len(loc_range) == 0:
                        pdb_to_scop[pdbc]['all'] = cls
                    else:
                        loc_split = loc_range.split('-')
                        try:
                            if len(loc_split) == 2:
                                start, end = [int(x) for x in loc_split]
                            elif len(loc_split) == 3:
                                start, end = int(''.join(loc_split[:2])), int(loc_split[-1])
                            elif len(loc_split) == 4:
                                start, end = int(''.join(loc_split[:2])), int(''.join(loc_split[2:]))
                        except ValueError:
                            continue
                        for i in range(start, end+1):
                            pdb_to_scop[pdbc][str(i)] = cls
    return pdb_to_scop

def get_scop(pdbc, resid, pdb_to_scop):
    pdbc = pdbc.split('_')[0] + pdbc.split('_')[1].lower()
    scop_dict = pdb_to_scop[pdbc]
    if 'all' in scop_dict:
        return scop_dict['all']
    resnum = resid[1:]
    return scop_dict.get(resnum, 'N/A')

if __name__=="__main__":

    with open('../COLLAPSE/data/datasets/pdb100_embeddings/pdb_embeddings.pkl', 'rb') as f:
        data = pickle.load(f)
        
    np.random.seed(SEED)
    # rand = np.random.choice(np.arange(data['embeddings'].shape[0]), size=1000000)
    # rand = np.arange(data['embeddings'].shape[0])
    emb_sample = data['embeddings']#[rand, :]
    # emb_sample = normalize(data['embeddings'])
    pdbs_sample = np.array(data['pdbs'])#[rand]
    pdb_chains_sample = [x[1:6] for x in pdbs_sample]
    # chains_sample = np.array(data['chains'])#[rand]
    resids_sample = np.array(data['resids'])#[rand]
    labels_sample = list(zip(pdbs_sample, resids_sample))
    
    scop_file = './data/dir.des.scope.2.08-stable.txt'

    pdb_to_scop = parse_scop_file(scop_file)
    
    sumsquares = []
    entropies = []
    folds = []
    fams = []
    sfams = []
    counts = []
    eps_list = [0.01, 0.05, 0.1, None]
    k_list = [20, 100, 1000, 5000, 10000, 20000, 30000, 40000, 50000, 100000]
    for k in k_list:
        start = time.time()
        print(f'Running for k={k}...')
        clusterer = MiniBatchKMeans(n_clusters=k, n_init=3, verbose=0, batch_size=4096, random_state=SEED)
        cluster_fit = clusterer.fit(emb_sample)
        print(f'clustering finished in {(time.time() - start)/60} minutes')
        sumsquares.append(cluster_fit.inertia_)
        cluster_idx = cluster_fit.labels_
        df = pd.DataFrame({'cluster': cluster_idx.astype(str), 'scop': [get_scop(p,r) for p,r in zip(pdbs_sample, resids_sample)]})
        df['scop_sf'] = ['.'.join(x.split('.')[:3]) if x is not None else None for x in df['scop']]
        df['scop_fold'] = ['.'.join(x.split('.')[:2]) if x is not None else None for x in df['scop']]
        avg_fams, avg_samples, entropy = calc_entropy(df)
        avg_sfs = df.value_counts(['cluster', 'scop_sf']).reset_index().groupby('cluster')['scop_sf'].count().mean()
        avg_folds = df.value_counts(['cluster', 'scop_fold']).reset_index().groupby('cluster')['scop_fold'].count().mean()
        print('WCSS', cluster_fit.inertia_, '| entropy', entropy, '| families/cluster', avg_fams, '| superfams/cluster', avg_sfs, '| folds/cluster', avg_folds, '| samples/cluster', avg_samples)

        entropies.append(entropy)
        folds.append(avg_folds)
        sfams.append(avg_sfs)
        fams.append(avg_fams)
        counts.append(avg_samples)
        
        serialize(cluster_fit, f'./data/pdb100_cluster_fit_{k}.pkl')
    
    res_df = pd.DataFrame({'k': k_list, 'Within-cluster sum of squares': sumsquares, 'SCOP entropy': entropies, 'Families per cluster': fams, 'Superfamilies per cluster': sfams, 'Folds per cluster': folds, 'Samples per cluster': counts})

    res_df.to_csv('./data/pdb100_cluster_tuning.csv', index=False)
    

