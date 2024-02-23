import numpy as np
import pandas as pd
import os
import argparse
import torch
from collapse.data import EmbedTransform
from atom3d.datasets import load_dataset, make_lmdb_dataset
import atom3d.util.file as fi
from collapse import initialize_model

parser = argparse.ArgumentParser()
parser.add_argument('data_dir', type=str)
parser.add_argument('out_dir', type=str)
parser.add_argument('--split_id', type=int, default=0)
parser.add_argument('--filetype', type=str, default='pdb')
parser.add_argument('--num_splits', type=int, default=1)
args = parser.parse_args()

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

all_files = fi.find_files(args.data_dir, args.filetype)
# names = [str(x).split('/')[-1].split('-')[1] for x in all_files]

model = initialize_model(device=device)
transform = EmbedTransform(model, device=device, include_hets=False, compute_res_graph=False)
dataset = load_dataset(args.data_dir, args.filetype, transform=transform)

if args.num_splits > 1:
    out_path = os.path.join(args.out_dir, f'tmp_{args.split_id}')
    os.makedirs(out_path, exist_ok=True)

    split_idx = np.array_split(np.arange(len(all_files)), args.num_splits)[args.split_id - 1]
    print(f'Processing split {args.split_id} with {len(split_idx)} examples...')
    dataset = torch.utils.data.Subset(dataset, split_idx)
else:
    os.makedirs(args.out_dir, exist_ok=True)
    out_path = args.out_dir

make_lmdb_dataset(dataset, out_path, serialization_format='pkl', filter_fn=lambda x: (x is None))
