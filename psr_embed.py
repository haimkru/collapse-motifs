import os
import torch
from atom3d.datasets import load_dataset, make_lmdb_dataset
from collapse.data import EmbedTransform, initialize_model

for split in ['train', 'val', 'test']:
    print(f'Processing {split} split')
    data_dir = f'/scratch/users/raphtown/atom3d_mirror/lmdb/PSR/splits/split-by-year/data/{split}'
    out_path = f'/scratch/users/aderry/psr_embedded/{split}'
    os.makedirs(out_path, exist_ok=True)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    model = initialize_model(device=device)
    transform = EmbedTransform(model, device=device, include_hets=False, compute_res_graph=False)
    dataset = load_dataset(data_dir, 'lmdb', transform=transform)

    make_lmdb_dataset(dataset, out_path, serialization_format='pkl', filter_fn=lambda x: (x is None))