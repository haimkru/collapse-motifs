import argparse
import datetime
import json
import os
import time
import tqdm

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset

from collapse.utils import serialize, deserialize
from collapse.models import MLP
from scipy.stats import spearmanr


class TFIDFDataset(Dataset):
    
    def __init__(self, dataset, label_map=None):
        self.data = dataset['data']
        self.labels = dataset['labels']
        self.ids = dataset['ids']
        self.label_map = label_map
    
    def __getitem__(self, i):
        if self.label_map:
            lab = self.label_map.get(self.labels[i])
        else:
            lab = self.labels[i]
        target, decoy = eval(self.ids[i])
        return self.data[i,:], lab, target, decoy
    
    def __len__(self):
        return len(self.data)


def compute_correlations(results):
    per_target = []
    for key, val in results.groupby(['target']):
        # Ignore target with 2 decoys only since the correlations are
        # not really meaningful.
        if val.shape[0] < 3:
            continue
        true = val['true'].astype(float)
        pred = val['pred'].astype(float)
        pearson = true.corr(pred, method='pearson')
        kendall = true.corr(pred, method='kendall')
        spearman = true.corr(pred, method='spearman')
        per_target.append((key, pearson, kendall, spearman))
    per_target = pd.DataFrame(
        data=per_target,
        columns=['target', 'pearson', 'kendall', 'spearman'])

    res = {}
    all_true = results['true'].astype(float)
    all_pred = results['pred'].astype(float)
    res['all_pearson'] = all_true.corr(all_pred, method='pearson')
    res['all_kendall'] = all_true.corr(all_pred, method='kendall')
    res['all_spearman'] = all_true.corr(all_pred, method='spearman')

    res['per_target_pearson'] = per_target['pearson'].mean()
    res['per_target_kendall'] = per_target['kendall'].mean()
    res['per_target_spearman'] = per_target['spearman'].mean()

    # print(
    #     '\nCorrelations (Pearson, Kendall, Spearman)\n'
    #     '    per-target: ({:.3f}, {:.3f}, {:.3f})\n'
    #     '    global    : ({:.3f}, {:.3f}, {:.3f})'.format(
    #     float(res["per_target_pearson"]),
    #     float(res["per_target_kendall"]),
    #     float(res["per_target_spearman"]),
    #     float(res["all_pearson"]),
    #     float(res["all_kendall"]),
    #     float(res["all_spearman"])))
    return res


def train_loop(model, loader, optimizer, device):
    model.train()

    loss_all = 0
    total = 0
    progress_format = 'train loss: {:6.6f}'
    with tqdm.tqdm(total=len(loader), desc=progress_format.format(0)) as t:
        for i, (x, y, _, _) in enumerate(loader):
            feature = x.to(device).to(torch.float32)
            label = y.to(device).to(torch.float32)
            # zero the parameter gradients
            optimizer.zero_grad()
            # forward + backward + optimize
            output = model(feature).squeeze(1)
            loss = F.mse_loss(output, label)
            loss.backward()
            loss_all += loss.item() * len(label)
            total += len(label)
            optimizer.step()
            # stats
            t.set_description(progress_format.format(np.sqrt(loss_all/total)))
            t.update(1)

    return np.sqrt(loss_all / total)


@torch.no_grad()
def test(model, loader, device):
    model.eval()

    losses = []

    targets = []
    decoys = []
    y_true = []
    y_pred = []

    for x, y, target, decoy in loader:
        feature = x.to(device).to(torch.float32)
        label = y.to(device).to(torch.float32)
        output = model(feature).squeeze()
        batch_losses = F.mse_loss(output, label, reduction='none')
        losses.extend(batch_losses.tolist())
        targets.extend(list(target))
        decoys.extend(list(decoy))
        y_true.extend(label.tolist())
        y_pred.extend(output.tolist())
        # print(label.shape, output.shape, id[0], id[1])
    # print([targets, decoys, y_true, y_pred])

    results_df = pd.DataFrame(
        np.array([targets, decoys, y_true, y_pred]).T,
        columns=['target', 'decoy', 'true', 'pred'],
        )

    corrs = compute_correlations(results_df)
    return np.sqrt(np.mean(losses)), corrs, results_df


def save_weights(model, weight_dir):
    torch.save(model.state_dict(), weight_dir)


def train(args, device, test_mode=False):
    # print("Training model with config:")
    # print(str(json.dumps(args.__dict__, indent=4)) + "\n")
    model_str = f'{args.learning_rate}_{args.num_hidden_layers}_{args.hidden_dim}_{args.dropout}_{args.weight_decay}'

    # Save config
    # with open(os.path.join(args.output_dir, 'config.json'), 'w') as f:
    #     json.dump(args.__dict__, f, indent=4)

    np.random.seed(args.random_seed)
    torch.manual_seed(args.random_seed)

    train_dataset = TFIDFDataset(deserialize(os.path.join(args.data_dir, f'psr_tfidf_train.pkl')))
    val_dataset = TFIDFDataset(deserialize(os.path.join(args.data_dir, f'psr_tfidf_val.pkl')))
    test_dataset = TFIDFDataset(deserialize(os.path.join(args.data_dir, f'psr_tfidf_test.pkl')))

    train_loader = DataLoader(train_dataset, args.batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, args.batch_size, shuffle=False)
    test_loader = DataLoader(test_dataset, args.batch_size, shuffle=False)

    for x,y,_,_ in train_loader:
        in_dim = x.shape[1]
        break

    model = MLP(in_dim, [args.hidden_dim]*args.num_hidden_layers, 1, dropout=args.dropout)
    # print(model)
    model.to(device)

    best_val_loss = np.Inf
    best_corrs = None

    optimizer = torch.optim.Adam(model.parameters(), lr=args.learning_rate, weight_decay=args.weight_decay)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=5)

    patience = 0
    for epoch in range(1, args.num_epochs+1):
        start = time.time()
        train_loss = train_loop(model, train_loader, optimizer, device)
        val_loss, corrs, val_df = test(model, val_loader, device)
        scheduler.step(val_loss)
        if val_loss < best_val_loss:
            print(f"\nSave model at epoch {epoch:03d}, val_loss: {val_loss:.4f}")
            save_weights(model, os.path.join(args.output_dir, model_str + '_best_weights.pt'))
            best_val_loss = val_loss
            best_corrs = corrs
            patience = 0
        else:
            patience += 1
        elapsed = (time.time() - start)
        # print('Epoch {:03d} finished in : {:.3f} s'.format(epoch, elapsed))
        print('\tTrain RMSE: {:.7f}, Val RMSE: {:.7f}, Per-target Spearman R: {:.7f}, Global Spearman R: {:.7f}'.format(
            train_loss, val_loss, corrs['per_target_spearman'], corrs['all_spearman']))
        if patience >= 10:
            print('Val loss did not improve for 10 epochs, stopping training')
            break

    if test_mode:
        
        model.load_state_dict(torch.load(os.path.join(args.output_dir, model_str + '_best_weights.pt')))
        rmse, corrs, test_df = test(model, test_loader, device)
        test_df.to_pickle(os.path.join(args.output_dir, model_str + '_test_results.pkl'))
        print('Test RMSE: {:.7f}, Per-target Spearman R: {:.7f}, Global Spearman R: {:.7f}'.format(
            rmse, corrs['per_target_spearman'], corrs['all_spearman']))
        test_file = os.path.join(args.output_dir, 'test_results.txt')
        with open(test_file, 'a+') as out:
            out.write('{}\t{:.7f}\t{:.7f}\t{:.7f}\n'.format(
                model_str, rmse, corrs['per_target_pearson'], corrs['all_pearson'], corrs['per_target_kendall'], corrs['all_kendall'], corrs['per_target_spearman'], corrs['all_spearman']))

    return best_val_loss, best_corrs['per_target_spearman'], best_corrs['all_spearman']


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str, default='./data')
    parser.add_argument('--output_dir', type=str, default='./data/psr_out')
    parser.add_argument('--mode', type=str, default='train',
                        choices=['train', 'test'])
    parser.add_argument('--learning_rate', type=float, default=1e-5)
    parser.add_argument('--num_epochs', type=int, default=100)
    parser.add_argument('--batch_size', type=int, default=64)
    parser.add_argument('--hidden_dim', type=int, default=2048)
    parser.add_argument('--num_hidden_layers', type=int, default=2)
    parser.add_argument('--dropout', type=float, default=0.5)
    parser.add_argument('--weight_decay', type=float, default=0.01)
    parser.add_argument('--random_seed', type=int, default=np.random.randint(0, 1000000))

    args = parser.parse_args()

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Set up output dir
    args.output_dir = os.path.join(args.output_dir, 'psr')
    assert args.output_dir != None
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Running mode {args.mode:} with seed {args.random_seed:} "
          f"and output dir {args.output_dir}")
    train(args, device, args.mode=='test')

# Best val results
# lr 1e-5, 2 layers, 2048 hidden dim, dropout 0.75
#  Val RMSE: 0.1954315, Per-target Spearman R: 0.5179840, Global Spearman R: 0.7553714