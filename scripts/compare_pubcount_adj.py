import pandas as pd
import os
import sys
sys.path.insert(0, r'D:\Bioinformatics\d1-benchmark')

labels = ['intogen2024', 'ot_all', 'ot_nolit', 'clingen_label']
nets = ['funmap', 'string_full700', 'string_phys700', 'intact', 'reactome', 'rna_coexp']

orig_runs = 'results/runs/'
adj_runs = 'results/runs/pubcount_adj/'

rows = []
for lbl in labels:
    for net in nets:
        orig_paths = [
            f'{orig_runs}end_to_end_{lbl}_{net}.tsv',
            f'{orig_runs}end_to_end_{net}_{lbl}.tsv',
            f'{orig_runs}orig_backfill/{lbl}_{net}.tsv',
        ]
        orig_path = None
        for p in orig_paths:
            if os.path.exists(p):
                orig_path = p
                break
        adj_path = f'{adj_runs}{lbl}_pubcount_adjusted_{net}.tsv'
        if orig_path is None or not os.path.exists(adj_path):
            print(f'  missing for {lbl} {net}: orig={orig_path}, adj={adj_path}')
            continue
        from d1.engine.evaluate import margins
        orig = pd.read_csv(orig_path, sep='\t')
        adj = pd.read_csv(adj_path, sep='\t')
        orig_m = margins(orig)
        adj_m = margins(adj)
        orig_d = orig_m[orig_m.alpha == 0.5].delta_auroc.mean()
        adj_d = adj_m[adj_m.alpha == 0.5].delta_auroc.mean()
        orig_sd = orig_m[orig_m.alpha == 0.5].delta_auroc.std(ddof=1)
        adj_sd = adj_m[adj_m.alpha == 0.5].delta_auroc.std(ddof=1)
        rows.append({
            'label': lbl, 'network': net,
            'orig_dauroc': orig_d, 'orig_sd': orig_sd,
            'adj_dauroc': adj_d, 'adj_sd': adj_sd,
            'diff': adj_d - orig_d,
            'orig_n_folds': len(orig_m[orig_m.alpha == 0.5]),
        })

df = pd.DataFrame(rows)
print(df.to_string(index=False))
print()
mean_diff = df['diff'].mean()
n_neg = (df['diff'] < 0).sum()
n_total = len(df)
print('Mean diff: ' + f'{mean_diff:.4f}')
print(f'Cells where adjusted < original: {n_neg} / {n_total}')
df.to_csv('results/tables/pubcount_adj_comparison.tsv', sep='\t', index=False, float_format='%.6f')
print('Wrote results/tables/pubcount_adj_comparison.tsv')