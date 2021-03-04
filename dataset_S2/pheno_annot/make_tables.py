import os
import pickle
import numpy as np
import pandas as pd

os.chdir('dataset_S2/pheno_annot')

os.makedirs('results', exist_ok=True)

results = {}

for result in ['results_annot', 'results_annot_contrast']:
    with open(f'{result}.pickle', 'rb') as file:
        results[result] = pickle.load(file)

databases = ['gwas', 'ensembl']#, 'genecard']
modes = ['persnp', 'pergene']
methods = ['pbs_windowed', 'xpehh', 'ihs']

for result in results:
    for database in databases:
        for mode in modes:
            for method in methods:
                if result == 'results_annot':
                    results[result][database][mode][method].to_csv(f'results/{database}_{mode}_{method}.csv', index=False)
                if result == 'results_annot_reverse' and method != 'ihs':
                    results[result][database][mode][method].to_csv(f'results/{database}_{mode}_{method}_contrast.csv', index=False)

