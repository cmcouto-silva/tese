## Import required libraries
import os
import numpy as np
import pandas as pd
import gseapy as gp
from more_itertools import consecutive_groups

# Set working directory
os.chdir('dataset_S1/enrichr')
# Create directory for outputs
os.makedirs('results', exist_ok=True)

# Function to split genes and associated info
def separate_genes(df, gene_col = 'GENE', keep_empty = False):
    if not isinstance(gene_col, list):
        gene_col = [gene_col]
    if not keep_empty:
        df = df[(df.GENE != '') | (df.GENE.isna())].reset_index(drop=True)
    df = (df.assign(**{col: df[col].str.split(',') for col in gene_col})
          .apply(pd.Series.explode)
          .reset_index(drop=True))
    return df

# Set common function for getting outliers
def get_top_genes(df, cutoff, mode = "persnp", pbsw_correction = True):

    df = df.copy()
    
    grp = ['CHR', 'POS', 'SNP']
    p_cols = ['PVALUE', 'LOG_PVALUE']
    gene_info_cols = ['GENE', 'GENE_ID', 'FUNCTION']
    
    for method in ['IHS', 'XPEHH', 'PBS']:
        if method in df.columns:
            df.rename(columns={method:'SCORE'}, inplace=True)
    
    if mode == 'persnp':
        df = df[df.SCORE >= df.SCORE.quantile(cutoff)]
        if pbsw_correction and cutoff > 0 and 'W_ID' in df.columns:
            for i, group in enumerate(consecutive_groups(df.W_ID)):
                df.loc[df.W_ID.isin(list(group)), 'W_GROUP'] = f'WG_{i}'
            df = (df.sort_values('SCORE', ascending=False)
                  .drop_duplicates('W_GROUP')
                  .drop(columns=['W_ID', 'W_GROUP'])
                  .reset_index(drop=True)
                  )
        
        df = df[grp + gene_info_cols + ['SCORE'] + p_cols]
        df = separate_genes(df, gene_info_cols)
        
        df = (df.loc[df.groupby('GENE')['SCORE'].idxmax()]
              .sort_values(['CHR','POS','GENE'])
              .reset_index(drop=True)
              .groupby(list(set(df.columns) - set(gene_info_cols)))
              .agg(lambda x: ','.join(x))
              .sort_values(['CHR', 'POS', 'GENE'])
              .reset_index()
              .loc[:, grp + gene_info_cols + ["SCORE"] + p_cols]
              )
    else:
        grp = grp[:-1]
        gene_info_cols = gene_info_cols[:-1]
        
        df = separate_genes(df, gene_info_cols)
        
        dup_genes = df.groupby('GENE')['CHR'].unique().apply(len)
        dup_genes = dup_genes[dup_genes>1].index
        
        df = df[~df.GENE.isin(dup_genes)]
        df = df.groupby(gene_info_cols, as_index=False).agg({
            'CHR': np.unique,
            'POS': lambda x: np.round(np.mean(x)),
            'SCORE': np.mean
            })
        
        df = (df.assign(POS = df.POS.astype(int))
              .assign(SNP = lambda x: x[['CHR','POS']].astype(str).agg(':'.join, axis=1))
              )
        
        df = (df[df.SCORE >= df.SCORE.quantile(cutoff)]
              .reset_index(drop=True)
              .groupby(list(set(df.columns) - set(gene_info_cols)))
              .agg(lambda x: ','.join(x))
              .reset_index()
              )
        
        df = (df.sort_values(['CHR', 'POS', 'GENE'])
              .loc[:, grp + gene_info_cols + ['SCORE']]
              .reset_index(drop=True)
              )
        
    genes = sorted(set(df.GENE.str.split(',').explode()))
    return genes
    
# Set Arguments
modes = ['persnp', 'pergene']
methods = ['pbs_windowed', 'xpehh', 'ihs']
cutoff_vals = [.999, .995, .99, .95]

# Load Data
file_path='../results_bin/{}/{}.feather'
results = {mode: {
    method: pd.read_feather(file_path.format(mode, method)) for method in methods
    } for mode in modes }

## Enrich Analysis
for mode in modes:
    print(f'  {mode}  '.center(25, '-'))
    print('-' * 25, end='\n\n')

    for method in methods:
        print(f'Method: \033[1m{method.upper()}\033[m')
        print('-'*25, end='\n')
    
        for cutoff in cutoff_vals:
            print(f'- pvalue <= {cutoff} ...')
            enrich_obj = gp.enrichr(
                gene_list=get_top_genes(results[mode][method], cutoff, mode),
                description=f'{method}_{mode}_pvalue{cutoff}',
                gene_sets=['GWAS_Catalog_2019', 'KEGG_2019_Human'],
                organism='Human',
                outdir=None,
                format='png'
            )
            enrich_res = enrich_obj.results.sort_values('Adjusted P-value').reset_index(drop=True)
            enrich_res.to_csv(f'results/{enrich_obj.descriptions}.txt', sep='\t', index=0, header=1)
        print()

