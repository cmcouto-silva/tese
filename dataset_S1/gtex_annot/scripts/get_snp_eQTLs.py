# Load libraries
import os
import sys
import json
import h5py
import allel
import pickle
import numpy as np
import pandas as pd
from time import sleep
from collections import defaultdict

# Set working directory
os.chdir('dataset_S1/gtex_annot')
sys.path.insert(0, os.path.join(os.getcwd(), 'scripts'))

# Load customized functions
from functions import get_top, separate_genes, get_eQTL, get_eqtl_data

# Load Results
methods = ['pbs_windowed', 'xpehh']
ds1_persnp = {method: pd.read_feather(f'../results_bin/persnp/{method}.feather') for method in methods}
ds1_persnp_top = {method: get_top(ds1_persnp[method], .999) for method in methods}

# Load Phased Data
vcf = h5py.File('../data/phased_adj_alleles/dataset_S1.h5', 'r')
pop = pd.read_csv('../data/populations/populations.csv')

# Manage Populations
samples = pd.DataFrame(list(map(lambda x: x.decode().split(':'), vcf['samples'])), columns=['FID','IID'])
if not (samples == pop[['FID','IID']]).all().all():
    print('Samples not in same order!')
pop_selection = pop.SelGroup.isin(['EastAsia', 'Mesoamerica', 'Amazonia'])

# Subset genotypes
gt = allel.GenotypeChunkedArray(vcf['calldata/GT'])
gt = gt.subset(sel1=pop_selection) 

# Subset populations
pop = pop[pop_selection].reset_index(drop=1)
samples = samples[pop_selection].reset_index(drop=1)

# Set target groups
groups = {
    'all': list(range(len(samples))),
    'amz': pop[pop.SelGroup == 'Amazonia'].index.tolist(),
    'mes': pop[pop.SelGroup == 'Mesoamerica'].index.tolist(),
    'eas': pop[pop.SelGroup == 'EastAsia'].index.tolist()
}

# Calculate Allele Frequency
groups_afreq = gt.count_alleles_subpops(groups, max_allele=1)
amz_ac = groups_afreq['amz'][:] / np.unique(groups_afreq['amz'][:].sum(axis=1))
mes_ac = groups_afreq['mes'][:] / np.unique(groups_afreq['mes'][:].sum(axis=1))
eas_ac = groups_afreq['eas'][:] / np.unique(groups_afreq['eas'][:].sum(axis=1))

# Make DataFrame with allele frequencies + SNP info
decode_array = np.vectorize(lambda s: s.decode())

df_freq = pd.DataFrame({
    'CHR': decode_array(vcf['variants/CHROM'][:]),
    'SNP': decode_array(vcf['variants/ID'][:]),
    'REF': decode_array(vcf['variants/REF'][:]),
    'ALT': decode_array(vcf['variants/ALT'][:,0]),
    'AMZ_REF': amz_ac[:,0], 'AMZ_ALT': amz_ac[:,1],
    'MES_REF': mes_ac[:,0], 'MES_ALT': mes_ac[:,1],
    'EAS_REF': eas_ac[:,0], 'EAS_ALT': eas_ac[:,1]
    })

# Create input for GTEx API and customized functions
ds1_persnp_eqtl_df = dict()
ds1_persnp_eqtl_dict = dict()

for method in methods:
    print('-'*50, end='\n')
    print(f'  {method.upper()}  '.center(50, '-'))
    print('-'*50, end='\n\n')
    
    eqtl_input_df = df_freq[df_freq.SNP.isin(ds1_persnp_top[method].SNP)].copy()
    eqtl_input_df.insert(4, 'TARGET', np.where(
        eqtl_input_df.AMZ_REF > eqtl_input_df.MES_REF,
        eqtl_input_df.REF, eqtl_input_df.ALT
    ))

    col_order = [
        'CHR', 'SNP', 'GENE',
        'REF', 'ALT', 'TARGET', 
        'AMZ_REF', 'AMZ_ALT',
        'MES_REF', 'MES_ALT',
        'EAS_REF', 'EAS_ALT'
    ]

    eqtl_input_df = eqtl_input_df.merge(ds1_persnp_top[method][['SNP','GENE']])
    eqtl_input_df = eqtl_input_df.loc[:, col_order]
    eqtl_input_df = separate_genes(eqtl_input_df)
    
    # Get eQTLs
    eqtl_dict = {}
    
    for idx,row in eqtl_input_df.iterrows():
        snpgene = row['SNP'] + '-' + row['GENE']
        print(f'  {snpgene} {idx+1}/{eqtl_input_df.shape[0]}  '.center(50, '-')); print('-'*50, end='\n\n')
        eqtl_dict[snpgene] = get_eQTL(row['SNP'], row['GENE'], row['TARGET'])
        print()
    
    snp_dict = defaultdict(dict)
    for snp in eqtl_dict:
        tissue_dict = defaultdict(dict)
        for tissue in eqtl_dict[snp]:
            df = get_eqtl_data(eqtl_dict[snp][tissue])
            df_gt = df.groupby('Genotype', as_index=False).mean()
            target = eqtl_dict[snp][tissue]['target']*2
            if target in df_gt.Genotype.tolist():
                if np.all(df_gt[df_gt.Genotype==target]['Data'].values > df_gt[df_gt.Genotype!=target]['Data'].values):
                    tissue_dict['up'][tissue] = df_gt
                elif np.all(df_gt[df_gt.Genotype==target]['Data'].values < df_gt[df_gt.Genotype!=target]['Data'].values):
                    tissue_dict['down'][tissue] = df_gt
            else:
                target_idx = np.array([target[0] in gt for gt in df_gt.Genotype])
                if df_gt[target_idx]['Data'].values > df_gt[~target_idx]['Data'].values:
                    tissue_dict['up'][tissue] = df_gt
                elif df_gt[target_idx]['Data'].values < df_gt[~target_idx]['Data'].values:
                    tissue_dict['down'][tissue] = df_gt
        snp_dict[snp] = tissue_dict
  
    eqtl_df = pd.DataFrame({
        'SNP': [k.split('-')[0] for k in eqtl_dict],
        'GENE': [k.split('-')[1] for k in eqtl_dict],
        'UP': [list(df['up'].keys()) for df in snp_dict.values()],
        'DOWN': [list(df['down'].keys()) for df in snp_dict.values()]
    })

    col_order_df = eqtl_input_df.columns.tolist() + ['UP', 'DOWN']
    eqtl_df = eqtl_df.merge(eqtl_input_df).loc[:, col_order_df]
    
    ds1_persnp_eqtl_df[method] = eqtl_df
    ds1_persnp_eqtl_dict[method] = eqtl_dict
    
    if method == 'pbs_windowed':
        method = 'pbsw'
    
    eqtl_df[['UP','DOWN']] = eqtl_df[['UP','DOWN']].agg(lambda x: x.str.join(','))
    eqtl_df.to_csv(f'eqtl_persnp_{method}.csv', index=False)
    
    with open(f'eqtl_persnp_{method}.json', 'w') as file:
        json.dump(eqtl_dict, file, indent=4)
    
# Write files
with open('eqtl_persnp_df.pickle', 'wb') as file:
    pickle.dump(ds1_persnp_eqtl_df, file)

with open('eqtl_persnp_dict.pickle', 'wb') as file:
    pickle.dump(ds1_persnp_eqtl_dict, file)

print('done')
