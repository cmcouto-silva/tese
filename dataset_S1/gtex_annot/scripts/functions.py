# Import required librarires
import requests
import numpy as np
import pandas as pd
from more_itertools import consecutive_groups

# Available tissues for Gtex_v7
tissues = [
    'Adipose_Subcutaneous',
    'Adipose_Visceral_Omentum',
    'Adrenal_Gland',
    'Artery_Aorta',
    'Artery_Coronary',
    'Artery_Tibial',
    'Brain_Amygdala',
    'Brain_Anterior_cingulate_cortex_BA24',
    'Brain_Caudate_basal_ganglia',
    'Brain_Cerebellar_Hemisphere',
    'Brain_Cerebellum',
    'Brain_Cortex',
    'Brain_Frontal_Cortex_BA9',
    'Brain_Hippocampus',
    'Brain_Hypothalamus',
    'Brain_Nucleus_accumbens_basal_ganglia',
    'Brain_Putamen_basal_ganglia',
    'Brain_Spinal_cord_cervical_c-1',
    'Brain_Substantia_nigra',
    'Breast_Mammary_Tissue',
    'Cells_EBV-transformed_lymphocytes',
    'Cells_Transformed_fibroblasts',
    'Colon_Sigmoid',
    'Colon_Transverse',
    'Esophagus_Gastroesophageal_Junction',
    'Esophagus_Mucosa',
    'Esophagus_Muscularis',
    'Heart_Atrial_Appendage',
    'Heart_Left_Ventricle',
    'Liver',
    'Lung',
    'Minor_Salivary_Gland',
    'Muscle_Skeletal',
    'Nerve_Tibial',
    'Ovary',
    'Pancreas',
    'Pituitary',
    'Prostate',
    'Skin_Not_Sun_Exposed_Suprapubic',
    'Skin_Sun_Exposed_Lower_leg',
    'Small_Intestine_Terminal_Ileum',
    'Spleen',
    'Stomach',
    'Testis',
    'Thyroid',
    'Uterus',
    'Vagina',
    'Whole_Blood'
    ]

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

# Function to get top values
def get_top(df, cutoff, mode = "snp", pbsw_correction = True):

    df = df.copy()
    
    grp = ['CHR', 'POS', 'SNP']
    p_cols = ['PVALUE', 'LOG_PVALUE']
    gene_info_cols = ['GENE', 'GENE_ID', 'FUNCTION']
    
    for method in ['IHS', 'XPEHH', 'PBS']:
        if method in df.columns:
            df.rename(columns={method:'SCORE'}, inplace=True)
    
    if mode == 'snp':
        df = df[df.SCORE >= df.SCORE.quantile(.999)]
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

        df['PVALUE'] = (1 - df.SCORE.rank(method='max') / (df.SCORE.size+1)).round(4)
        df['LOG_PVALUE'] = -np.log10(df['PVALUE'])
        
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
        
    return df


def get_eQTL(snp, gene, target=None, p=0.01, datasetId='gtex_v7'):
    
    res = {}
    
    for tissue in tissues:
        print(f'> Getting eQTL for tissue {tissue}...')
        server = 'https://gtexportal.org/rest/v1/'
        ext = f'association/dyneqtl?gencodeId={gene}&variantId={snp}&tissueSiteDetailId={tissue}&datasetId={datasetId}'
        r = requests.get(server+ext, headers={"Accept" : "application/json"})
        if not r.ok:
            print(f'--- Request for SNP \33[4;33m"{snp}" and "{gene}"\33[m returned \33[1;4;31man error!\33[m ---')
            continue
        decoded = r.json()
        r.close()
        if 'pValue' not in decoded.keys():
            print(f"> {decoded['message']}")
            continue
        if decoded['pValue'] <= p:
            if target is not None:
                decoded['target'] = target
            res[tissue] = decoded
        else:
            continue
            
    return res


def get_eqtl_data(data):
    alleles = data['variantId'].split('_')[-3:-1]
    genotypes = np.where (
        np.array(data['genotypes']) == 0, alleles[0]+alleles[0],
        np.where(np.array(data['genotypes']) == 1, alleles[0]+alleles[1],
        alleles[1]+alleles[1])
        )
    df = pd.DataFrame(dict(Data=data['data'], Genotype=genotypes))
    return df


