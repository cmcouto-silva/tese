# SETTINGS

#### Load requried libraries

import os
import gzip
import shutil
import pickle
import numpy as np
import pandas as pd
import urllib.request
import datatable as dt
from collections import defaultdict
from more_itertools import consecutive_groups

import goatools
from goatools.mapslim import mapslim
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

#### Set working directory

os.chdir('dataset_S2/goatools')


#### Functions

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
    df = df[df.GENE_ID!="NA"]
    
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
    
    genes = df.query("GENE_ID!=''").GENE_ID.str.split(',').explode()
    genes = sorted(set(genes[genes!='NA'].astype(int)))
    return genes


# Download requried files
def download_file(url, out):
    if not os.path.exists(out):
        if url.endswith('.gz'):
            urllib.request.urlretrieve(url, out + '.gz')
            with gzip.open(out + '.gz', 'rb') as file_gz:
                with open(out, 'wb') as file_txt:
                    shutil.copyfileobj(file_gz, file_txt)
            os.remove(out + '.gz')    
        else:
            urllib.request.urlretrieve(url, out)


# Map multiple go terms to their ancestors
def apply_mapslim(assoc, only_direct=True):
    assoc_slim = {}
    for gene_id, go_terms in assoc.items():
        gene_direct_anc = set()
        gene_covered_anc = set()
        gene_all_anc = set()
        for go_term in go_terms:
            if go_term not in go_dag:
                continue
            term_direct_anc, term_all_anc = mapslim(go_term, go_dag, goslim_dag)
            gene_all_anc |= term_all_anc
            gene_covered_anc |= (term_all_anc - term_direct_anc)
        gene_direct_anc = gene_all_anc - gene_covered_anc
        if only_direct:
            assoc_slim[gene_id] = gene_direct_anc
        else:
            assoc_slim[gene_id] = gene_all_anc
    return assoc_slim


# Get a nested dict object
def nested_dict():
    return defaultdict(nested_dict)


# Download required data

os.makedirs('data', exist_ok=True)

# Download required file
download_file('http://geneontology.org/ontology/go-basic.obo', 'data/go-basic.obo')
download_file('http://www.geneontology.org/ontology/subsets/goslim_generic.obo', 'data/goslim_generic.obo')
download_file('ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz', 'data/gene2go.txt')

# Write human association file to run via command line
if not os.path.exists('data/association.txt'):
    gene2go = dt.fread('data/gene2go.txt')
    assoc_human = dt.fread('data/gene2go.txt')[dt.f['#tax_id']==9606, ['GeneID', 'GO_ID']]
    assoc_human = assoc_human.to_pandas().groupby('GeneID', as_index=False).agg(';'.join)
    assoc_human.to_csv('data/association.txt', sep='\t', index=False, header=False)

# Load required data (ontologies + association)
go_dag = GODag('data/go-basic.obo')
goslim_dag = GODag('data/goslim_generic.obo')
assoc = Gene2GoReader('data/gene2go.txt').get_ns2assc()

# Get GO slim terms
assoc_slim = {k: apply_mapslim(assoc[k]) for k in assoc}

# Check out how many annotated genes by GO namespace
for nspc, id2gos in assoc.items():
    print("{NS} {N:,} annotated human genes".format(NS=nspc, N=len(id2gos)))

## Load data

# Set Arguments
modes = ['persnp', 'pergene']
methods = ['pbs_windowed', 'xpehh', 'ihs']
cutoff_vals = [.999, .995, .99, .95]

# Load results
results = {mode: {
    method: pd.read_feather('../results_bin/{}/{}.feather'.format(mode, method)) for method in methods
    } for mode in modes }

# Load gene reference
ref = pd.read_csv('../data/annot/gene_keys.annot').dropna().reset_index(drop=True)
ref['ENTREZID'] = ref.ENTREZID.astype(int)

# Set common gene background
background_genes = ref.ENTREZID.sort_values().tolist()

# Running GO Enrichment Analysis

fields = ['GO','NS','depth','name','ratio_in_study','ratio_in_pop','p_uncorrected','p_fdr','p_fdr_bh']

go_all = nested_dict()
go_slim = nested_dict()

# goeaobj = GOEnrichmentStudy(background_genes, assoc['BP'], go_dag, methods=['fdr','fdr_bh'], propagate_counts = False)
# goeaobj_slim = GOEnrichmentStudy(background_genes, assoc['BP'], goslim_dag, methods=['fdr','fdr_bh'], propagate_counts = False)

goeaobj = GOEnrichmentStudyNS(background_genes, assoc, go_dag, methods=['fdr','fdr_bh'], propagate_counts = True)
goeaobj_slim = GOEnrichmentStudyNS(background_genes, assoc, goslim_dag, methods=['fdr','fdr_bh'], propagate_counts = True)

for mode in modes:
    for method in methods:
        for cutoff in cutoff_vals:
            print(f'{mode}_{method}_{cutoff}')
            target_genes = get_top_genes(results[mode][method], cutoff, mode)
            
            goea_all = goeaobj.run_study(target_genes, prt=None)
            goea_all_tbl = pd.DataFrame([goea_all[i].get_field_values(fields) for i in range(len(goea_all))], columns=fields)

            goea_slim = goeaobj_slim.run_study(target_genes, prt=None)
            goea_slim_tbl = pd.DataFrame([goea_slim[i].get_field_values(fields) for i in range(len(goea_slim))], columns=fields)

            go_all[mode][method][cutoff] = goea_all_tbl.copy()
            go_slim[mode][method][cutoff] = goea_slim_tbl.copy()


## Immune GO Terms

go_immune_term = go_dag['GO:0002376']

obodag_immune = {'GO:0002376': go_dag['GO:0002376']}

# # Only direct terms
# for child in [go.id for go in go_immune_term.children]: #go_immune_term.get_all_children():
#     obodag_immune[child] = go_dag[child]

# All immune terms
for child in [term for term in go_immune_term.get_all_children()]: #go_immune_term.get_all_children():
    obodag_immune[child] = go_dag[child]

# Verify results
#[x.name for x in obodag_immune.values()]

go_immune = nested_dict()

goeaobj_immune = GOEnrichmentStudyNS(background_genes, assoc, obodag_immune, methods=['fdr','fdr_bh'], propagate_counts = True)

for mode in modes:
    for method in methods:
        for cutoff in cutoff_vals:
            target_genes = get_top_genes(results[mode][method], cutoff, mode)
            goea_immune = goeaobj_immune.run_study(target_genes, prt=None)
            goea_immune_tbl = pd.DataFrame([goea_immune[i].get_field_values(fields) for i in range(len(goea_immune))], columns=fields)
            go_immune[mode][method][cutoff] = goea_immune_tbl.copy()


# Write binary files
os.makedirs('results', exist_ok=True)
results = {'go_all': go_all, 'go_slim': go_slim, 'go_immune': go_immune}
for k in results:
    with open(f'{k}.pickle', 'wb') as fh:
        pickle.dump(results[k], fh)


# Write csv files
for k in results:
    for mode in modes:
        for method in methods:
            for cutoff in cutoff_vals:
                results[k][mode][method][cutoff].to_csv(f'results/{k}_{mode}_{method}_{cutoff}.csv', index=False)

