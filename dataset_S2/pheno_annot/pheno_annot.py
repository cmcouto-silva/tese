# SETTINGS

#### Load required libraries

import os
import re
import pickle
import numpy as np
import pandas as pd
import requests, sys
from time import sleep
from collections import defaultdict
from more_itertools import consecutive_groups

# from selenium.webdriver import Chrome
# from selenium.webdriver.common.by import By
# from selenium.webdriver.chrome.options import Options
# from selenium.webdriver.support.wait import WebDriverWait
# from selenium.webdriver.common.action_chains import ActionChains
# from selenium.webdriver.support import expected_conditions as EC

#### Set working directory

os.chdir('dataset_S2/pheno_annot')

#### Functions

# Function to split genes and related info
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
    
    genes = sorted(set(df.query("GENE!=''").GENE.str.split(',').explode()))
    return genes


# Get a nested dict object
def nested_dict():
    return defaultdict(nested_dict)

def gene2gwas(gene, verbose=True):
    if verbose: print(f'Getting annotation for gene \033[1m{gene}\033[m...')
    url = f'https://www.ebi.ac.uk/gwas/api/search/downloads?q=ensemblMappedGenes:{gene}&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=&genotypingfilter[]=&traitfilter[]=&dateaddedfilter=&facet=association&efo=true'
    efotraits = pd.read_csv(url, sep='\t')[['MAPPED_GENE', 'MAPPED_TRAIT', 'REPORTED GENE(S)', 'DISEASE/TRAIT']]
    if efotraits.shape[0]==0:
        if verbose:
            print(f'> Gene {gene} \033[1;4;31mnot available\033[m or doens\'t have associated traits on GWASCatalog!')
        return None
    efotraits = efotraits[efotraits.MAPPED_GENE.str.replace(' ', '').str.split(',|;|-| ').apply(lambda x: gene in x)]
    efotraits = (
        efotraits[['MAPPED_GENE','MAPPED_TRAIT']]
        .assign(MAPPED_TRAIT = efotraits.MAPPED_TRAIT.str.split(','))
        .explode('MAPPED_TRAIT')
        .assign(MAPPED_TRAIT = lambda x: x['MAPPED_TRAIT'].str.strip())
        .groupby('MAPPED_TRAIT', as_index=False).count()
        .rename({'MAPPED_GENE': 'N_TRAITS'}, axis=1)
        .assign(GENE=gene)
        .set_index('GENE')
        .sort_values(['N_TRAITS', 'MAPPED_TRAIT'], ascending=[False, True])
        )
    return(efotraits)

def gene2ensembl(gene, verbose=True):
    if verbose: print(f'Getting annotation for gene \033[1m{gene}\033[m...')
    server = "http://rest.ensembl.org"
    ext = f"/phenotype/gene/homo_sapiens/{gene}?include_associated=1" 
    try:
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        decoded = r.json()
        r.close()
        df = (
            pd.DataFrame(decoded).groupby('description', as_index=False)['source'].count()
              .query("description!='ClinVar: phenotype not specified'")
              .rename(columns={'description': 'Phenotype', 'source': 'N'})
              .sort_values(['N', 'Phenotype'], ascending=[False, True])
              .assign(Gene=gene)
              .set_index('Gene')
            )
    except:
        print(f"Couldn't reach gene {gene}.")
        df = None
    return df

def gene2genecard(gene, sleep_time=1, verbose=True):
    opts = Options()
    opts.headless = True
    browser = Chrome(options=opts)
    browser.maximize_window()
    actions = ActionChains(browser)
    if verbose: print(f'Getting phenotypes for gene \033[1m{gene}\033[m...')

    browser.get(f'https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}')

    if len(browser.find_elements(By.XPATH, '//*[@id="cardPage"]')) == 0:
        if verbose: print(f'> Gene \033[1;41m{gene}\033[m not available on GeneCards!\n')
        browser.close()
        return None
    
    if len(browser.find_elements(By.XPATH, '//*[@id="function-phenotypes-from-gwas"]')) == 0:
        if verbose: print(f'> Gene \033[1;33m{gene}\033[m has no GWAS information!\n')
        browser.close()
        return None
    
    gwastable_loc = browser.find_element_by_xpath('//*[@id="function-phenotypes-from-gwas"]')
    actions.move_to_element(gwastable_loc).perform()

    table_size = browser.find_element_by_xpath('//*[@id="function-phenotypes-from-gwas"]/div/div/div[2]/div/div[3]/form/small[2]/em')
    non_decimal = re.compile(r'[^\d.]+')
    table_size = int(non_decimal.sub('', table_size.text))
    table_size

    if table_size > 5:
        see_all = '//*[@id="function-phenotypes-from-gwas"]/div/div/div[2]/div/div[1]/a[1]'
        see_less = '//*[@id="function-phenotypes-from-gwas"]/div/div/div[2]/div/div[1]/a[2]'
        gwastable_loc_buttom = browser.find_element_by_xpath(see_all)
        browser.execute_script("arguments[0].click();", gwastable_loc_buttom)
        WebDriverWait(browser, 3).until(EC.element_to_be_clickable((By.XPATH, see_less)))
        sleep(sleep_time)

    gwas_tbl = browser.find_element_by_css_selector('#function-phenotypes-from-gwas > div > div > div.ng-scope > div > table')
    gwas_tbl = pd.read_html(gwas_tbl.get_attribute('outerHTML'))[0]
    gwas_tbl['SNP IDs'] = gwas_tbl['SNP IDs'].apply(lambda x: x[:x.find('See all')].strip().replace(' ', ','))
    gwas_tbl['Gene Relation'] = gwas_tbl['Gene Relation'].str.replace(' ', '')
    
    if table_size != gwas_tbl.shape[0]:
        if verbose: print('> \033[1;43mNot all traits included\033[m. Please consider increasing the sleeping time.')
        browser.close()
        return None
    
    gwas_tbl = gwas_tbl.assign(Phenotype = gwas_tbl.Phenotype.str.split(',')).explode('Phenotype')
    
    gwas_tbl = (gwas_tbl.groupby('Phenotype', as_index=False).agg({
        'Gene Relation': lambda x: ','.join(np.unique(x)),
        'Best Score': np.max,
        'Mean Score': np.mean,
        '# of Snps': lambda x : sum(x) * 0,
        'SNP IDs': lambda x: ','.join(np.unique(sum(map(lambda x: x.split(','), x), [])))
    }))
    
    gwas_tbl = gwas_tbl.assign(Gene=gene).set_index('Gene')
    gwas_tbl['# of Snps'] = gwas_tbl['SNP IDs'].str.count(',') + 1
    
    print(f'> Phenotype successfully aquired!\n')
    browser.close()
    return gwas_tbl


def gene2phenotype(gene, database='gwas', verbose=True, sleep_time=1):
    """
    Database can ben "gwas", "ensembl", or "genecard".
    """
    if database not in ['gwas', 'ensembl', 'genecard']:
        raise ValueError('database must be "gwas", "ensembl", or "genecard".')
    if database == "gwas":
        efolist = gene2gwas(gene, verbose)
    elif database == "ensembl":
        efolist = gene2ensembl(gene, verbose)
    else:
        efolist = gene2genecard(gene, sleep_time, verbose)
    return efolist


def efolist2efotable(efolist):
    df = pd.concat(efolist).reset_index()
    if "MAPPED_TRAIT" in df.columns:
        df.rename(columns={"GENE": "Gene", "MAPPED_TRAIT": "Phenotype"}, inplace=True)
    df = df.groupby('Phenotype', as_index=False).agg({'Gene': lambda x: ','.join(x)})
    df = (df.assign(Genes = df.Gene, N = df.Gene.str.count(',') + 1)
         .drop('Gene', 1)
         .sort_values(['N', 'Phenotype'], ascending=[False, True])
                    .reset_index(drop=True))
    return df


## Load data

# Set Arguments
cutoff = .999
modes = ['persnp', 'pergene']
methods = ['pbs_windowed', 'xpehh', 'ihs']

# Load results
results = {mode: {
    method: pd.read_feather('../results_bin/{}/{}.feather'.format(mode, method)) for method in methods
    } for mode in modes }

os.getcwd()

# Load gene reference
ref = pd.read_csv('../data/annot/gene_keys.annot').dropna().reset_index(drop=True)
ref['ENTREZID'] = ref.ENTREZID.astype(int)

# Get top genes per mode-method
results_topgenes = defaultdict(dict)

for mode in modes:
    for method in methods:
            results_topgenes[mode][method] = get_top_genes(results[mode][method], cutoff, mode)

## Annotation

# Using GWAS Database
gwas_annot =  defaultdict(dict)

for mode in modes:
    print('-'* 50)
    print(f' Mode "{mode}" '.center(50, '-'))
    print('-'* 50, end='\n\n')
    for method in methods:
        print(f' Method {method.upper()} '.center(50, '-'))
        gwas_annot[mode][method] = efolist2efotable([gene2phenotype(gene, database='gwas') for gene in results_topgenes[mode][method]])
        print()

# Using Ensembl Database
ensembl_annot =  defaultdict(dict)

for mode in modes:
    print('-'* 50)
    print(f' Mode "{mode}" '.center(50, '-'))
    print('-'* 50, end='\n\n')
    for method in methods:
        print(f' Method {method.upper()} '.center(50, '-'))
        ensembl_annot[mode][method] = efolist2efotable([gene2phenotype(gene, database='ensembl') for gene in results_topgenes[mode][method]])
        print()

# # Using GeneCard Database
# genecard_annot =  defaultdict(dict)

# for mode in modes:
#     print('-'* 50)
#     print(f' Mode "{mode}" '.center(50, '-'))
#     print('-'* 50, end='\n\n')
#     for method in methods:
#         print(f' Method {method.upper()} '.center(50, '-'))
#         genecard_annot[mode][method] = efolist2efotable([gene2phenotype(gene, database='genecard', sleep_time=3) for gene in results_topgenes[mode][method]])
#         print()

results_annot = dict(gwas = gwas_annot, ensembl = ensembl_annot)#, genecard = genecard_annot)

with open('results_annot.pickle', 'wb') as file:
    pickle.dump(results_annot, file)

