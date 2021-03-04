## Setup functions

import re
import numpy as np
import pandas as pd
import requests, sys
from time import sleep
from collections import defaultdict
from more_itertools import consecutive_groups

#from selenium.webdriver import Chrome
#from selenium.webdriver.common.by import By
#from selenium.webdriver.chrome.options import Options
#from selenium.webdriver.support.wait import WebDriverWait
#from selenium.webdriver.common.action_chains import ActionChains
#from selenium.webdriver.support import expected_conditions as EC


def get_gene_studies(gene, gene_keys, set_style=True):
    
    gene = gene_keys[gene]
    headers = {'content-type': 'application/x-www-form-urlencoded'}
    params = f'ids={gene}&fields=generif'

    r = requests.post('http://mygene.info/v3/gene', data=params, headers=headers)
    df = pd.DataFrame(r.json()[0]['generif'])
    r.close()

    if set_style:
        df = df.assign(Title = df.text)
        df = df[['Title']].style.set_properties(**{'text-align': 'left', 'font-size': '12px', 'padding': '20px 10px', 'color': 'white'})
    return df


def separate_genes(df, gene_col = 'GENE', keep_empty = False):
    if not isinstance(gene_col, list):
        gene_col = [gene_col]
    if not keep_empty:
        df = df[(df.GENE != '') | (df.GENE.isna())].reset_index(drop=True)
    df = (df.assign(**{col: df[col].str.split(',') for col in gene_col})
          .apply(pd.Series.explode)
          .reset_index(drop=True))
    return df


def rm_undef_genes(df, pattern='^LINC|^LOC|^MIR\\d|^C.*or.*|^$'):
    undef = df.GENE.str.contains(pattern)
    return df[~undef]


def get_top_genes(df, cutoff, mode = "persnp", pbsw_correction = True, out = 'genes'):
    
    if mode not in ['persnp', 'pergene']:
        raise ValueError("mode must be 'persnp' or 'pergene'")
    if out not in ['dataframe','genes']:
        raise ValueError("out must be 'dataframe' or 'genes'")
                         
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
    
    if out == 'dataframe':
        out = df
    else:
        out = sorted(set(df.query("GENE!=''").GENE.str.split(',').explode()))
        
    return out


def get_genes_under_selection(results, genes):
    return {k: sorted(set(v) & set(genes)) for k,v in results.items()}


def get_intersection(data, thres=0.01, target='gene'):
    
    intersection = set()
    target = target.upper()
    
    if isinstance(data, dict):
        data = [df for df in data.values()]

    # SNP
    if target == 'SNP':
        for df in data:
            df = get_top(df, thres)
            if len(intersection) == 0:
                intersection.update(df[target])
            else:
                intersection.intersection_update(df[target])
    # GENE
    if target == 'GENE':
        for df in data:
            df = get_top_genes(df, thres)
            if len(intersection) == 0:
                intersection.update(df[target])
            else:
                intersection.intersection_update(df[target])
    
    return intersection        


def open_genecards(genes):
    for gene in genes:
        webbrowser.open(f'https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}')


def set_trait_label(efotraits, ref_labels='manuscript/studies/traits.xlsx', ignore_traits=[]):
    efotraits_labels = pd.read_excel('manuscript/studies/traits.xlsx')
    
    efotraits_tbl = (pd.concat(efotraits.values())[['MAPPED_TRAIT']]
                     .drop_duplicates()
                     .rename({'MAPPED_TRAIT':'Efotrait'}, axis=1)
                     .sort_values('Efotrait')
                     .drop_duplicates()
                    )
    
    efotraits_labeled = {}
    efotraits_tbl = efotraits_tbl.merge(efotraits_labels, how='outer').sort_values('Efotrait')
    
    if efotraits_tbl.Label.isna().any():
        print("There's NA in the \"Label\" column. Your file has been updated, please fix it.")
        efotraits_tbl.to_excel(ref_labels, index=False)
        return None
    
    for k in efotraits:
        merge = efotraits_tbl.merge(efotraits[k], left_on='Efotrait', right_on='MAPPED_TRAIT').drop('MAPPED_TRAIT', axis=1)
        merge = merge.assign(Label=merge.Label.str.split(','))
        merge.Label = merge.Label.apply(lambda lst: [x.strip() for x in lst])
        merge = merge.explode('Label').groupby('Label', as_index=False).sum()
        merge = merge[merge.Label!='undef'].sort_values('N', ascending=False).reset_index(drop=True)
        if len(ignore_traits) > 0:
            merge = merge[~merge.Label.isin(ignore_traits)]
        if merge.size > 0:
            efotraits_labeled[k] = merge
            
    return efotraits_labeled


def get_pmid_data(genes, efotraits):
    pmid_dict = {}
    for gene in efotraits_labeled:
        pmid_dict[gene] = {}
        pmid_dict[gene]['traits'] = efotraits_labeled[gene]
        pmid_dict[gene]['ihs_S1'] = gene in genes['ihs_s1']
        pmid_dict[gene]['xpehh_S1'] = gene in genes['xpehh_s1']
        pmid_dict[gene]['pbs_S1'] = gene in genes['pbs_s1']
        pmid_dict[gene]['ihs_S2'] = gene in genes['ihs_s2']
        pmid_dict[gene]['xpehh_S2'] = gene in genes['xpehh_s2']
        pmid_dict[gene]['pbs_S2'] = gene in genes['pbs_s2']
        
    return pmid_dict


def make_tableS2(gene_sel, pmid, inplace=True):
    uniq_genes = np.sort(np.unique(sum([v for v in gene_sel.values()], [])))
    col_names = ['PMID', 'PBS_S1', 'XP-EHH_S1', 'iHS_S1',
            'PBS_S2', 'XP-EHH_S2', 'iHS_S2', 'N_GENES', 'GENES']

    row = pd.Series([
        pmid,
        len(gene_sel['pbs_s1']),
        len(gene_sel['xpehh_s1']),
        len(gene_sel['ihs_s1']),
        len(gene_sel['pbs_s2']),
        len(gene_sel['xpehh_s2']),
        len(gene_sel['ihs_s2']),
        len(uniq_genes),
        ','.join(uniq_genes)
        ], index=col_names, name=pmid)
    
    try:
        global tableS2
        tableS2 = tableS2.append(row).drop_duplicates()
    except:
        tableS2 = pd.DataFrame(row).T.drop_duplicates()
    
    return tableS2


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
              .query("Phenotype != 'Annotated by HGMD'")
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
        gwastable_loc_buttom.click()
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
    if database == "gwas":
        efolist = gene2gwas(gene, verbose)
    elif database == "ensembl":
        efolist = gene2ensembl(gene, verbose)
    else:
        efolist = gene2genecard(gene, sleep_time, verbose)
    return efolist


def efolist2efotable(efolist, gene_sel):
    df = pd.concat(efolist).reset_index()
    if "MAPPED_TRAIT" in df.columns:
        df.rename(columns={"GENE": "Gene", "MAPPED_TRAIT": "Phenotype"}, inplace=True)
    df = df.groupby('Phenotype', as_index=False).agg({'Gene': lambda x: ','.join(x)})
    df = (df.assign(Genes = df.Gene, N = df.Gene.str.count(',') + 1)
         .drop('Gene', 1)
         .sort_values(['N', 'Phenotype'], ascending=[False, True])
                    .reset_index(drop=True))
    df = include_methods(df, gene_sel)
    return df


def include_methods(efotable, gene_sel):
    df = efotable.copy()
    pheno_genes = df.Genes.str.split(',')
    gene_sel_uniq = set(sum(gene_sel.values(), []))
    
    ref = pd.DataFrame()
    for gene in sorted(gene_sel_uniq):
        for method, genes in gene_sel.items():
            ref.loc[gene, method] = gene in genes 
    
    genes_by_method = []
    for genes in pheno_genes:
        genes_by_method.append(ref[ref.index.isin(genes)].sum())
    genes_by_method = pd.DataFrame(genes_by_method)
    genes_by_method = genes_by_method.assign(sum_stat = genes_by_method.sum(1)).astype(int)
    df = df.join(genes_by_method)
    df = (df.assign(n_genes=df.N, genes=df.Genes)
          .drop(['Genes', 'N'], 1))
    return df

