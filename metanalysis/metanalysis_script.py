# Setup

### Setup libraries

import os
import sys
import requests
import numpy as np
import pandas as pd
pd.set_option('max_colwidth', None)

sys.path.insert(0, 'metanalysis')
from custom_functions import *

### Setup data

# Getting data
results = {
    'ho': {
        'pbs': pd.read_feather('dataset_S1/results_bin/persnp/pbs_windowed.feather'),
        'xpehh': pd.read_feather('dataset_S1/results_bin/persnp/xpehh.feather'),
        'ihs': pd.read_feather('dataset_S1/results_bin/persnp/ihs.feather')
    },
    'illumina': {
        'pbs': pd.read_feather('dataset_S2/results_bin/persnp/pbs_windowed.feather'),
        'xpehh': pd.read_feather('dataset_S2/results_bin/persnp/xpehh.feather'),
        'ihs': pd.read_feather('dataset_S2/results_bin/persnp/ihs.feather')
    }
}

for geno_platform in results:
    for method in results[geno_platform]:
        results[geno_platform][method].dropna(inplace=True)
        results[geno_platform][method].reset_index(drop=True, inplace=True)

# Load gene_keys for retrieving papers
annot_S1 = pd.read_csv('dataset_S1/data/annot/gene_keys.annot')
annot_S2 = pd.read_csv('dataset_S2/data/annot/gene_keys.annot')
gene_keys = pd.concat([annot_S1, annot_S2]).dropna().drop_duplicates()
gene_keys['ENTREZID'] = gene_keys['ENTREZID'].astype(int)
gene_keys = gene_keys.set_index('SYMBOL').to_dict()['ENTREZID']

# Filtering data
cutoff = .995
results_outliers = dict(
    pbs_s1 = get_top_genes(results['ho']['pbs'], cutoff),
    xpehh_s1 = get_top_genes(results['ho']['xpehh'], cutoff),
    ihs_s1 = get_top_genes(results['ho']['ihs'], cutoff),
    
    pbs_s2 = get_top_genes(results['illumina']['pbs'], cutoff),
    xpehh_s2 = get_top_genes(results['illumina']['xpehh'], cutoff),
    ihs_s2 = get_top_genes(results['illumina']['ihs'], cutoff)
)

# Metanalysis

# Set object for storing all results:

#studies = {}
#phenotype_dataset = 'ensembl' # available options: 'gwas', 'ensembl', 'genecard'

for phenotype_dataset in ['gwas', 'ensembl']:

    studies = {}
    print(f'INITIALIZING analysis with {phenotype_dataset.upper()} dataset')
    print('-'*40)
    print()

    ### Amorim et al. (2015) ─ PLoS ONE

    pmid = 25849546

    written_on_paper = "ABLIM3, ACSS2, AKAP6, ANKRD26, ARHGEF10, ATIC, BCAT1, C20orf111, C2orf73, CBLN1, CNTN4, CNTNAP5, COL22A1, CPA5, CRTC3, CWH43, DCUN1D4, DHCR7, EPHB4, FAM188B, FKBP6, GALNT16, GLIS3, GLRB, GPC6, GRIK2, HLA-DPA1, HMG20B, HSF2, IQGAP1, KIAA1598, KLHL29, LRRC66, MASTL, MPST, NAALADL2, NRG1, NRP2, NSUN5, PARK2, PKIB, PPP2R2C, RABGGTB, RAD51B, RBFOX1, RBM9, RBMS3, RFX3, ROBO2, SCP2, SGCB, SHISA6, SPAG16, SPATA13, SPATA18, ST6GAL1".split(',')
    written_on_tableS2 = "MAST2,LAMC1,TTLL4,CCDC108,ITSN2,QRFRP,SH3TC2,TRIP6,ZAN,ZAN,ZAN,ECD,ECD,TTC18,PUS3,FGD6,FGD6,ZNF862,LIMK2".split(',')
    genes = written_on_paper + written_on_tableS2
    genes = set([gene.strip() for gene in genes])

    # Get gene intersection with our analysis
    gene_sel = get_genes_under_selection(results_outliers, genes)
    gene_sel_uniq = sorted(set(sum(gene_sel.values(), [])))
    make_tableS2(gene_sel, pmid)

    efolist = [gene2phenotype(gene, phenotype_dataset, sleep_time=3) for gene in gene_sel_uniq]
    efotable = efolist2efotable(efolist, gene_sel)

    studies[pmid] = dict(genes=gene_sel, efolist=efolist, efotable=efotable)

    ### Bergey et al. (2018) ─ PNAS

    pmid = 30413626

    url = 'https://www.pnas.org/highwire/filestream/836471/field_highwire_adjunct_files/1/pnas.1812135115.sd01.xlsx'
    sheets = ['Dataset S1a PBS Outlier SNPs', 'Dataset S1c PBS Outlier Genes', 'Dataset S1f Bayenv Outlier SNPs', 'Dataset S1g Bayenv Outlier Genes']

    bergey_data = pd.read_excel(url, sheets, skiprows=1)

    S1a = set(bergey_data['Dataset S1a PBS Outlier SNPs'][:52].dropna().Gene)
    S1c = set(bergey_data['Dataset S1c PBS Outlier Genes'][:49].dropna().Gene.str.replace('"', ''))
    S1f = set(bergey_data['Dataset S1f Bayenv Outlier SNPs'].Gene)
    S1g = set(bergey_data['Dataset S1g Bayenv Outlier Genes'].Gene)

    genes = sorted(set.union(S1a, S1c, S1f, S1g))

    # Get gene intersection with our analysis
    gene_sel = get_genes_under_selection(results_outliers, genes)
    gene_sel_uniq = sorted(set(sum(gene_sel.values(), [])))
    make_tableS2(gene_sel, pmid)

    efolist = [gene2phenotype(gene, phenotype_dataset, sleep_time=3) for gene in gene_sel_uniq]
    efotable = efolist2efotable(efolist, gene_sel)

    studies[pmid] = dict(genes=gene_sel, efolist=efolist, efotable=efotable)

    ### Harrison et al. (2019) ─ Nat Ecol Evol

    pmid = 31358949

    df = pd.read_excel('https://static-content.springer.com/esm/art%3A10.1038%2Fs41559-019-0947-6/MediaObjects/41559_2019_947_MOESM9_ESM.xlsx', engine='openpyxl', skiprows=1)
    df = df.rename({'HUGO':'Gene'}, axis=1).iloc[:, :-1]

    genes = df[(df.Twa_PBS >= df.Twa_PBS.quantile(.99)) | (df.iHS_Standardized_Batwa.abs() >= df.iHS_Standardized_Batwa.abs().quantile(.99))]
    genes = genes.Gene.unique()

    # Get gene intersection with our analysis
    gene_sel = get_genes_under_selection(results_outliers, genes)
    gene_sel_uniq = sorted(set(sum(gene_sel.values(), [])))
    make_tableS2(gene_sel, pmid)

    efolist = [gene2phenotype(gene, phenotype_dataset, sleep_time=3) for gene in gene_sel_uniq]
    efotable = efolist2efotable(efolist, gene_sel)

    studies[pmid] = dict(genes=gene_sel, efolist=efolist, efotable=efotable)

    ### Hsieh et al. (2016)

    pmid = 28962010

    # Select target genes from paper
    genes_tableS4 = ['ABCC8', 'ADAR', 'ADCK3', 'APOBEC3', 'APOBEC3A', 'APOBEC3B', 'APOBEC3BAS1', 'APOBEC3C', 'APOBEC3F', 'APOBEC3G', 'ARF1', 'ATP2B2', 'BRF2NCR3LG1', 'BRK1', 'BTNL10', 'C1orf145', 'C1orf35', 'C2CD4D', 'C3orf36', 'C5orf52', 'CDC42BPAWNT9A', 'CDRT15', 'CHRNB2', 'CLINT1NDUFA4', 'COX10', 'COX10-AS1', 'CYFIP2SOX30', 'DPYD', 'DPYD-AS1', 'DUSP5P1MIR4776-2', 'ERLIN2', 'FANCD2OS', 'GALNTL5XRCC2UNC5DZNF703', 'GAREMLARGEAPOBEC3A_B', 'GHRL', 'GHRLOS', 'GJC2', 'GPR124', 'GUK1', 'HAVCR1ITK', 'HIST3H2A', 'HIST3H2BB', 'HIST3H3', 'HS3ST3B1', 'IBA57', 'IBA57AS1', 'IKZF2AGAP1FANCD2', 'IRAK2TATDN2', 'KCNJ11', 'KCNN3ACBD3', 'LIN9', 'LINC00852', 'LINGO4', 'LOC100132111', 'LOC728024', 'LSM11', 'MGC12916MEP1B', 'MIR3620', 'MIR4666A', 'MIR4776-1', 'MIR5008', 'MIR6742', 'MIR885PLCXD2', 'MIXL1', 'MRPL55', 'NBPF18P', 'OBSCN', 'OTOGESRRBCDRT15P1', 'PARP1ITPKB', 'PHF14PRKAG2', 'PHLDB2SRPRB', 'PPP1R2P3', 'PRKAG2-AS1', 'PROSC', 'PSEN2', 'RAB6B', 'RASA2TSPAN5SGCD', 'RHOU', 'RNF187', 'RORC', 'S100A10', 'S100A11', 'SEC13', 'SHE', 'SLCO2A1EPHB1ZBTB38', 'TCHHL1IL6R', 'TDRD10', 'TDRKH', 'THEM4', 'THEM5', 'THG1L', 'TIMD4', 'TRIM11', 'TRIM17', 'UBE2Q1', 'USH1C', 'VHL', 'WNT3A']
    genes_tableS5 = ['AXDND1', 'HLA-DPA1', 'HLA-DPB1ZNF44', 'LAMC1', 'LAMC2', 'RNU5F-1FLNBHLA-DOA', 'ZNF442', 'ZNF563', 'ZNF799']
    genes = sorted(set(genes_tableS4) | set(genes_tableS5))

    # Get gene intersection with our analysis
    gene_sel = get_genes_under_selection(results_outliers, genes)
    gene_sel_uniq = sorted(set(sum(gene_sel.values(), [])))
    make_tableS2(gene_sel, pmid)

    efolist = [gene2phenotype(gene, phenotype_dataset, sleep_time=3) for gene in gene_sel_uniq]
    efotable = efolist2efotable(efolist, gene_sel)

    studies[pmid] = dict(genes=gene_sel, efolist=efolist, efotable=efotable)

    ### Jarvis et al. (2012)

    pmid = 22570615

    lsbl = pd.read_excel('https://journals.plos.org/plosgenetics/article/file?type=supplementary&id=info:doi/10.1371/journal.pgen.1002641.s010', skiprows=2)
    lsbl_genes = set(lsbl.Genes.dropna().str.split(';').explode())

    xpehh = pd.read_excel('https://journals.plos.org/plosgenetics/article/file?type=supplementary&id=info:doi/10.1371/journal.pgen.1002641.s012', skiprows=2)
    xpehh_genes = set(xpehh[xpehh.XPEHH_raw >= 2].Genes.dropna().str.split(';').explode())

    ihs = pd.read_excel('https://journals.plos.org/plosgenetics/article/file?type=supplementary&id=info:doi/10.1371/journal.pgen.1002641.s013', skiprows=2)
    ihs_genes = set(ihs.Genes.dropna().str.split(';').explode())

    genes = set.union(lsbl_genes, xpehh_genes, ihs_genes)

    # Get gene intersection with our analysis
    gene_sel = get_genes_under_selection(results_outliers, genes)
    gene_sel_uniq = sorted(set(sum(gene_sel.values(), [])))
    make_tableS2(gene_sel, pmid)

    efolist = [gene2phenotype(gene, phenotype_dataset, sleep_time=3) for gene in gene_sel_uniq]
    efotable = efolist2efotable(efolist, gene_sel)

    studies[pmid] = dict(genes=gene_sel, efolist=efolist, efotable=efotable)

    ### Lachance et al. (2012) ─ Cell

    pmid = 22840920

    r = requests.get('https://www.cell.com/cms/10.1016/j.cell.2012.07.009/attachment/6f1acaf0-5727-4cb3-9e08-08d24f4668cb/mmc3.xlsx', headers={'User-Agent': 'Mozilla/5.0'})
    mmc3 = r.content
    r.close()

    lsbl = pd.read_excel(mmc3, sheet_name=['LSBL (Pygmy)', 'LSBL (Hadza)', 'LSBL (Sandawe)'], engine='openpyxl', skiprows=2)
    lsbl_set = set.union(*[set(df.iloc[:25, -1].str.split(';').explode()) for k,df in lsbl.items()])

    aim_clusters = pd.read_excel(mmc3, sheet_name=['AIM clusters (Pygmy)', 'AIM clusters (Hadza)', 'AIM clusters (Sandawe)'], skiprows=2)
    aim_clusters_set = set.union(*[set(df.dropna().iloc[:, -1].str.replace(' ', '').str.split(',').explode()) for df in aim_clusters.values()])

    genes = sorted(lsbl_set | aim_clusters_set)

    # Get gene intersection with our analysis
    gene_sel = get_genes_under_selection(results_outliers, genes)
    gene_sel_uniq = sorted(set(sum(gene_sel.values(), [])))
    make_tableS2(gene_sel, pmid)

    efolist = [gene2phenotype(gene, phenotype_dataset, sleep_time=3) for gene in gene_sel_uniq]
    efotable = efolist2efotable(efolist, gene_sel)

    studies[pmid] = dict(genes=gene_sel, efolist=efolist, efotable=efotable)

    ### López Herráez et al. (2009) ─ PLoS ONE

    pmid = 19924308

    df = pd.read_excel('https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0007888.s012', skiprows=1)
    df = df.rename({'HGNC gene ID': 'GENE', 'Population ID(rank)': 'POPULATION'}, axis=1)
    df['POPULATION'] = df.POPULATION.apply(lambda x: re.sub('\(.*?\)', '', x))

    pop = set(df.POPULATION.apply(lambda x: x.split(',')).sum())
    pop.difference_update({'Karitiana', 'Surui', 'Mbutipygmy', 'Biakapygmy'})
    pop = '|'.join(sorted(pop)).replace(' ', '')

    df = df[~df.POPULATION.apply(lambda x: bool(re.search(pop, x)))].reset_index(drop=True)
    df = df[~df.GENE.isna()]

    genes = sorted(set(df.GENE.apply(lambda x: x.split(',')).sum()))

    # Get gene intersection with our analysis
    gene_sel = get_genes_under_selection(results_outliers, genes)
    gene_sel_uniq = sorted(set(sum(gene_sel.values(), [])))
    make_tableS2(gene_sel, pmid)

    efolist = [gene2phenotype(gene, phenotype_dataset, sleep_time=3) for gene in gene_sel_uniq]
    efotable = efolist2efotable(efolist, gene_sel)

    studies[pmid] = dict(genes=gene_sel, efolist=efolist, efotable=efotable)

    ### Miligliano et al. (2013)

    pmid = 24297229
    # obs: the supplementary material link for this is corrupted

    # Select target genes from paper
    genes_on_table2 = set(
        'DAPP1, CSDA, CCNE1, NRIP1, LHCGR, IMMP2L, FOXO3, COL9A3, PDGFRA, MAS1, SHH, BCL2L11, NOBOX'.replace(' ', '').split(',') + 
        'COPS2, NR4A3, MED30, FOXE1, SULT1B1, TRIP4, CREB1'.replace(' ', '').split(',') +
        'KCNMA1, GUCY1A3, CSDA, NRIP1, KCNJ8, PTGTS1, CACNA1C, TRPA1, CIRBP, GCLC, AGTR1, SOD2, TNFRSF11A, ECE1, AVPR1A, CLPB'.replace(' ', '').split(',') +
        'RPS6KA1, JAK2, CSNK2A2, IGF2R, GHITM, PLCG1'.replace(' ', '').split(',')
    )

    genes_on_paper_corpus = set([
        'ZFAT', 'RNU4ATAC', 'NAALADL2', 'GHRHR', 'IGF2R', 'IGF2BP2', 'IGF2BP3', 'CSNK2A1', 'IGF2R', 'RPS6KA1', 'JAK2', 'CSNK2A2', 'PLCG1', 'GHITM',
        'TRIP4', 'IYD', 'FOXE1', 'MSTN', 'FOCS2', 'PPARA', 'CLPB', 'DMRT1', 'DMRT3', 'DAPP1', 'CSMD1', 'CDC42SE2'
    ])

    genes = genes_on_table2 | genes_on_paper_corpus

    # Get gene intersection with our analysis
    gene_sel = get_genes_under_selection(results_outliers, genes)
    gene_sel_uniq = sorted(set(sum(gene_sel.values(), [])))
    make_tableS2(gene_sel, pmid)

    efolist = [gene2phenotype(gene, phenotype_dataset, sleep_time=3) for gene in gene_sel_uniq]
    efotable = efolist2efotable(efolist, gene_sel)

    studies[pmid] = dict(genes=gene_sel, efolist=efolist, efotable=efotable)

    ### Perry et al. (2014) ─ PNAS

    pmid = 25136101

    # Select target genes from paper

    topgenes_table = 'https://www.pnas.org/highwire/filestream/616994/field_highwire_adjunct_files/2/pnas.1402875111.sd03.xlsx'
    batwa_sheet = 'Batwa iHS outlier regions'
    baka_sheet = 'Baka iHS outlier regions'

    df = pd.concat([
        pd.read_excel(topgenes_table, batwa_sheet, skiprows=2, na_values='NOT_FOUND', engine='openpyxl'),
        pd.read_excel(topgenes_table, baka_sheet, na_values='NOT_FOUND', engine='openpyxl')
    ])

    df.columns = ['POS','GENE']
    df.dropna(inplace=True)

    genes = set((df.assign(GENE=df.GENE.str.split(','))
    .explode('GENE')
    .GENE
    ))

    # Get gene intersection with our analysis
    gene_sel = get_genes_under_selection(results_outliers, genes)
    gene_sel_uniq = sorted(set(sum(gene_sel.values(), [])))
    make_tableS2(gene_sel, pmid)

    efolist = [gene2phenotype(gene, phenotype_dataset, sleep_time=3) for gene in gene_sel_uniq]
    efotable = efolist2efotable(efolist, gene_sel)

    studies[pmid] = dict(genes=gene_sel, efolist=efolist, efotable=efotable)

    ### Scheinfeldt et al. (2019) ─ PNAS

    pmid = 30782801

    target_pops = ['Boni','Dahalo','Elmolo','Hadza','Ogiek','Sabue','Sandawe','Sengwer','WRHG','Wata','Yaaku']

    df1 = pd.read_csv('https://www.pnas.org/highwire/filestream/849885/field_highwire_adjunct_files/1/pnas.1817678116.sd01.csv')
    df2 = pd.read_csv('https://www.pnas.org/highwire/filestream/849885/field_highwire_adjunct_files/2/pnas.1817678116.sd02.csv')
    df3 = pd.read_csv('https://www.pnas.org/highwire/filestream/849885/field_highwire_adjunct_files/3/pnas.1817678116.sd03.csv')

    def clean_data(df):
        df = (df.assign(gene=df.genes_within_100kb.str.replace(' ', '').str.split(';'))
        .explode('gene').query('gene!= ""')#.sort_values('gene')
        .dropna(subset=['gene'])
        )
        df = df[df['population grouping'].isin(target_pops)]
        return df

    df1 = clean_data(df1)
    df2 = clean_data(df2)
    df3 = clean_data(df3)

    df2['iHS_score'] = df2['iHS_score'].abs()
    df3['xp-clr'] = df3['xp-clr'].astype(float)

    genes_d = df1[df1.D >= df1.D.quantile(.99)].gene.unique()
    genes_ihs = df2[df2.iHS_score >= df2.iHS_score.quantile(.99)].gene.unique()
    genes_xpclr = df3[df3['xp-clr'] >= df3['xp-clr'].quantile(.99)].gene.unique()

    genes = np.concatenate([genes_d, genes_ihs, genes_xpclr])

    # Get gene intersection with our analysis
    gene_sel = get_genes_under_selection(results_outliers, genes)
    gene_sel_uniq = sorted(set(sum(gene_sel.values(), [])))
    make_tableS2(gene_sel, pmid)

    efolist = [gene2phenotype(gene, phenotype_dataset, sleep_time=3) for gene in gene_sel_uniq]
    efotable = efolist2efotable(efolist, gene_sel)

    studies[pmid] = dict(genes=gene_sel, efolist=efolist, efotable=efotable)

    ### Make tables

    # Get main efotable
    efotable = (
        pd.concat({study:studies[study]['efotable'] for study in studies})
        .reset_index(level=0)
        .reset_index(drop=1)
        .rename(columns={'level_0': 'PMID', 'genes': 'Genes'})
        .assign(PMID = lambda x: x['PMID'].astype(str))
        .filter(regex='^[A-Z]')
        .groupby('Phenotype', as_index=False)
        .agg(lambda x: ','.join(sorted(set(x))))
        .assign(N= lambda x: x['Genes'].str.count(',')+1)
        .sort_values('N', ascending=False)
        .reset_index(drop=True)
    )

    # Get genes from main efotable
    genes = sorted(set(sum(efotable.Genes.str.split(','), [])))
    gene_sel = get_genes_under_selection(results_outliers, genes)
    gene_sel_uniq = sorted(set(sum(gene_sel.values(), [])))

    # Create a Pandas Excel writer using XlsxWriter as the engine
    writer = pd.ExcelWriter(f'metanalysis/efotable_by_pmid_{phenotype_dataset}.xlsx', engine='xlsxwriter')

    # Write each dataframe to a different worksheet
    for study in studies:
        df = studies[study]['efotable']
        df.to_excel(writer, sheet_name=str(study), index=False)
        writer.sheets[str(study)].set_column(0, 0, 70)
        writer.sheets[str(study)].autofilter(0, 0, df.shape[0]-1, df.shape[1]-1)
    # Close the Pandas Excel writer and output the Excel file
    writer.save()

    # Get main table
    efotable_main = include_methods(efotable, gene_sel)
    efotable_main['n_pmid'] = efotable_main.PMID.str.count(',') + 1
    efotable_main['pmid'] = efotable_main.PMID
    efotable_main.drop(columns='PMID', inplace=True)

    # Sum genes in methods ignoring iHS results
    efotable_main.insert(efotable_main.columns.to_list().index('sum_stat')+1, 'sum_stat2', efotable_main.filter(regex='^pbs|xpehh').sum(1))

    # Order results
    efotable_main = efotable_main.sort_values(['sum_stat2', 'Phenotype'], ascending=[False, True])

    # Write main table
    writer = pd.ExcelWriter(f'metanalysis/efotable_{phenotype_dataset}.xlsx', engine='xlsxwriter')
    efotable_main.to_excel(writer, index=False)
    writer.sheets['Sheet1'].set_column(0, 0, 70)
    writer.sheets['Sheet1'].autofilter(0, 0, efotable_main.shape[0]-1, efotable_main.shape[1]-1)
    writer.save()

