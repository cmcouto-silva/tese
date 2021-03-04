# Globals ---------------------------------------------------------------------

datasets = ['dataset_S1', 'dataset_S2']
container: 'docker://cmcoutosilva/tese:latest'

# Rules -----------------------------------------------------------------------	

rule all:
	input:
		# EHH-based Analysis
		expand("{dataset}/rehh/bin/{bin}.RDS", dataset=datasets, bin=['haplohh_list','scanhh_list']),
		expand("{dataset}/rehh/results/{out}.csv", dataset=datasets, out=['ihs_amz','xpehh_amz_mes','xpehh_amz_eas']),
		# PBS Analysis
		expand('{dataset}/pbs/figures/Amazonia_Mesoamerica_EastAsia{mode}.png', dataset=datasets, mode=['', '_mean']), # rm [:1]
		[f'{dataset}/pbs/results/{file[0]}/Amazonia_Mesoamerica_EastAsia.{file[1]}' for 
		file in [('RDS','Rds'),('XLS','xls')] for dataset in datasets # rm [:1]
		],
 		# Binary files
 		expand("{dataset}/results_bin/{mode}/{method}.feather",
		 dataset=datasets, mode=['persnp','pergene'],
		 method=['ihs','xpehh','xpehh_contrast','pbs','pbs_windowed','pbs_windowed_contrast']
		 ),
		# Phenotype Annotation
		expand("{dataset}/pheno_annot/results_annot{pat}.pickle", pat=['','_contrast'], dataset=datasets),
		expand("{dataset}/pheno_annot/results/{db}_{mode}_{method}.csv", dataset=datasets, db=['gwas','ensembl'], mode=['persnp','pergene'], method=['ihs','xpehh','pbs_windowed']),
		# Enrichr
		expand("{dataset}/enrichr/results/{method}_{mode}_pvalue{p}.txt", method=['pbs_windowed','xpehh','ihs'], mode=['persnp','pergene'],
		dataset=datasets, p=['0.999', '0.995', '0.99', '0.95']),
		# GTEx
		expand("{dataset}/gtex_annot/eqtl_persnp_{out}.pickle", dataset=datasets, out=['df','dict']),
		expand("{dataset}/gtex_annot/eqtl_persnp_{method}.{ext}", dataset=datasets, method=['pbsw','xpehh'], ext=['json','csv']),
		# WebGestalt
		expand("{dataset}/webgestalt/ORA/{mode}/{method}/pvalue_{p}/pvalue_{p}.RDS", dataset=datasets,
		mode=['persnp','pergene'], method=['pbs_windowed', 'xpehh', 'ihs'], p=['0.999', '0.995', '0.99', '0.95']),
		expand("{dataset}/webgestalt/GSEA/gsea.RDS", dataset=datasets),
		# GOATOOLS
		expand("{dataset}/goatools/go_{cat}.pickle", dataset=datasets, cat=['all','slim','immune']),
		expand("{dataset}/goatools/results/go_{cat}_{mode}_{method}_{p}.csv", dataset=datasets,
		cat=['all','slim','immune'], mode=['persnp','pergene'], method=['pbs_windowed','xpehh','ihs'], p=['0.999', '0.995', '0.99', '0.95']),
		# Metanalysis
		expand("metanalysis/efotable_{pmid}{mode}.xlsx", pmid=['by_pmid_',''], mode=['ensembl','gwas']),
		# Results
		expand([f'figures/{{ext}}/{fig}.{{ext}}' for fig in ["ds1_density", "ds1_ihs_pergene", "ds1_ihs_persnp", "ds1_intersection", "ds1_pbs_eqtl_score", "ds1_pbs_eqtl", 
		"ds1_pbsw_both", "ds1_pbsw_pergene", "ds1_pbsw_persnp", "ds1_pergene_density", "ds1_persnp_density", "ds1_xpehh_both", "ds1_xpehh_eqtl_score", "ds1_xpehh_eqtl", "ds1_xpehh_pergene", 
		"ds1_xpehh_persnp", "ds2_density", "ds2_ihs_pergene", "ds2_ihs_persnp", "ds2_intersection", "ds2_pbs_eqtl_score", "ds2_pbs_eqtl", "ds2_pbsw_both", "ds2_pbsw_pergene", "ds2_pbsw_persnp", 
		"ds2_pergene_density", "ds2_persnp_density", "ds2_xpehh_both", "ds2_xpehh_eqtl_score", "ds2_xpehh_eqtl", "ds2_xpehh_pergene", "ds2_xpehh_persnp", "overall_eqtl_count_split", 
		"overall_eqtl_count", "overall_eqtl_score", "pbs_intersection", "pbs_pergene_intersection.995", "pbs_pergene_intersection.999", "pbs_persnp_intersection.995", "pbs_persnp_intersection.999"]],
		ext=['pdf', 'png']),
		[f'tables/tex/{file}.tex' for file in ["ds1_ensembl_pergene", "ds1_ensembl_persnp", "ds1_fumagwas_ihs_persnp", "ds1_go_all_persnp", "ds1_gwas_pergene", "ds1_gwas_persnp",
		"ds1_ihs_pergene", "ds1_ihs_persnp", "ds1_natives", "ds1_pbsw_pergene", "ds1_pbsw_persnp", "ds1_selgroups", "ds1_xpehh_pergene", "ds1_xpehh_persnp", "ds2_ensembl_pergene",
		"ds2_ensembl_persnp", "ds2_fumagwas_pbs_pergene", "ds2_fumagwas_pbs_persnp", "ds2_go_all_pergene", "ds2_go_immune_pergene", "ds2_gwas_pergene", "ds2_gwas_persnp", "ds2_ihs_pergene",
		"ds2_ihs_persnp", "ds2_kegg_ora_persnp", "ds2_natives", "ds2_pbsw_pergene", "ds2_pbsw_persnp", "ds2_selgroups", "ds2_xpehh_pergene", "ds2_xpehh_persnp", "gsea", "metanalysis_ensembl",
		"metanalysis", "pbsw_pergene"]],
		[f'tables/xlsx/{file}' for file in ["ds1_pheno_tables.xlsx", "ds2_pheno_tables.xlsx", "enrichment.xls", "eqtl_score.xls", "ihs.xls", "pbs.xls", "populations.xls", "xpehh.xls"]]


rule setup_env:
	conda:
		"env/dependencies.yaml"
	script:
		"env/pkgs.R"


rule rehh:
	input:
		script = "{dataset}/rehh/scripts/rehh.R",
		data = "{dataset}/data/phased_adj_alleles/{dataset}.vcf.gz"
	conda:
		"env/dependencies.yaml"
	output:
		expand("{{dataset}}/rehh/bin/{bin}.RDS", bin=['haplohh_list','scanhh_list']),
		expand("{{dataset}}/rehh/results/{out}.csv", out=['ihs_amz','xpehh_amz_mes','xpehh_amz_eas'])
	shell:
		"Rscript {input.script}"


rule pbs:
	input:
		script = '{dataset}/pbs/scripts/pbs_script.R',
		scripts = expand('{{dataset}}/pbs/scripts/pbs_{script}.R', script=['setup','pipeline']) 
	output:
		expand('{{dataset}}/pbs/figures/Amazonia_Mesoamerica_EastAsia{mode}.png', mode=['', '_mean']),
		[f'{{dataset}}/pbs/results/{file[0]}/Amazonia_Mesoamerica_EastAsia.{file[1]}' for file in [('RDS','Rds'),('XLS','xls')]]
	conda:
		"env/dependencies.yaml"
	shell:
		"Rscript {input.script}"


rule results_bin:
	input:
		script = '{dataset}/results_bin/convert2feather.R',
		# iHS & XP-EHH
		rehh = expand("{{dataset}}/rehh/results/{out}.csv", out=['ihs_amz','xpehh_amz_mes','xpehh_amz_eas']),
		# PBS
		pbs = [f'{{dataset}}/pbs/results/{file[0]}/Amazonia_Mesoamerica_EastAsia.{file[1]}' for file in [('RDS','Rds'),('XLS','xls')]]
	output:
		expand("{{dataset}}/results_bin/{mode}/{method}.feather", mode=['persnp','pergene'], method=['ihs','xpehh','xpehh_contrast','pbs','pbs_windowed','pbs_windowed_contrast'])
	conda:
		"env/dependencies.yaml"
	shell:
		"Rscript {input.script}"


rule enrichr:
	input: 
		script = "{dataset}/enrichr/scripts/enrichr.py",
		bin_files = expand("{{dataset}}/results_bin/{mode}/{method}.feather", mode=['persnp','pergene'], method=['ihs','xpehh','xpehh_contrast','pbs','pbs_windowed','pbs_windowed_contrast'])
	output:
		expand("{{dataset}}/enrichr/results/{method}_{mode}_pvalue{p}.txt", method=['pbs_windowed','xpehh','ihs'], mode=['persnp','pergene'], p=['0.999', '0.995', '0.99', '0.95']),
	conda:
		"env/dependencies.yaml"
	shell:
		"python {input.script}"


rule get_eQTL:
	input: 
		script = "{dataset}/gtex_annot/scripts/get_snp_eQTLs.py",
		bin_files = expand("{{dataset}}/results_bin/{mode}/{method}.feather", mode=['persnp','pergene'], method=['ihs','xpehh','xpehh_contrast','pbs','pbs_windowed','pbs_windowed_contrast'])
	output:
		expand("{{dataset}}/gtex_annot/eqtl_persnp_{out}.pickle", out=['df','dict']),
		expand("{{dataset}}/gtex_annot/eqtl_persnp_{method}.{ext}", method=['pbsw','xpehh'], ext=['json','csv'])
	conda:
		"env/dependencies.yaml"
	shell:
		"python {input.script}"


rule webgestalt:
	input: 
		script = "{dataset}/webgestalt/scripts/webgestalt.R",
		bin_files = expand("{{dataset}}/results_bin/{mode}/{method}.feather", mode=['persnp','pergene'], method=['ihs','xpehh','xpehh_contrast','pbs','pbs_windowed','pbs_windowed_contrast'])
	output:
		expand("{{dataset}}/webgestalt/ORA/{mode}/{method}/pvalue_{p}/pvalue_{p}.RDS", 
		mode=['persnp','pergene'], method=['pbs_windowed', 'xpehh', 'ihs'], p=['0.999', '0.995', '0.99', '0.95']),
		"{dataset}/webgestalt/GSEA/gsea.RDS"
	conda:
		"env/dependencies.yaml"
	shell:
		"Rscript {input.script}"


rule goatools:
	input: 
		script = "{dataset}/goatools/scripts/goatools_script.py",
		bin_files = expand("{{dataset}}/results_bin/{mode}/{method}.feather", mode=['persnp','pergene'], method=['ihs','xpehh','xpehh_contrast','pbs','pbs_windowed','pbs_windowed_contrast'])
	output:
		expand("{{dataset}}/goatools/go_{cat}.pickle", cat=['all','slim','immune']),
		expand("{{dataset}}/goatools/results/go_{cat}_{mode}_{method}_{p}.csv", cat=['all','slim','immune'], mode=['persnp','pergene'], method=['pbs_windowed','xpehh','ihs'], p=['0.999', '0.995', '0.99', '0.95'])
	conda:
		"env/dependencies.yaml"
	shell:
		"python {input.script}"


rule metanalysis:
	input:
		script = "metanalysis/metanalysis_script.py",
		bin_files = expand("{dataset}/results_bin/{mode}/{method}.feather", dataset=datasets,
		mode=['persnp','pergene'], method=['ihs','xpehh','xpehh_contrast','pbs','pbs_windowed','pbs_windowed_contrast'])
	output:
		expand("metanalysis/efotable_{pmid}{mode}.xlsx", pmid=['by_pmid_',''], mode=['ensembl','gwas'])
	conda:
		"env/dependencies.yaml"
	shell:
		#"jupyter nbconvert --execute --inplace --to notebook {input.script}"
		"python {input.script}"


rule pheno_annot:
	input:
		script = expand("{{dataset}}/pheno_annot/{file}.py", file=['pheno_annot','pheno_annot_contrast', 'make_tables']),
		bin_files = expand("{dataset}/results_bin/{mode}/{method}.feather", dataset=datasets,
		mode=['persnp','pergene'], method=['ihs','xpehh','xpehh_contrast','pbs','pbs_windowed','pbs_windowed_contrast'])
	output:
		expand("{{dataset}}/pheno_annot/results_annot{pat}.pickle", pat=['','_contrast']),
		expand("{{dataset}}/pheno_annot/results/{db}_{mode}_{method}.csv", db=['gwas','ensembl'], mode=['persnp','pergene'], method=['ihs','xpehh','pbs_windowed'])
	conda:
		"env/dependencies.yaml"
	shell:
		"""
		python {input.script[0]}
		python {input.script[1]}
		python {input.script[2]}
		"""


rule summarise_all:
	input:
		# Script
		"scripts/summarise_results.R",
		[f'scripts/{file}' for file in ["density_plot.R", "efotrait2DT.R", "enrichment.R", "gene2efotrait.R", "get_intersection.R", "gtex_plot.R",
		"ihs.R", "manhattan.R", "pbs.R", "pheno_tables.R", "pop_tables.R", "summarise_results.R", "xpehh.R"]],
		# Phenotype Annotation
		expand("{dataset}/pheno_annot/results_annot{pat}.pickle", pat=['','_contrast'], dataset=datasets),
		expand("{dataset}/pheno_annot/results/{db}_{mode}_{method}.csv", dataset=datasets, db=['gwas','ensembl'], mode=['persnp','pergene'], method=['ihs','xpehh','pbs_windowed']),
		# Enrichr
		expand("{dataset}/enrichr/results/{method}_{mode}_pvalue{p}.txt", method=['pbs_windowed','xpehh','ihs'], mode=['persnp','pergene'],
		dataset=datasets, p=['0.999', '0.995', '0.99', '0.95']),
		# GTEx
		expand("{dataset}/gtex_annot/eqtl_persnp_{out}.pickle", dataset=datasets, out=['df','dict']),
		expand("{dataset}/gtex_annot/eqtl_persnp_{method}.{ext}", dataset=datasets, method=['pbsw','xpehh'], ext=['json','csv']),
		# WebGestalt
		expand("{dataset}/webgestalt/ORA/{mode}/{method}/pvalue_{p}/pvalue_{p}.RDS", dataset=datasets,
		mode=['persnp','pergene'], method=['pbs_windowed', 'xpehh', 'ihs'], p=['0.999', '0.995', '0.99', '0.95']),
		expand("{dataset}/webgestalt/GSEA/gsea.RDS", dataset=datasets),
		# GOATOOLS
		expand("{dataset}/goatools/go_{cat}.pickle", dataset=datasets, cat=['all','slim','immune']),
		expand("{dataset}/goatools/results/go_{cat}_{mode}_{method}_{p}.csv", dataset=datasets,
		cat=['all','slim','immune'], mode=['persnp','pergene'], method=['pbs_windowed','xpehh','ihs'], p=['0.999', '0.995', '0.99', '0.95']),
		# Metanalysis
		expand("metanalysis/efotable_{pmid}{mode}.xlsx", pmid=['by_pmid_',''], mode=['ensembl','gwas'])
	output:
		figures = expand([f'figures/{{ext}}/{fig}.{{ext}}' for fig in ["ds1_density", "ds1_ihs_pergene", "ds1_ihs_persnp", "ds1_intersection", "ds1_pbs_eqtl_score", "ds1_pbs_eqtl", 
		"ds1_pbsw_both", "ds1_pbsw_pergene", "ds1_pbsw_persnp", "ds1_pergene_density", "ds1_persnp_density", "ds1_xpehh_both", "ds1_xpehh_eqtl_score", "ds1_xpehh_eqtl", "ds1_xpehh_pergene", 
		"ds1_xpehh_persnp", "ds2_density", "ds2_ihs_pergene", "ds2_ihs_persnp", "ds2_intersection", "ds2_pbs_eqtl_score", "ds2_pbs_eqtl", "ds2_pbsw_both", "ds2_pbsw_pergene", "ds2_pbsw_persnp", 
		"ds2_pergene_density", "ds2_persnp_density", "ds2_xpehh_both", "ds2_xpehh_eqtl_score", "ds2_xpehh_eqtl", "ds2_xpehh_pergene", "ds2_xpehh_persnp", "overall_eqtl_count_split", 
		"overall_eqtl_count", "overall_eqtl_score", "pbs_intersection", "pbs_pergene_intersection.995", "pbs_pergene_intersection.999", "pbs_persnp_intersection.995", "pbs_persnp_intersection.999"]],
		ext=['pdf', 'png']),
		tables_tex = [f'tables/tex/{file}.tex' for file in ["ds1_ensembl_pergene", "ds1_ensembl_persnp", "ds1_fumagwas_ihs_persnp", "ds1_go_all_persnp", "ds1_gwas_pergene", "ds1_gwas_persnp",
		"ds1_ihs_pergene", "ds1_ihs_persnp", "ds1_natives", "ds1_pbsw_pergene", "ds1_pbsw_persnp", "ds1_selgroups", "ds1_xpehh_pergene", "ds1_xpehh_persnp", "ds2_ensembl_pergene",
		"ds2_ensembl_persnp", "ds2_fumagwas_pbs_pergene", "ds2_fumagwas_pbs_persnp", "ds2_go_all_pergene", "ds2_go_immune_pergene", "ds2_gwas_pergene", "ds2_gwas_persnp", "ds2_ihs_pergene",
		"ds2_ihs_persnp", "ds2_kegg_ora_persnp", "ds2_natives", "ds2_pbsw_pergene", "ds2_pbsw_persnp", "ds2_selgroups", "ds2_xpehh_pergene", "ds2_xpehh_persnp", "gsea", "metanalysis_ensembl",
		"metanalysis", "pbsw_pergene"]],
		tables_xlsx = [f'tables/xlsx/{file}' for file in ["ds1_pheno_tables.xlsx", "ds2_pheno_tables.xlsx", "enrichment.xls", "eqtl_score.xls", "ihs.xls", "pbs.xls", "populations.xls", "xpehh.xls"]]
	conda:
			"env/dependencies.yaml"
	script:
			"scripts/summarise_results.R"

