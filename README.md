# Tese: Influência da seleção natural em populações nativas de diferentes ecorregiões americanas

### Descrição

Este repositório contém os códigos necessários para reprodução da minha tese.
Para compilar o arquivo da tese em pdf, ou visualizar o pdf já compilado, favor verificar o repositório https://github.com/cmcouto-silva/tese_latex.

### Requerimentos

Para o gerenciamento do fluxo de trabalho (workflow) foi utilizado [Snakemake](https://snakemake.readthedocs.io/en/stable/).

Como input para o workflow, espera-se os seguintes arquivos:

```
.
├── dataset_S1
│   └── data
│       ├── annot
│       │   ├── annovar_refGene.annot
│       │   ├── annovar_refGene_SNP.annot
│       │   ├── gene_dist.annot
│       │   ├── gene_dist_collapsed.annot
│       │   └── gene_keys.annot
│       ├── phased
│       │   ├── dataset_S1.bed
│       │   ├── dataset_S1.bim
│       │   └── dataset_S1.fam
│       ├── phased_adj_alleles
│       │   ├── dataset_S1.h5
│       │   ├── dataset_S1.vcf.gz
│       │   └── dataset_S1.vcf.gz.tbi
│       └── populations
│           ├── populations.csv
│           └── populations.xls
└── dataset_S2
    └── data
        ├── annot
        │   ├── annovar_refGene.annot
        │   ├── annovar_refGene_SNP.annot
        │   ├── gene_dist.annot
        │   ├── gene_dist_collapsed.annot
        │   └── gene_keys.annot
        ├── phased
        │   ├── dataset_S2.bed
        │   ├── dataset_S2.bim
        │   └── dataset_S2.fam
        ├── phased_adj_alleles
        │   ├── dataset_S2.h5
        │   ├── dataset_S2.vcf.gz
        │   └── dataset_S2.vcf.gz.tbi
        └── populations
            └── populations.csv
```

Dado os múltiplos softwares e bibliotecas necessárias para gerar todos os resultados da tese, a maneira mais prática de lidar com estas dependências sem alterar diretamente o sistema operacional é através de *containers*. 
Aqui, forneço duas opções para lidar com as dependências: 1) via containerização com [Docker](https://hub.docker.com/repository/docker/cmcoutosilva/tese) e 2) via ambiente virtual com [Conda](https://docs.conda.io/en/latest/). 

Para instalar o Snakemake, recomenda-se instalar primeiramente o **Miniconda** ([vide aqui](https://docs.conda.io/en/latest/miniconda.html) como instalar o Miniconda).

Uma vez instalado o Miniconda, crie um novo ambiente virtual para Snakemake com o seguinte comando:

```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake singularity -y
```

### Executando o Workflow

Para cada dataset, executamos as análises de seleção baseadas no desequilíbrio de ligação (REHH) e diferenciação populacional (Fst - PBS) e salvamos os resultados em arquivos binários (results_bin), que por sua vez são utilizados nas demais análises. Para cada uma das análises, há um script na pasta correspondente ao dataset, e uma vez que todas as análises estejam completas, executa-se a regra "sumarise_all", que compila todos estes resultados e gera as figuras e tabelas apresentadas na tese e disponíveis no [repositório de dados Mendeley](http://dx.doi.org/10.17632/gztff7wmjt.1) (DOI:
10.17632/gztff7wmjt.1).

O esquema das análises está representado na figura abaixo e é orquestrado via Snakemake (vide [Snakefile](https://github.com/cmcouto-silva/tese/blob/main/Snakefile)), cuja instalação das dependências, gerenciamento de inputs e outputs, bem como paralelização das computações é configurada automaticamente.

![directed acyclic graph](dag.svg)

Para executar o workflow, ative o ambiente virtual criado para rodar o Snakemake:

    conda activate snakemake

Em seguida, podemos executar o workflow via Docker (necessário a instalação do Docker):

    snakemake -j1 --use-singularity

ou via Conda:

```
snakemake -j1 --use-conda --until setup_env
snakemake -j1 --use-conda
```

**Nota:** Recomenda-se que substitua o parâmetro `j1` por um número maior de núcleos, como, por exemplo, `-j11` (todos os testes que realizei foi com este parâmetro). Ao possibilitar o uso de múltiplas _threads_, Snakemake configura automaticamente o workflow para ser computado com pararelização.

Para saber quantos núcleos tem disponível no seu computador, no ambiente virtual do Snakemake, execute o seguinte comando:

    python -c "import multiprocessing; print(multiprocessing.cpu_count())"

### Resultados

Todas as figuras e tabelas resultantes deste workflow estão disponíveis no [repositório de dados Mendeley](http://dx.doi.org/10.17632/gztff7wmjt.1).
Vale ressaltar, contudo,
que os dados de entrada (_e.g._ Plink, VCF) possuem limitações éticas e por esse motivo não estão inclusos no repositório.

### Observações

Se for utilizar o [Docker container](https://hub.docker.com/repository/docker/cmcoutosilva/tese), certifique-se de [instalar](https://docs.docker.com/get-docker/) o Docker primeiro.

Se for utilizar o arquivo de configuração para o ambiente conda, note que se faz necessário executar o comando `snakemake -j1 --use-conda --until setup_env` uma primeira e única vez, para a correta instalação de todas as dependências do workflow, visto que há pacotes do R não inclusos nos repositórios do Conda.
