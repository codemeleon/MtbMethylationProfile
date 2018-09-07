# Introduction

This is pipeline written in **Python** programming language. It is used for the analysis of *Mycobacterium tuberculosis* Methylation profile data generated using PacBio Sequencing system. The input data for the pipeline were generated using SMRT analysis tool from PacBio sequencing data.

# Objectives

- Studying m6A Methylation profile
- Phylogeny analysis of the provided sample
- Selected gene analysis
- SNPs analysis

# Additional data

- [Reference Genome H37Rv](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz) in fasta format
- Outgroup sequence [NC_015848](https://www.ncbi.nlm.nih.gov/genome/?term=NC_015848)
- [ProteinTable](https://www.ncbi.nlm.nih.gov/genome/proteins/166?genome_assembly_id=159857)

# External Tools Require by Script
- [MUSCLE](https://www.drive5.com/muscle/)
- [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/)
- [BLAT](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [SMRT Link](https://www.pacb.com/support/software-downloads/)


# External  Python Modules  Required

All the package can be installed using pip command

- Biopython
- ruffus
- pandas
- ete3
- seaborn
- click

# Generating Initial Data

## Generating methylation profile

Methylation profile was generate using SMRT Link tool (Base modification and Motif Analysis option). H37Rv was used as reference genome for  *Mycobacterium tuberculosis*.
**Identify basemods** in **ADVANCED ANALYSIS PARAMETERS** option of the analysis, **m6A,m4C** was added.

## Generating Genome based on mapped reads

Resequencing from Analysis option was selected with H37Rv as reference genome. **quiver** was selected from **Algorithm** dropdown option in **ADVANCED ANALYSIS PARAMETERS**.

# Running Steps

- Methylation profile gff was downloaded from SMRT Link analysis tool and compressed in modifications.gff.gz
- Genome generated by mapping reads to reference genome was downloaded from SMRT link analysis tool in fastq format and convert to consensus.fasta.
- __Operon_locations.tsv__ (operon file) was generated by uploading H37Rv reference fasta to [Sortberry](http://www.softberry.com/berry.phtml?topic=fgenesb&group=programs&subgroup=gfindb). _Mycobacterium tuberculosis_ H37Rv was selected as closes organism.
- Gene which need to be  studied, sequences were kept in individual file (.ffn) in __selected_genes__ folder. __selected_genes__ folder was transferred into __common_files__ folder.
- Outgroup genome fasta (.fna), reference genome fasta (.fna), __Operon_locations.tsv__ and __selected_genes__ folder were move __common_files__.
- Each sample folder must have two folders, 5_BASE_MODIFICATIONS and Analysis.5_BASE_MODIFICATIONS should contain modifications.gff.gz and Analysis folder should contain consensus.fasta. All the sample folders were transferred to __raw_data__ folder
- __Results__ folder need to be created.
- 'python pipline.py' was executed


# Folders organisation
```
Project_folder
	|_raw_data
	|		|_sample1
  |   |		|_5_BASE_MODIFICATIONS
  |   |   |		|_modifications.gff.gz
  |   |   |_Analysis
  |   |   		|_consensus.fasta
  |   |_sample2
	|   |		|_5_BASE_MODIFICATIONS
  |   |   |		|_modifications.gff.gz
  |   |   |_Analysis
  |   |   		|_consensus.fasta
  |   :
  |		:
	|_Results [Folder]
	|_common_folder
	|_pipeline.py
```

# Running the script

python pipline.py


# System tested

Ubuntu 16.04 64bit
