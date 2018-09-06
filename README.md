# Introdoction

This is pipeline written in **Python** programming language. It is used for the analysis of *Mycobacterium tuberculosis* Methylation profile data generated using PacBio Sequencing system. The input data for the pipeline were generated using SMRT analysis tool from PacBio sequencing data.

# Additional data

- [Reference Genome H37Rv](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz) in fasta format
- [ProteinTable](https://www.ncbi.nlm.nih.gov/genome/proteins/166?genome_assembly_id=159857)
-
- Fasta sequences of genese need to be analysed


# Generating Initial Data

## Generating methylation profile

Methylation profile was generate using SMRT Analysis tool (Base modification and Motif Analysis option). H37Rv was used as reference genome for  *Mycobacterium tuberculosis*.
**Identify basemods** in **ADVANCED ANALYSIS PARAMETERS** option of the analysis **m6A,m4C** was added.

## Generating Genome based on mapped reads

Resequencing option was selected with H37Rv as reference genome. **quiver** was selected from **Algorithm** dropdown option in **ADVANCED ANALYSIS PARAMETERS**.   

# Objective

- Studying m6A Methylation profile
- Phylogeny analysis of the provided sample
- Selected gene analysis
- SNPs analysis

# External Tools Require by Script
- [MUSCLE](https://www.drive5.com/muscle/)
- [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/)
- [BLAT](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)

# External  Python Modules  Required

All the package can be installed using pip command

- Biopython
- ruffus
- pandas
- ete3
- seaborn
- click



# Running Steps

- Download methylation profile identfied using SMRT tools
- Download genome generated (fastq file) using resequencing analysis in SMRT tools and convert that in consensus.fasta
- Put all the gene sequences need to be analysed in separated file (*.ffn) in __selected_genes__ folder and put __selected_genes__ folder in __common_files__ folder.
- Download an outgroup genome (NC_015848.a.fna) and put it in the __common_files__. In our analysis [NC_015848](https://www.ncbi.nlm.nih.gov/genome/?term=NC_015848) is used.
- Protein table was uploaded to [Some Server]() to indetify the operons. Operons file was downloaded and saved as __Operon_locations.tsv__ in __common_files__ folder.
- Outgroup genome was file __NC_015848.1.fna__ file saved in __common_files__
- Each sample folder must have two folders, 5_BASE_MODIFICATIONS and Analysis.5_BASE_MODIFICATIONS should contain modifications.gff.gz and Analysis folder should contain consensus.fasta. All the sample folders were transferred to __raw_data__ folder
- __Results__ folder need to be created.


# Folders organisation
```
Project_folder
	|_raw_data
	|_Results [Folder]
	|_common_folder
	|_pipeline.py
```

# Running the script

python pipline.py


# System tested

Ubuntu 16.04 64bit
