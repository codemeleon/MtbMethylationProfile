#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
from pylab import *
import seaborn as sns
from ruffus import *
from ruffus.combinatorics import *
from Bio import SeqIO, AlignIO
import pandas as pd
import click
from glob import glob
from os import path, system
import os
from ete3 import Tree, TreeStyle

##
#mugsy -p mugsy 12366_10.fasta 12366_11.fasta 12366_12.fasta 12366_13.fasta 12366_1.fasta 12366_2.fasta 12366_3.fasta 12366_7.fasta 12366_8.fasta 12366_9.fasta 9099_1.fasta 9099_2.fasta 9099_3.fasta 9099_4.fasta


## Clobal variables
# prot_file = None
# genf = None

def context(txt):
    return txt.split(";")[1].split("=")[1]

def frac(txt):
    return txt.split(";")[3].split("=")[1]


def seq_ext(seqfile, seqid, start, end, strand):
    seq = SeqIO.to_dict(SeqIO.parse(seqfile, "fasta"))
    seq = seq[seqid].seq[start:end]
    if strand == '-':
        seq = seq.reverse_complement()
    return seq

def msa(infile, outfile):
    system("muscle -in %s -out %s" %(infile, outfile))


def operontab(inputfiles, outputfile, upstrm):
    data = pd.read_table(inputfiles[0],
                         compression="gzip", header=None, comment="#")
    data = data[data[2] == "m6A"]
    data["context"] = list(map(context, data[8]))
    data["frac"] = list(map(frac, data[8]))
    data = data.drop([0, 1, 2, 3, 5, 7], axis=1)
    data = data.rename(columns={4: 'position',
                                6: 'strand',
                                8: 'info'})
    data['operon'] = ''
    data['operon_strand'] = ''
    operon = pd.read_table(inputfiles[1], header=None)
    range_operon = {'start_min': operon.groupby([1])[7].min(),
                    'end_max': operon.groupby([1])[9].max()}
    range_operon = pd.DataFrame.from_dict(range_operon)
    strands = pd.pivot_table(operon, index=[1], columns=[5], aggfunc=len)[0]
    strands[pd.isnull(strands)] = 0
    range_operon['strand'] = '-'
    range_operon.loc[strands['+'] > 0, 'strand'] = '+'
    range_operon = range_operon.reset_index()
    range_operon_shape = range_operon.shape
    for i, row in range_operon.iterrows():
        tstart, tend = 0, 0
        if row['strand'] == '+':
            tend = row['start_min']
            if i == 0:
                if row['start_min'] < upstrm:
                    tstart = 0
                else:
                    tstart = row['start_min']-upstrm
            else:
                if row['start_min'] - range_operon.ix[i-1]['end_max'] < upstrm:
                    tstart = range_operon.ix[i-1]['end_max']
                else:
                    tstart = row['start_min'] - upstrm
        else:
            tstart = row['end_max'] + 1
            if i == range_operon_shape[0] - 1:
                tend = row['end_max'] + upstrm + 1
            else:
                if ((range_operon.ix[i+1]['start_min'] - row['end_max']) <
                        upstrm):
                    tend = range_operon.ix[i + 1]['start_min']
                else:
                    tend = row['end_max'] + upstrm + 1

        selected_indexes = data[(data['position'] >= tstart) &
                                (data['position'] < tend)].index
        for ix in selected_indexes:
            data.loc[ix, 'operon'] = str(row[1])
            data.loc[ix, 'operon_strand'] = row['strand']

    data.to_csv(outputfile, sep="\t", index=False)




def genictab(infiles, outfile):

    data = pd.read_table(infiles[0], compression="gzip", header=None, comment="#")
    data = data[data[2] == "m6A"]
    data["context"] = list(map(context, data[8]))
    data["frac"] = list(map(frac, data[8]))
    del data[8]
    protein_table = pd.read_table(infiles[1])
    protein_table["Locus"] = list(map(str.lower, protein_table["Locus"]))
    genf = infiles[2]
    if genf:
        genes = list(pd.read_table(genf, header=None)[0])
        genes = list(map(str.lower, genes))
        #print(genes)
        protein_table = protein_table[protein_table["Locus"].isin(genes)]
    columns = ["gene", "start", "end", "mod_site", "gene_strand", "mod_strand",
               "frac_mod", "context"]
    geneinfo = {"gene": [], "start": [], "end": [], "mod_site": [],
                "gene_strand": [], "mod_strand": [], "frac_mod": [],
                "context": []}
    for i, rowi in protein_table.iterrows():
        ttable = data[(data[3] >= rowi[2]) & (data[3] < rowi[3])]
        for j, rowj in ttable.iterrows():
            geneinfo["gene"].append(rowi["Locus"])
            geneinfo["start"].append(rowi["Start"])
            geneinfo["end"].append(rowi["Stop"])
            geneinfo["mod_site"].append(rowj[3])
            geneinfo["gene_strand"].append(rowi["Strand"])
            geneinfo["mod_strand"].append(rowj[6])
            geneinfo["frac_mod"].append(rowj["frac"])
            geneinfo["context"].append(rowj["context"])
    geneinfo = pd.DataFrame.from_dict(geneinfo)
    geneinfo = geneinfo[columns]
    geneinfo.to_csv(outfile, index=False, sep="\t")



def intergenictab(infiles, outfile):

    data = pd.read_table(infiles[0],
                         compression="gzip", header=None, comment="#")
    data = data[data[2] == "m6A"]
    data["context"] = list(map(context, data[8]))
    data["frac"] = list(map(frac, data[8]))
    data = data.drop([0, 1, 2, 4, 5, 7, 8], axis=1)
    data['leftgene'] = ['-'] * data.shape[0]
    data['leftgene_strand'] = ['-']*data.shape[0]
    data['rightgene'] = ['-'] * data.shape[0]
    data['rightgene_strand'] = ['-'] * data.shape[0]
    protein_table = pd.read_table(infiles[1],
                                  usecols=["Start", "Stop", "Strand", "Locus"])
    drop_rows = []
    for i, row in protein_table.iterrows():
        drop_rows += list(data.loc[((data[3] < row['Stop']) &
                                    (data[3] >= row['Start'])), :].index)

    data = data.drop(drop_rows, axis=0)
    for i, row in protein_table.iterrows():
        if i == 0:
            data.loc[(data[3] < row['Start']),
                     ['leftgene', 'leftgene_strand',
                      'rightgene', 'rightgene_strand']
                     ] = ['-', '-', row['Locus'], row['Strand']]
        else:
            data.loc[((data[3] >= protein_table.loc[i-1, 'Stop']) &
                      (data[3] < row['Start'])),
                     ['leftgene', 'leftgene_strand', 'rightgene',
                      'rightgene_strand']
                     ] = [protein_table.loc[i-1, 'Locus'],
                          protein_table.loc[i-1, 'Strand'],
                          row['Locus'],
                          row['Strand']]
            if i == protein_table.shape[1] - 1:
                data.loc[(data[3] > row['Stop']),
                         ['leftgene', 'leftgene_strand',
                          'rightgene', 'rightgene_strand']
                         ] = [row['Locus'],
                              row['Strand'],
                              '-',
                              '-']
    data.to_csv(outfile, index=False, sep="\t")




# @click.command()
# @click.option("-sd", help="Sample's directory", type=str, default=None,
              # show_default=True)
# @click.option("-pf", help="Protein File", type=str, default=None,
#               show_default=True)


sd = "../"
pf = "../ProteinTable166_159857.txt"
gf = "../genes"
gs = "all" #gene, 5utr
of = "../Operon_locations.tsv"
upstrm = 100
rgd="../selected_gene"
outgroup = "../NC_015848.1.fna" #  Name should be the same as sequence name
def run(sd, pf, gf, of, upstrm, rgd, outgroup):
    """
    This source code analyses data from victor.
    """
    if not sd:
        exit("Sample data folder not given. Exiting....")
    if not path.isdir(sd):
        exit("Given sample directory path is not folder. Exiting . . . .")
    if not pf:
        exit("Protein tab file not given. Exiting.....")
    if not path.isfile(pf):
        exit("Give protein tab file path doesn't exist or it is not a file"
             " Exiting....... ")
    global prot_file, genes
    # prot_file = pf
    # gen = gf
    files = glob("%s/*/5_BASE_MODIFICATIONS/modifications.gff.gz" % sd)
    if not files:
        exit("Given sample diectroy doesn't have valid data. Exiting....")
    expt = "anmol"
    @mkdir("Results")
    @transform(files, formatter(
        ".+/(?P<filebase>\w+)/5_BASE_MODIFICATIONS/modifications.gff.gz"),
        add_inputs(pf, gf),
    "Results/{filebase[0]}.genictab")
    def genictab_run(inputfiles, outputfiles):
        genictab(inputfiles, outputfiles)

    @follows(genictab_run)
    @transform(files, formatter(
        ".+/(?P<filebase>\w+)/5_BASE_MODIFICATIONS/modifications.gff.gz"),
        add_inputs(pf),
    "Results/{filebase[0]}.intergenictab")
    def intergenictab_run(inputfiles, outputfiles):
        intergenictab(inputfiles, outputfiles)

    @follows(genictab_run)
    @transform(files, formatter(
        ".+/(?P<filebase>\w+)/5_BASE_MODIFICATIONS/modifications.gff.gz"),
        add_inputs(of),
    "Results/{filebase[0]}.operontab", upstrm)
    def operontab_run(inputfiles, outputfiles, upstrm):
        operontab(inputfiles, outputfiles, upstrm)

    gene_files = glob("%s/*.ffn" % rgd)
    consensus_fasta = glob("%s/*/Analysis/consensus.fasta" % sd)
    consensus_fasta.append(outgroup)
    if gene_files:
        @follows(genictab_run)
        @product(consensus_fasta, formatter(
            ".+/(?P<basename>\w+)/Analysis/consensus.fasta"),
            gene_files, formatter(".+/(?P<basename>\w+).ffn"),
            "Results/{basename[0][0]}__{basename[1][0]}.psl")
        def blat_run(inputfiles, outputfile):
            system("blat -noHead %s %s %s" %(inputfiles[0], inputfiles[1],
                                            outputfile))
        @collate(blat_run, formatter(".+/*__(?P<basename>\w+).psl$"),
                 "Results/{basename[0]}.ffn")
        def gene_seq_groupping(inputfiles, outfile):
            ofile = open(outfile,"w")
            gene_file_path = "%s/%s" %(rgd, path.split(outfile)[1])
            for rec in SeqIO.parse(gene_file_path,"fasta"):
                ofile.write(">%s\n%s\n" %(rec.id, rec.seq))
            for inf in inputfiles:
                sample_name = path.split(inf)[1].split("__")[0]
                samp_consensus = "%s/%s/Analysis/consensus.fasta" % (sd,sample_name)
                df = pd.read_table(inf, header=None)
                df = df.ix[0].values.tolist()
                gen_seq = seq_ext(samp_consensus, df[13], df[15], df[16], df[8])
                ofile.write(">%s\n%s\n" %(sample_name, gen_seq))
                # print(df, "anmol", gen_seq)
                # continue
            ofile.close()
        @transform(gene_seq_groupping, suffix(".ffn"), ".faa")
        def translate(infile, outfile):
            with open(outfile, "w") as fout:
                for rec in SeqIO.parse(infile, "fasta"):
                    fout.write(">%s\n%s\n" %(rec.id, rec.seq.translate()))

        @transform(gene_seq_groupping, suffix(".ffn"), "_aln.ffn")
        def msa_nuc(infile, outfile):
            msa(infile, outfile)

        @transform(translate, suffix(".faa"), "_aln.faa")
        def msa_aa(infile, outfile):
            msa(infile, outfile)
    @follows(operontab_run)
    @merge(consensus_fasta, "Results/consensus.fasta")
    def genome_groupping(inputfiles, outputfile):
        print(inputfiles)
        with open(outputfile, "w") as fout:
            for fl in inputfiles[:-1]:
                sample_name = fl.split("/")[-3]
                for rec in SeqIO.parse(fl, "fasta"):
                    fout.write(">%s\n%s\n" %(sample_name, rec.seq))
            for rec in SeqIO.parse(inputfiles[-1], 'fasta'):
                fout.write(">%s\n%s\n" %(rec.id, rec.seq))
    @transform(genome_groupping, suffix(".fasta"), "_mafft.fasta")
    def mafft(inputfile, outputfile):
        system("mafft --thread -1 --auto %s > %s" %(inputfile, outputfile))

    @transform(mafft, suffix(".fasta"), "_snp.fasta")
    def snpaln(inputfile, outputfile):
        alignment = {}
        for rec in AlignIO.read(inputfile, 'fasta'):
            alignment[rec.id] = list(str(rec.seq).upper())
        alignment = pd.DataFrame(alignment)
        alignment = alignment[alignment.apply(lambda x: len(set(x)),
                                              axis=1) != 1]
        with open(outputfile, "w") as fout:
            for col in alignment.columns:
                fout.write(">%s\n%s\n" %(col, ''.join(alignment[col])))
    @transform(snpaln, suffix("_snp.fasta"), ("_snp_trimmed.fasta",
                                             "Tree.phy", "_snp_count", "_tree.png"))
    def phylo_snp_count(inputfile, outputfiles):
        # http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0038969
        sequences = {}
        for rec in SeqIO.parse(inputfile, "fasta"):
            sequences[rec.id] = list(str(rec.seq))
        sequences = pd.DataFrame(sequences)
        x = sequences.apply(lambda x: list(x).count("-"), axis=1)
        sequences = sequences[x == 0]
        snp_count = open(outputfiles[2], "w")
        snp_count.write("Number of SNPs: %d\n" % len(sequences))
        minval = sequences.apply(lambda x: min(list(x.value_counts())), axis=1)
        snp_count.write("Present in only sample: %s" % sum(minval==1))
        snp_count.close()
        with open(outputfiles[0], "w") as fout:
            for k in sequences:
                fout.write(">%s\n%s\n" % (k, ''.join(sequences[k])))
        system("raxmlHPC-PTHREADS-SSE3 -s %s -m GTRCAT -f a -T 30 -n"
               " temp_tree -x 13 -p 13 -N 1000 -o NC_015848.1" % outputfiles[0])
        system("mv RAxML_bestTree.temp_tree %s" %(outputfiles[1]))
        system("rm RAxML_*")
        tree = Tree(outputfiles[1])
        tree.ladderize()
        ts = TreeStyle()
        _ = tree.render(outputfiles[3], tree_style=ts)

    @transform(snpaln, suffix(".fasta"), (".mat", ".png"))
    def distance_matrix(inputfile, outputfiles):
        alignment = {}
        for rec in SeqIO.parse(inputfile, 'fasta'):
            alignment[rec.id] = list(rec.seq)
        alignment = pd.DataFrame(alignment)
        ix = alignment.columns
        distances = {}
        for c1 in ix:
           distances[c1]=[]
           for c2 in ix:
               distances[c1].append(sum(alignment[c1] != alignment[c2]))
        distances = pd.DataFrame(distances)
        distances = distances[ix]
        distances.index = ix
        distances.to_csv(outputfiles[0], index=False, sep="\t")
        sns.clustermap(distances).savefig(outputfiles[1])

    @transform(files,formatter(
        ".+/(?P<filebase>\w+)/5_BASE_MODIFICATIONS/modifications.gff.gz"),
        add_inputs(mafft), "Results/{filebase[0]}.change_aln")
    def aln_mod_site(inputfiles, outputfile):
        data = pd.read_table(inputfiles[0],
                             compression="gzip", header=None, comment="#")
        data = data.loc[data[2] == "m6A", 4]
        seqid = path.split(outputfile)[1].split('.cha')[0]
        seq = ''
        for rec in SeqIO.parse(inputfiles[1],'fasta'):
            if rec.id == seqid:
                seq = str(rec.seq)
        aln_pos = {'orgpos':[], 'alnpos':[]}
        aln = 1
        for i, c in enumerate(seq,1):
            if c != '-':
                if aln in data:
                    aln_pos['orgpos'].append(aln)
                    aln_pos['alnpos'].append(i)
                aln += 1
        aln_pos = pd.DataFrame(aln_pos)
        aln_pos.to_csv(outputfile, index=False, sep="\t")

        # pass

    # TODO: Add Mafft, Mugsy and Mauve Features
    # GoodToKnow Mauve messes up with original files
    # def align_genome():

        # #TODO: use mafft
        # pass
    # def phylogeny():
        # #TODO: based on SNPs
        # pass


    # # TODO: Add check fro protein file
    pipeline_run(verbose=9, multiprocess=1)




if __name__ == '__main__':
    run(sd, pf, gf, of, upstrm,rgd, outgroup)
