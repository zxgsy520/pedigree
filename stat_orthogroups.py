#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "1.2.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.lstrip().rstrip('\n')

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith("gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)

    seq = ''
    for line in fp:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()

        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            if seq!='':
                yield seq.split('\n')
            seq = "%s\n" % line
            continue
        seq += line

    if seq!='':
        yield seq.split('\n')
    fp.close()


def get_sample(file):

    name = file.split('/')[-1]

    if '--' in name:
        name = name.split('--')[1].split('.bam')[0]
    else:
        name = name.split('.')[0:-1]
        name = '.'.join(name)

    return name.strip('.')


def stat_gene(files):

    data = {}

    for file in files:
        name = get_sample(file)
        n = 0
        for seqid, seq in read_fasta(file):
            if len(seq)<=2:
                continue
            n += 1
        data[name] = n

    return data


def read_head(heads):

    data = OrderedDict()
    n = 0

    for i in heads:
        if ( i != 'Orthogroup') and (i != 'Total'):
            data[i] = n
        n += 1

    return data


def get_gene(string):

    string = string.replace(' ', '')

    return string.split(',')


def read_orthogroup(file):

    data = {}

    for line in read_tsv(file, '\t'):
        if 'Orthogroup' in line[0]:
            head = read_head(line)
            continue

        group = []
        for i in head:
            genes = line[head[i]]
            try:
                n = int(genes)
            except:
                if not genes or genes=='':
                    n = 0
                else:
                    n = len(get_gene(genes))
            group.append(n)
        data[line[0]] = group

    return data, head


def counts(list, x):

    n = 0
    for i in list:
        if i == x:
            n += 1
    return n


def stat_gene_family(data, head, gene_dict):

    fam_dict = OrderedDict()
    #print('#Sample\tSingle\tMulti\tUnique\tOther\tUnclustered')
    for i in head:
        fam_dict[i] = [0, 0, 0, 0, 0]

    for i in data:
        line = data[i]
        gmax = max(line)
        gmin = min(line)
        if gmax==1 and gmin==1:
            for i in fam_dict:
                fam_dict[i][0] += 1
            continue
        if gmin>=1 and gmax>=2:
            for i in fam_dict:
                fam_dict[i][1] += line[head[i]-1]
            continue

        nc = counts(line, 0)
        if nc==len(line)-1:
            for i in fam_dict:
                fam_dict[i][2] += line[head[i]-1]
            continue
        if nc>=1:
            for i in fam_dict:
                fam_dict[i][3] += line[head[i]-1]
            continue

    for i in fam_dict:
        fam_dict[i][4] = gene_dict[i] - sum(fam_dict[i])

    return fam_dict


def stat_family_gene(data, head, gene_dict):

    fam_dict = OrderedDict()
    #print('#Sample\tGenes Number\tGenes Number In Families\tUnclustered Genes\tFamily Number\tUnique Families Number\tAverage Genes Number Per Family')
    for i in head:
        fam_dict[i] = [gene_dict[i], 0, 0, 0, 0, 0]

    for i in data:
        line = data[i]
        nc = counts(line, 0)
        for i in fam_dict:
            fam_dict[i][1] += line[head[i]-1]
            if line[head[i]-1] >0:
                fam_dict[i][3] += 1
            if nc==len(line)-1 and line[head[i]-1]>=0:
                fam_dict[i][4] += line[head[i]-1]

    for i in fam_dict:
        fam_dict[i][2] = fam_dict[i][0] - fam_dict[i][1]
        fam_dict[i][5] = fam_dict[i][1]*1.0/fam_dict[i][3]

    return fam_dict


def stat_orthogroups(proteins, groups):

    gene_dict = stat_gene(proteins)
    data, head = read_orthogroup(groups)
    fam_dict = stat_gene_family(data, head, gene_dict)
    fgene = stat_family_gene(data, head, gene_dict)

    fa = open('group.type.tsv', 'w')
    fa.write('#Sample\tSingle\tMulti\tUnique\tOther\tUnclustered\n')
    for i in fam_dict:
        line = fam_dict[i]
        fa.write('{0}\t{1:,}\t{2:,}\t{3:,}\t{4:,}\t{5:,}\n'.format(i, line[0], line[1], line[2], line[3], line[4]))
    fa.close()

    fg = open('group.stat.tsv', 'w')
    fg.write('#Sample\tGenes Number\tGenes Number In Families\tUnclustered Genes\tFamily Number\tUnique Families Number\tAverage Genes Number Per Family\n')
    for i in fgene:
        line = fgene[i]
        fg.write('{0}\t{1:,}\t{2:,}\t{3:,}\t{4:,}\t{5:,}\t{6:.2f}\n'.format(i, line[0], line[1], line[2], line[3], line[4], line[5]))
    fg.close()


def add_hlep_args(parser):

    parser.add_argument('-p', '--proteins', nargs='+', metavar='FILE', type=str, required=True,
        help='Input protein sequence file, format(fasta,fastq,fa.gz')
    parser.add_argument('-g', '--groups', metavar='STR', type=str, required=True,
        help='Input protein clustering results')

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
name:
    stat_orthogroups: Perform genome family analysis based on protein clustering results.
attention:
    stat_orthogroups -p protein*.faa -g Orthogroups.tsv
    stat_orthogroups -p protein*.faa -g Orthogroups.GeneCount.tsv
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    stat_orthogroups(args.proteins, args.groups)


if __name__ == "__main__":

    main()
