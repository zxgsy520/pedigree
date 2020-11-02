#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_fasta_ann(file):

    LOG.info("reading message from %r" % file)

    if file.endswith("gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)

    data = {}
    for line in fp:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()

        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            line = line.strip(">")
            seqid = line.split()[0]
            try:
                ann = line.split(":")[1].strip()
            except:
                ann = " ".join(line.split()[1::])
            ann = ann.split("[")[0].strip()
            data[seqid] = ann
    fp.close()
    return data


def read_ann(files):

    data = {}

    for file in files:
        sample = file.split("/")[-1].split(".")
        sample = ".".join(sample[0:-1])
        ann = read_fasta_ann(file)
        data[sample] = ann

    return data


def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.lstrip().rstrip('\n')

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


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


def annotate_gene(genes, data):

    ann = []

    for i in genes:
        if i not in data:
            LOG.info("Gene %s does not exist" % i)
            continue
        ann.append("%s:%s" % (i, data[i]))
    return ann


def annotate_groups(group, files):

    data = read_ann(files)

    for line in read_tsv(group, '\t'):
        if 'Orthogroup' in line[0]:
            head = read_head(line)
            print('\t'.join(line))
            continue

        lines = [line[0]]
        for i in head:
            genes = get_gene(line[head[i]])
            anns = annotate_gene(genes, data[i])
            lines.append(';'.join(anns))
        print('\t'.join(lines))

    return 1


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
    annotate_groups: Annotated gene family clustering results.
attention:
    annotate_groups -p protein*.faa -g Orthogroups.tsv >annotate_groups.tsv
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    annotate_groups(args.groups, args.proteins)


if __name__ == "__main__":

    main()
