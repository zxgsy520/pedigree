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


def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.lstrip().rstrip('\n')

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_gff_anno(file):

    data = {}

    for line in read_tsv(file, "\t"):
        if line[2] != "CDS":
            continue
        contents = line[-1].split(";")
        geneid = contents[0].split("=")[1]
        annotate = ";".join(contents[3::])
        data[geneid] = annotate

    return data


def get_gene(string):

    string = string.replace(' ', '')

    return string.split(',')


def read_head(heads):

    data = OrderedDict()
    n = 0

    for i in heads:
        if ( i != 'Orthogroup') and (i != 'Total'):
            data[i] = n
        n += 1

    return data


def unique_gene_annotation(file, gff, sample):

    data = read_gff_anno(gff)

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

        if (group[head[sample]-1]) == 0 or (sum(group)-group[head[sample]-1]>0):
            continue

        string = line[0]
        for i in get_gene(line[head[sample]]):
            annotate = ""
            if i in data:
                annotate = data[i]
            string += "\t%s:%s" % (i, annotate)

        print(string)


def add_hlep_args(parser):

    parser.add_argument('-g', '--groups', metavar='STR', type=str, required=True,
        help='Input protein clustering results.')
    parser.add_argument('--gff', metavar='FILE', type=str, required=True,
        help='Input sample annotation file(gff).')
    parser.add_argument('-s', '--sample', metavar='STR', type=str, required=True,
        help='Input sample name.')

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
    unique_gene_annotation: Obtain unique genes from the results of gene family clustering.
attention:
    unique_gene_annotation -g Orthogroups.tsv --gff genome.gff -s name >unique_gene_annotation.tsv

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    unique_gene_annotation(args.groups, args.gff, args.sample)


if __name__ == "__main__":

    main()
