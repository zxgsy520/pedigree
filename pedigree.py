#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import gzip
import hashlib
import logging
import argparse

from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.lstrip().rstrip("\n")

        if not line:
            continue

        yield line.split(sep)

    fh.close()


def read_merge_annotate(file):

    r = {}

    for line in read_tsv(file, "\t"):
        if line[0].startswith("#"):
            continue
        r[line[0]] = line

    return r


def read_head(heads):

    data = OrderedDict()
    n = 0

    for i in heads:
        if "." in i:
            i = i.split(".")[0]
        if n >= 1:
            data[i] = n
        n += 1

    return data


def get_prefix(file):

    file = file.split("/")[-1]

    if "." in file:
        prefix = file.split(".")[0]
    else:
        prefix = file

    return prefix


def get_gene_count(group):

    genes = []
    counts = []
    n = 0

    for i in group:
        if n >= 1:
            if i:
                i = i.replace(" ", "").split(",")
                genes.append(i)
                counts.append(len(i))
            else:
                genes.append([])
                counts.append(0)
        n += 1

    return genes, counts
########-------------Unique gene annotation-------------------------------#####
def ref_unassigned(file):

    data = {}
    n = 0
    for line in read_tsv(file, "\t"):
        n += 1
        if n == 1:
            group = read_head(line)
            continue
        for i in group:
            site = group[i]
            if i not in data:
                data[i] = []
            if not line[site]:
                continue
            data[i].append(line[site])

    return data


def get_unigene(args):

    data = ref_unassigned(args.input)

    funcs = {}
    for file in args.function:
        sample = get_prefix(file)
        funcs[sample] = read_merge_annotate(file)

    fo = open("stat_unique_gene.xls", "w")
    fo.write("#Sample\tUnique gene count\n")

    print("#Sample\tGene_Id\tStrand\tStart\tEnd\tGene_Length(bp)\tLocation\tGene_Name\tRefseq_Description\tPfam_Id\tPfam_Description\tTIGRFAMs_Id\tTIGRFAMs_Description\tCOG_Id\tCOG_Description\tCOG_Type\tKO_Id\tKO_Description\tPathway\tGO_Id\tGO_Description")
    for i in data:
        genes = data[i]
        fo.write("%s\t%s\n" % (i, len(genes)))
        if i not in funcs:
            LOG.info("Sample %s has no annotation information" % i)
            temp = []
        else:
            temp = funcs[i]
        for j in genes:
            if j not in temp:
                print("%s\t%s" % (i, j))
            else:
                print("%s\t%s" % (i, "\t".join(temp[j])))

    fo.close()

    return 0


def add_get_unigene_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input protein clustering results(Orthogroups_UnassignedGenes.tsv).")
    parser.add_argument("-f", "--function", nargs="+", metavar="FILE", type=str,
        help="Input gene function integration(*.merge.annotate.xls).")

    return parser

#######-----------------------------------------------------------------######

########-------------Core gene annotation-------------------------------#####
def ref_core_orthogroups(file):

    data = {}
    n = 0
    fo = open("shared_gene.tsv", "w")

    for line in read_tsv(file, "\t"):
        n += 1
        if n == 1:
            fo.write("%s\n" % "\t".join(line))
            continue
        genes, counts = get_gene_count(line)
        #print(counts)
        if min(counts) == 0:
            continue
        fo.write("%s\n" % "\t".join(line))
        gene = genes[0][0]
        data[line[0]] = gene
    fo.close()

    return data


def get_core_gene(args):

    data = ref_core_orthogroups(args.input)
    funcs = {}

    for file in args.function:
        funcs.update(read_merge_annotate(file))

    print("#Group_Id\tRef_Gene_Id\tStrand\tStart\tEnd\tGene_Length(bp)\tLocation\tGene_Name\tRefseq_Description\tPfam_Id\tPfam_Description\tTIGRFAMs_Id\tTIGRFAMs_Description\tCOG_Id\tCOG_Description\tCOG_Type\tKO_Id\tKO_Description\tPathway\tGO_Id\tGO_Description")
    for i in data:
        gene = data[i]
        if gene in funcs:
            print("%s\t%s" % (i, "\t".join(funcs[gene])))
        else:
            print("%s\t%s" % (i, gene))
    return 0


def add_core_gene_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input protein clustering results(Orthogroups.tsv).")
    parser.add_argument("-f", "--function", nargs="+", metavar="FILE", type=str,
        help="Input gene function integration(*.merge.annotate.xls).")

    return parser
#######-----------------------------------------------------------------######
def add_pedigree_parser(parser):

    subparsers = parser.add_subparsers(
        title='command',
        dest='commands')
    subparsers.required = True

    get_unigene_parser = subparsers.add_parser("get_unigene", help="Obtain unique genes and annotation results")
    get_unigene_parser = add_get_unigene_args(get_unigene_parser)
    get_unigene_parser.set_defaults(func=get_unigene)

    core_gene_parser = subparsers.add_parser("get_core_gene", help="Obtain core genes or common genes")
    core_gene_parser = add_core_gene_args(core_gene_parser)
    core_gene_parser.set_defaults(func=get_core_gene)

    return parser


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
name:
pedigree：Tools for microbial pan-genome analysis.
URL：https://github.com/zxgsy520/pedigree

version: %s
contact:  %s <%s>\
        """ % (__version__, " ".join(__author__), __email__))

    parser = add_pedigree_parser(parser)
    args = parser.parse_args()

    args.func(args)

    return parser.parse_args()


if __name__ == "__main__":

    main()
