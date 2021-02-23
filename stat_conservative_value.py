#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

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


def stat_gene(file):

    for line in read_tsv(file, "\t"):
        if 'Orthogroup' in line[0]:
            continue

        cv = 0
        gene = 0

        for i in line[1::]:
            i = i.strip().replace(" ", "")
            if not i:
                continue
            cv += 1
            gene += len(i.split(","))
            #gene += 1
        yield cv, gene


def stat_conservative_value(files):

    r = {}

    for file in files:
        for cv, gene in stat_gene(file):
            if cv not in r:
                r[cv] = 0
            r[cv] += 1

    print("#Conservative Value\tGenome Number")
    for i in sorted(r.keys()):
        print("%d\t%d" % (i, r[i]))

    return 0


def add_hlep_args(parser):

    parser.add_argument('-i', '--input', nargs='+', metavar='FILE', type=str, required=True,
        help='Input gene cluster table')

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
    stat_conservative_value: Statistical conservative value
attention:
    stat_conservative_value.py -i Orthogroups.tsv Orthogroups_UnassignedGenes.tsv >stat_conservative_value.tsv
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    stat_conservative_value(args.input)


if __name__ == "__main__":

    main()
