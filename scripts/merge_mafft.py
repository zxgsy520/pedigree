#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_fasta(file):

    '''Read fasta file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)

    seq = []
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            line = line.strip(">")
            if len(seq) == 2:
                yield seq
            seq = []
            seq.append(line.split()[0])
            continue
        if len(seq) == 2:
            seq[1] += line
        else:
            seq.append(line)

    if seq == 2:
        yield seq
    fp.close()


def merge_mafft(files):

    data = {}
    for file in files:
        for seqid, seq in read_fasta(file):
            seqid = seqid.split("|")[0]
            seq = seq.replace(" ", "").upper()
            if seqid not in data:
                data[seqid] = ""
            data[seqid] += seq

    for i in data:
        print(">%s\n%s" % (i, data[i]))

    return 0


def add_hlep_args(parser):

    parser.add_argument('-m', '--mafft', nargs='+', metavar='FILE', type=str, required=True,
        help='Input mafft aligned file')

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
    merge_mafft.py: Merge files aligned with mafft
attention:
    merge_mafft.py -m *.mafft >all_aligned.fasta
    merge_mafft.py -m *.mafft-gb >all_aligned.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    merge_mafft(args.mafft)


if __name__ == "__main__":

    main()
