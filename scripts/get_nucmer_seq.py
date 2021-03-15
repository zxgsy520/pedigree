#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.2.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def rm_dup(string, sep=" "):

    string = string.strip()
    r = ""

    for i in range(len(string)):
        if string[i]==sep and string[i+1] == sep and i<len(string)-1:
            continue
        r += string[i]

    return r


def read_coords(file):

    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("="):
            continue
        if "[S1]" in line:
            continue

        line = line.replace("|", "")
        line = rm_dup(line).split()

        if len(line) < 10:
            continue
   

        yield line


def get_max_len(file):

    max_len = 0

    for line in read_coords(file):
        if int(line[5]) <= max_len:
            continue
        max_len = int(line[5])
        seq_pos = [line[10], int(line[2]), int(line[3])]

    return seq_pos


def read_fasta(file):

    '''Read fasta file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

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

    if len(seq) == 2:
        yield seq
    fp.close()


def complement(seq):

    cdict = {"A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    seq = list(seq.upper())
    nseq = ""
    for i in seq:
        nseq += cdict[i]

    return nseq


def reverse_complement(seq):

    seq = seq[::-1]

    return complement(seq)


def get_nucmer_seq(coords, file):

    seq_pos = get_max_len(coords)
    sample = file.split('/')[-1].split('.')[0]

    for seqid, seq in read_fasta(file):
        if seqid not in seq_pos:
            continue
        direct = "+"
        start, end = seq_pos[1], seq_pos[2]
        if start >= end:
            direct = "-"
            start, end = end, start
        seq = seq[start-1:end]

        if direct == "-":
            seq = reverse_complement(seq)

        print(">%s|%s-%s\n%s" % (sample, end, start, seq))

    return 0


def add_hlep_args(parser):

    parser.add_argument('-c', '--coords', metavar='FILE', type=str, required=True,
        help='Input the nucmer comparison result file')
    parser.add_argument('-g', '--genome', metavar='FILE', type=str, required=True,
        help='Input genome file.')

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
    get_nucmer_seq.py: Get the target sequence from the nucmer alignment result
attention:
    get_nucmer_seq.py -c group.coords -s genome.fasta >seq.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    get_nucmer_seq(args.coords, args.genome)


if __name__ == "__main__":

    main()
