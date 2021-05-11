#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_fasta(file):

    '''Read fasta file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".faa") or file.endswith(".fna"):
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


def read_fastq(file):

    '''Read fastq file'''
    if file.endswith("fastq.gz") or file.endswith(".fq.gz"):
        fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
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
        if not seq:
            seq.append(line.split()[0].strip("@"))
            continue
        seq.append(line)
        if len(seq) == 4:
            yield seq
            seq = []

    fp.close()


def filter_protein(file, minlen, suffix="fa"):

    if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
        fh = read_fastq(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fasta.gz") or file.endswith(".fa.gz") or file.endswith(".faa"):
        fh = read_fasta(file)
    else:
        raise Exception("%r file format error" % file)

    name = file.split("/")[-1].split(".")[0:-1]
    name = ".".join(name)
    fo = open("%s.%s" % (name, suffix), "w")

    for line in fh:
        seqid, seq = line[0], line[1]
        seqid = seqid.split("|")[-1]
        seqlen = len(seq)
        if seqlen <= minlen:
            LOG.info("%s\t%s" % (seqid, seqlen))
            continue
        if "." in seq or "*" in seq:
            continue
        fo.write(">%s\n%s\n" % (seqid, seq))

    fo.close()


def filter_proteins(files, minlen, suffix="fa"):

    for file in files:
        filter_protein(file, minlen, suffix)
    return 0


def add_hlep(parser):

    parser.add_argument('input', nargs='+', metavar='STR',type=str,
        help='Input file.')
    parser.add_argument("--minlen", metavar="INT", type=int, default=30,
        help="Filter the length value of the protein, default=30")
    parser.add_argument("-s", "--suffix", metavar="STR", type=str, default="fa",
        help="Set the output suffix, default=fa")

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
    filter_proteins.py Filter abnormal proteins

attention:
    filter_proteins.py reads.fastq -s fa
    filter_proteins.py reads.fasta -s fa
''')
    args = add_hlep(parser).parse_args()

    filter_proteins(args.input, args.minlen, args.suffix)


if __name__ == "__main__":

    main()
