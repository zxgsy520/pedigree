#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
import numpy as np
from itertools import combinations
from scipy.optimize import curve_fit

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


def read_groups(file):

    head = []
    r = {}

    for line in read_tsv(file, "\t"):
        if 'Orthogroup' in line[0]:
            head = line[1::]
            continue
        temp = {}
        for i in range(1, len(head)+1, 1):
            if not line[i]:
                continue
            temp[head[i-1]] = line[i].replace(" ", "").split(",")
        r[line[0]] = temp

    return head, r


def all_true(list1, list2):

    jc = True
    for i in list1:
        if i not in list2:
            jc = False
            break
    return jc


def only_true(list1, list2):

    jc = False
    for i in list1:
        if i in list2:
            jc = True
            break
    return jc


def stat_gene(data, samples):

    core = 0
    pan = 0

    for line in data.values():
        if all_true(samples, line):
            core += 1
            pan += 1
            continue
        if only_true(samples, line):
            pan += 1

    return core, pan


def simulation_gene(data, head, number):

    cores = []
    pans = []

    for samples in combinations(head, number):
        samples = list(samples)
        core, pan = stat_gene(data, samples)

        cores.append(core)
        pans.append(pan)

    return sum(cores)*1.0/len(cores), sum(pans)*1.0/len(pans)


def stat_r2(y, py):

    my = sum(y)/len(y)
    sst= sum((y-my)**2)
    ssw = sum((y-py)**2)

    return 1-(ssw/sst)


def pan_func(x, a, b, c):
    
    return a*(x**b)+c


def pan_show(x, p_fit):

    a, b, c = p_fit.tolist()

    return a*(x**b)+c


def simu_pan(x, pan, gmax=31):

    pan_bounds = ([1000, -0.5, 0], [10000, 1, 10000])
    p_fit, pcov = curve_fit(pan_func, x, pan, bounds=pan_bounds)
    a, b, c = p_fit.tolist()
    pan_p = pan_show(x, p_fit)
    r2 = stat_r2(pan, pan_p)

    if max(x) <= gmax:
        x = np.array(range(1, gmax, 1))
        pan_p = pan_show(x, p_fit)

    print("#Pan-genome model: P={0:.6f}*n**{1:.6f}+{2:.6f}".format(a, b, c))
    print("#R2: {0:.6f}%".format(r2*100.0))

    return pan_p, x


def core_func(x, a, b, c):

    return a*np.exp(b*x)+c


def core_show(x, p_fit):

    a, b, c = p_fit.tolist()

    return a*np.exp(b*x)+c


def simu_core(x, core, gmax=31):

    core_bounds = ([500, -0.5, 0], [10000, 0, 10000])
    p_fit, pcov = curve_fit(core_func, x, core, bounds=core_bounds)
    a, b, c = p_fit.tolist() 
    core_p = core_show(x, p_fit)
    r2 = stat_r2(core, core_p)

    if max(x) <= gmax:
        x = np.array(range(1, gmax, 1))
        core_p = core_show(x, p_fit)
    print("#Core-genome model: C={0:.6f}*e**{1:.6f}n+{2:.6f}".format(a, b, c))
    print("#R2: {0:.6f}%".format(r2*100.0))

    return core_p, x


def plot_dual_axis(x, core, pan, x2, core_p, pan_p):

    fig = plt.figure(figsize=(8, 5.5))
    ax = fig.add_subplot(111)

    plt.grid(True, which='minor', axis='both', lw=1.5, color='#E5C700', alpha=0.3)
    plt.grid(True, which='major', axis='both', lw=2, color='#E2BDD5', alpha=0.3)
    ax.tick_params(axis='both', which='both', color='#ffffff', length=0, width=2)
    font1 = {'family': 'Times New Roman', 'weight': 'normal', 'color': '#212121', 'size': 16}

    ax.set_xlabel('Genes', font1)
    ax.set_ylabel('Genome number', font1)
    ax.plot(x, core, "o", color='#212121', markersize=8, label='Core gene')
    ax.plot(x2, core_p, "-.", color='#212121')
    font2 = {'family': 'Times New Roman', 'weight': 'normal', 'color': '#212121', 'size': 16}
    
    ax.plot(x, pan, "o", color='#FFA000', markersize=8, label='Pan gene')
    ax.plot(x2, pan_p, "-.", color='#FFA000')

    #plt.legend(loc='center right', bbox_to_anchor=(2.0, 0.64) ,frameon=False)
    plt.legend(loc='center right', frameon=False)
    plt.savefig('core_pan_gene.png', dpi=700)
    plt.savefig('core_pan_gene.pdf')


def simulation_pan_genome(file, gmax=6):

    head, data = read_groups(file)
    x = []
    cores = []
    pans = []

    for i in range(1, len(head)+1, 1):
        core, pan = simulation_gene(data, head, i)
        x.append(i)
        cores.append(core)
        pans.append(pan)

    x = np.array(x)
    cores = np.array(cores)
    pans = np.array(pans)
    pans_p, nx = simu_pan(x, pans, gmax)
    cores_p, nx = simu_core(x, cores, gmax)

    print("Genome number\tCore genes\tPredict core genes\tPan genes\tPredict pan genes")
    for i in x:
        print("{0}\t{1:.2f}\t{2:.2f}\t{3:.2f}\t{4:.2f}".format(i, cores[i-1], cores_p[i-1], pans[i-1], pans_p[i-1]))

    plot_dual_axis(x, cores, pans, nx, cores_p, pans_p)
        
    return 0


def add_hlep_args(parser):

    parser.add_argument('input', metavar='FILE', type=str,
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
    simulation_pan_genome.py: Simulation statistics of core gene number and pan gene number
attention:
    simulation_pan_genome.py Orthogroups.tsv >stat_pan_gene.tsv
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    simulation_pan_genome(args.input)


if __name__ == "__main__":

    main()
