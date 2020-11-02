#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse
import matplotlib
matplotlib.use('Agg')

from collections import OrderedDict
from matplotlib import pyplot as plt

LOG = logging.getLogger(__name__)

__version__ = "1.2.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


#PCOLOR = ["#FF6630", "#D9C125", "#18F04E", "#4A7EDB", "#FF40FF"]
#PCOLOR = ["#FFDE87", "#B3DB65", "#56F0A8", "#4285DB", "#B754FF"]
PCOLOR = ["#FFBC7D", "#BDDB5C", "#4CF0E6", "#7939DB", "#FF6642"]

def str2float(string):

    string = string.replace(',', '')

    return float(string)


def read_group_type(file):

    data = OrderedDict()
    samples = []

    for line in open(file, 'r'):
        line = line.strip()
        if not line:
            continue

        line = line.split('\t')
        if "Sample" in line[0]:
            types = line[1::]
            continue
        samples.append(line[0].replace('_', ' '))

        for i in range(len(types)):
            if types[i] not in data:
                data[types[i]] = []
            data[types[i]].append(str2float(line[i+1]))

    return data, samples


def inter_axis_area(data):

    inter_dice = OrderedDict()

    for i in data:
        lines = []
        for j in data[i]:
            lines.append(j)
            lines.append(j)
        inter_dice[i] = lines

    return inter_dice


def add_list(list1, list2):

    lists = []

    for i in range(len(list1)):
        lists.append(list1[i]+list2[i])
    return lists


def produce_x(data, samples):

    x = [i*2 for i in range(len(samples))]
    x2=[]
    for i in range(len(x)):
        a=x[i]-0.5
        b=x[i]+0.5
        x2.append(a)
        x2.append(b)

    lines = [i*0 for i in range(len(samples))]
    for i in data:
        lines = add_list(lines, data[i])
        data[i] = lines

    return x, x2, data


def plot_bar(data, samples, prefix, height):

    inter_color = [(253/255.0,222/255.0,194/255.0), (188/255.0,205/255.0,225/255.0), (239/255.0,194/255.0,188/255.0), (160/255.120,150/255.0,200/255.0), (215/255.0,215/255.0,215/255.0)]
    x, x2, data =  produce_x(data, samples)
    inter_dice = inter_axis_area(data)

    if len(samples) <=15:
        fig = plt.figure(figsize=[10,8.5], facecolor=(239/255.0,248/255.0,253/255.0))
    else:
        fig = plt.figure(figsize=[14,8.5], facecolor=(239/255.0,248/255.0,253/255.0))
    ax1=fig.add_subplot(1,1,1,facecolor=(239/255.0,248/255.0,253/255.0))
    ax1.spines['top'].set_visible(False)#去掉上边边框
    ax1.spines['bottom'].set_visible(False)#去掉下方边边框
    ax1.spines['right'].set_visible(False)#去掉右边边框
    ax1.spines['left'].set_visible(False)#去掉左边边框
    ax1.grid(True, 'major', 'y', ls='--', lw=.5, c='black', alpha=.3)
    ax1.xaxis.set_major_formatter(plt.FuncFormatter(''.format)) #X轴不显示刻度
    ax1.xaxis.set_minor_formatter(plt.FuncFormatter(''.format))

    plt.subplots_adjust(left=0.1, right=0.8, top=0.95, bottom=height) #设置页边距
    #plt.yticks(range(0, 3000, 500), fontsize=10)

    dlen = len(inter_dice)
    types = list(inter_dice.keys())

    if len(samples) <=5:
        n = 0
        for i in types[::-1]:
            ax1.stackplot(x2, inter_dice[i],color=inter_color[n], alpha=1)
            n += 1

    n = 0
    for i in types[::-1]:
        ax1.bar(x, data[i], width=1.6, label=i ,color=PCOLOR[n])
        n += 1

#    ax1.plot([min(x),max(x)],[0,0],color='black',linewidth=5)
    for a,b in zip(x, samples):
       ax1.text(a,-2,b,ha='right',va='top',fontsize=12,rotation=45)

    plt.legend(loc='center right', bbox_to_anchor=(1.22, 0.5),frameon=False)
#    plt.legend(loc='center right')

    plt.ylabel('Gene Number', fontsize=14)
    plt.savefig("%s.png" % prefix, dpi=500)
    plt.savefig("%s.pdf" % prefix)


def plot_args(parser):
    parser.add_argument('gtype',
        help='Input the classification statistics of the gene family.')
    parser.add_argument('-b', '--bottom', metavar='FLOAT', type=float, default=0.25,
        help='Set footer height, default=0.25.')
    parser.add_argument('-p', '--prefix', metavar='STR', type=str, default='out',
        help='Output result prefix, default=out.')

    return parser


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
name:
    bar_group_type.py.py -- Map gene family distribution.
attention:
    bar_group_type.py.py group.type.xls -o name.group.type

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = plot_args(parser)
    args = parser.parse_args()

    data, samples  = read_group_type(args.gtype)
    plot_bar(data, samples, args.prefix, args.bottom)


if __name__ == "__main__":
    main()
