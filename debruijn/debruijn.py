#!/bin/env python3
# -*- coding: utf-8 -*-

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""


from random import randint
import argparse
from operator import itemgetter
import os
import sys
import networkx as nx
import matplotlib.pyplot as plt
import random
import statistics
random.seed(9001)


__author__ = "IMBERT Pierre"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["IMBERT Pierre"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "IMBERT Pierre"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq):
    """
    Read a fastq file and return the sequences in this file.

    Parameter
    ---------
    fastq: str
        name of the file

    Return
    ------
    a generator of sequences with yield
    """
    with open(fastq, "r") as filin:
        for _ in filin:
            yield next(filin)
            next(filin)
            next(filin)


def cut_kmer(sequence, kmer_size):
    """
    Generate kmer of size k of the sequence.

    Parameter
    ---------
    sequence: str
        a specific sequence
    kmer_size: int
        size of the kmer

    Return
    ------
    a generator of kmer for a specific sequence
    """
    for i in range(len(sequence)-kmer_size):
        yield sequence[i:kmer_size+i]


def build_kmer_dict(fastq, kmer_size):
    """

    Paramater
    ---------
    fastq: str
        name of the file
    kmer_size: int
        size of the kmer

    Return
    ------
    kmer_dico: dictionary
        keys: kmer, values: nombre d'occurence of the kmer in the fastq file
    """
    kmer_dico = {}
    for sequence in read_fastq(fastq):  # get sequences
        for kmer in cut_kmer(sequence, kmer_size):  # get kmer
            if kmer not in kmer_dico.keys():
                kmer_dico[kmer] = 0
            kmer_dico[kmer] += 1
    return kmer_dico


def build_graph(kmer_dico):
    """
    Creation of a graph of kmer.

    Parameter
    ---------
    kmer_dico: dictionary
        dictionary of kmer

    Return
    ------
    directed_graph: nx graph
        a graph
    """
    directed_graph = nx.DiGraph()
    for kmer, poids in kmer_dico.items():
        directed_graph.add_edge(kmer[:-1], kmer[1:], weight=poids)

    # Affichage du graphe: a ne pas faire sur le jeu de données entier
    # plt.subplot(111)
    # nx.draw(directed_graph, with_labels=True, font_weight='bold')
    # plt.savefig("graph")
    return directed_graph


def save_graph(directed_graph):
    """
    Save graph in a file.

    Parameter
    ---------
    directed_graph: nexgraph
        a graph
    """
    plt.subplot(111)
    nx.draw(directed_graph, with_labels=True, font_weight='bold')
    plt.savefig("graph")


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Kmer dictionary
    kmer_dico = build_kmer_dict(args.fastq_file, args.kmer_size)
    build_graph(kmer_dico)


if __name__ == '__main__':
    main()
