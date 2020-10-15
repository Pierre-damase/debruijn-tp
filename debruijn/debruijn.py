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
            yield next(filin).strip()
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

    Example:
      k-mer ACTG see 10 times

      Node ACT ------> Node CTG
              weight=10

    Parameter
    ---------
    kmer_dico: dictionary
        dictionary of kmer

    Return
    ------
    graph: nx graph
        a graph
    """
    graph = nx.DiGraph()
    for kmer, poids in kmer_dico.items():
        graph.add_edge(kmer[:-1], kmer[1:], weight=poids)

    # Affichage du graphe: a ne pas faire sur le jeu de données entier
    plt.subplot(111)
    nx.draw(graph, with_labels=True, font_weight='bold')
    plt.savefig("graph")
    return graph


def save_graph(graph):
    """
    Save graph in a file.

    Parameter
    ---------
    graph: nexgraph
        a graph
    """
    plt.subplot(111)
    nx.draw(graph, with_labels=True, font_weight='bold')
    plt.savefig("graph")


def get_starting_nodes(graph):
    """
    Parameter
    ---------
    graph: nexgraph
        a graph

    Return
    ------
    starting_nodes: list
        list of starting nodes
    """
    nodes, starting_nodes = list(graph.nodes()), []
    for node in nodes:
        if not list(graph.predecessors(node)):
            starting_nodes.append(node)
    return starting_nodes


def get_sink_nodes(graph):
    """
    Parameter
    ---------
    graph: nexgraph
        a graph

    Return
    ------
    sink_nodes: list
        list of output nodes
    """
    nodes, sink_nodes = list(graph.nodes()), []
    for node in nodes:
        if not list(graph.successors(node)):
            sink_nodes.append(node)
    return sink_nodes


def get_contigs(graph, starting_nodes, sink_nodes):
    """
    Déterminer chaque chemin d'un graphe entre un starting & sink nodes.

    Parameter
    ---------
    graph: nexgraph
        a graph
    starting_nodes: list
        list of input nodes
    sink_nodes: list
        list of output nodes

    Return
    ------
    contigs: list
        list of tupple (contig, contig size)
    """
    contigs = []
    for starting in starting_nodes:
        for sink in sink_nodes:
            tmp = list(nx.all_simple_paths(graph, starting, sink))

            if tmp:
                contig = tmp[0][0]
                for i in range(1, len(tmp[0])):
                    contig += tmp[0][i][-1]
                contigs.append((contig, len(contig)))

    return contigs


def save_contigs(output_file, contigs):
    """

    Parameter
    --------
    output_file: str
        an output file
    contigs: list
        list of tupple (contig, contig size)
    """
    with open(output_file, "w") as filout:
        for index, element in enumerate(contigs):
            tmp = ">contig_" + str(index) + " len=" + str(element[1]) + "\n"
            filout.write(tmp)
            filout.write(fill(element[0]))
            filout.write("\n")


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def std():
    pass


def path_average_weight():
    pass


def remove_paths():
    pass


def select_best_path():
    pass


def solve_bubble():
    pass


def simplify_bubbles():
    pass


def solve_entry_tips():
    pass


def solve_out_tips():
    pass


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

    # De Bruijn graph
    graph = build_graph(kmer_dico)

    # Contigs
    contigs = get_contigs(graph,
                          get_starting_nodes(graph), get_sink_nodes(graph))

    # Save contigs
    save_contigs(args.output_file, contigs)

if __name__ == '__main__':
    main()
