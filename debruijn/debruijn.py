#!/bin/env python3
# -*- coding: utf-8 -*-

#    This program is free software: you can redistribute it and/or modify
#    itunder the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html
#
#    Usage
#    -----
#      python -m debruijn -i data/eva71_two_reads.fq -k 21 -o contigs.fasta

"""Perform assembly based on debruijn graph."""


from random import randint
import argparse
# from operator import itemgetter
import os
import sys
import random
import statistics
import matplotlib.pyplot as plt
import networkx as nx
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
            try:
                yield next(filin).strip()
                next(filin)
                next(filin)
            except StopIteration:
                return


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
    for i in range(len(sequence)-kmer_size+1):
        try:
            yield sequence[i:kmer_size+i]
        except StopIteration:
            return


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
    graph: nexgraph
    """
    graph = nx.DiGraph()
    for kmer, poids in kmer_dico.items():
        graph.add_edge(kmer[:-1], kmer[1:], weight=poids)

    return graph


def save_graph(graph):
    """
    Save graph in a file.

    Affichage du graphe: a ne pas faire sur le jeu de données entier

    Parameter
    ---------
    graph: nexgraph
    """
    plt.subplot(111)
    nx.draw(graph, with_labels=True, font_weight='bold')
    plt.savefig("graph")


def get_starting_nodes(graph):
    """
    Parameter
    ---------
    graph: nexgraph

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


def save_contigs(contigs, output_file):
    """
    Save contigs to an output file.

    Parameter
    --------
    contigs: list
        list of tupple (contig, contig size)
    output_file: str
        an output file
    """
    with open(output_file, "w") as filout:
        for index, element in enumerate(contigs):
            filout.write(">contig_{} len={}\n".format(index, element[1]))
            filout.write("{}\n".format(fill(element[0])))


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def std(values):
    """
    Compute the standard deviation.

    Parameter
    ---------
    values: list
        list of value

    Return
    ------
    std: float
    """
    return statistics.stdev(values)


def path_average_weight(graph, path):
    """
    Compute the average weight of a given path.

    Parameter
    ---------
    graph: nexgraph
    path: list

    Return
    ------
    weight: float
    """
    weights = 0
    for i in range(len(path)-1):
        weights += graph[path[i]][path[i+1]]["weight"]
    return weights / (len(path)-1)


def remove_paths(graph, paths, delete_entry_node, delete_sink_node):
    """
    Clear the graph.

    Parameter
    ---------
    graph: nexgraph
    paths: list
        a list of paths
    delete_entry_node: boolean
        True entry node delete
    delete_sink_node: boolean
        True sink node delete

    Parameter
    ---------
    cleared_graph: nexgrap
    """
    for path in paths:
        if delete_entry_node:
            graph.remove_node(path[0])
        if delete_sink_node:
            graph.remove_node(path[-1])
        for node in path[1:-1]:
            if node:
                graph.remove_node(node)
    return graph


def select_best_path(graph, paths, path_length, path_weight,
                     delete_entry_node=False, delete_sink_node=False):
    """
    Clear the graph.

    Parameter
    ---------
    graph: nexgraph
    paths: list
        a list of paths
    path_length: list
        length of each path
    path_weight: list
        weight of each path
    delete_entry_node: boolean
        True entry node delete
    delete_sink_node: boolean
        True sink node delete

    Parameter
    ---------
    cleared_graph: nexgrap
    """
    # Get index of the "plus léger" path
    min_weight = [
        i for i, weight in enumerate(path_weight) if weight == min(path_weight)
    ]
    if len(min_weight) == 1:
        return remove_paths(graph, [paths[min_weight[0]]],
                            delete_entry_node, delete_sink_node)

    # Get index of the "plus petit" path
    min_length = [
        i for i, length in enumerate(path_length) if length == min(path_length)
    ]
    if len(min_length) == 1:
        return remove_paths(graph, [paths[min_length[0]]],
                            delete_entry_node, delete_sink_node)

    # Select random
    rand = randint(0, 1)
    return remove_paths(graph, [paths[rand]],
                        delete_entry_node, delete_sink_node)


def solve_bubble(graph, predecessor, successor):
    """
    Clear the graph of a bubble.

    Parameter
    ---------
    graph: nexgraph
    predecessor:
        predecessor node
    successor:
        successor node

    Return
    ------
    cleared_graph: nexgrap
    """
    tmp = list(nx.all_simple_paths(graph, predecessor, successor))

    while len(tmp) > 1:
        paths = [tmp[0], tmp[1]]
        path_weight = [
            path_average_weight(graph, tmp[0]),
            path_average_weight(graph, tmp[1])
        ]
        path_length = [len(tmp[0])-1, len(tmp[1])-1]

        graph = select_best_path(graph, paths, path_length, path_weight)

        tmp = list(nx.all_simple_paths(graph, predecessor, successor))

    return graph


def find_bubble(path1, path2):
    """
    Permet de déterminer le noeud d'entrée & de sortie d'une bulle, si présente.
    """
    entry, out = 0, 0
    flag = 0
    for i, _ in enumerate(path1):
        if flag == 0:
            if path1[i] not in path2:
                entry = path1[i-1]
                flag += 1
        else:
            if path1[i] in path2:
                out = path1[i]
                break
    return entry, out


def simplify_bubbles(graph):
    """
    Simplify the De Bruijn's graph.

    Parameter
    ---------
    graph: nexgraph

    Return
    ------
    cleared_graph: nexgraph
    """
    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)

    for predecessor in starting_nodes:
        for successor in sink_nodes:
            paths = list(nx.all_simple_paths(graph, predecessor, successor))
            while len(paths) > 1:
                entry, out = find_bubble(paths[0], paths[1])
                graph = solve_bubble(graph, entry, out)
                paths = list(nx.all_simple_paths(graph, predecessor, successor))

    return graph


def solve_entry_tips(graph, starting_nodes):
    """
    graph: nexgraph
    starting_nodes: list
    """
    pointes = []

    for starting in starting_nodes:
        successor = list(graph.successors(starting))[0]
        predecessors = list(graph.predecessors(successor))

        while len(predecessors) < 2:
            successor = list(graph.successors(successor))[0]
            predecessors = list(graph.predecessors(successor))

            if successor in get_sink_nodes(graph):
                break

        if successor not in get_sink_nodes(graph):
            pointes.append((starting, successor))

    while len(pointes) > 1:
        tmp = [
            list(nx.all_simple_paths(graph, pointes[0][0], pointes[0][1])),
            list(nx.all_simple_paths(graph, pointes[1][0], pointes[1][1]))
        ]
        path1, path2 = tmp[0][0], tmp[1][0]
        paths = [path1, path2]
        path_weight = [
            path_average_weight(graph, path1),
            path_average_weight(graph, path2)
        ]
        path_length = [len(path1)-1, len(path2)-1]
        graph = select_best_path(graph, paths, path_length, path_weight,
                                 delete_entry_node=True)
        if pointes[0][0] in graph.nodes():
            pointes.pop(1)
        else:
            pointes.pop(0)

    return graph


def solve_out_tips(graph, sink_nodes):
    """
    graph: nexgraph
    sink_nodes:
    """
    pointes = []

    for sink in sink_nodes:
        predecessor = list(graph.predecessors(sink))[0]
        successors = list(graph.successors(predecessor))

        while len(successors) < 2:
            predecessor = list(graph.predecessors(predecessor))[0]
            successors = list(graph.successors(predecessor))

            if predecessor in get_starting_nodes(graph):
                break

        if not predecessor in get_starting_nodes(graph):
            pointes.append((sink, predecessor))

    while len(pointes) > 1:
        tmp = [
            list(nx.all_simple_paths(graph, pointes[0][1], pointes[0][0])),
            list(nx.all_simple_paths(graph, pointes[1][1], pointes[1][0]))
        ]
        path1, path2 = tmp[0][0], tmp[1][0]
        paths = [path1, path2]
        path_weight = [
            path_average_weight(graph, path1),
            path_average_weight(graph, path2)
        ]
        path_length = [len(path1)-1, len(path2)-1]
        graph = select_best_path(graph, paths, path_length, path_weight,
                                 delete_sink_node=True)
        if pointes[0][0] in graph.nodes():
            pointes.pop(1)
        else:
            pointes.pop(0)

    return graph


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

    # Simplification of De Bruijn's graph
    graph = simplify_bubbles(graph)
    graph = solve_entry_tips(graph, get_starting_nodes(graph))
    graph = solve_out_tips(graph, get_sink_nodes(graph))

    # Contigs
    contigs = get_contigs(graph,
                          get_starting_nodes(graph), get_sink_nodes(graph))

    # Save contigs
    save_contigs(contigs, args.output_file)


if __name__ == '__main__':
    main()
