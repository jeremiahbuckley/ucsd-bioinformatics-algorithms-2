#! /bin/python3

import sys
import time

class GraphNodes:

    def __init__(self):
        self.prefixNodes = []
        self.followingNodes = []
'''
builds a dictionary like this:
[ACGT] -> [AACGT] -> ["ACGTA", "ACGTC", ACGTG", etc...]
          [CACGT] -> ["ACGTA", "ACGTC
          [GACGT]
          [TACGT]

'''

def print_graph(graph):
    for kvp in graph.items():
        #print(kvp[0])
        for kmer in kvp[1].prefixNodes:
            follow_str = "None"
            if len(kvp[1].followingNodes) > 0:
                follow_str =  ",".join(kvp[1].followingNodes)
            print("{0} -> {1}".format(kmer, follow_str))

def build_fragment_graph(kmers, kmer_len):
    graph = {}

    for i in range(len(kmers)):
        kmer = kmers[i]
        kmer_as_key = kmer[1:len(kmer)]
        kmer_as_lookup = kmer[0:len(kmer)-1]

        if kmer_as_key not in graph:
            graph[kmer_as_key] = GraphNodes()
        if kmer not in graph[kmer_as_key].prefixNodes:
            graph[kmer_as_key].prefixNodes.append(kmer)

        if kmer_as_lookup not in graph:
            graph[kmer_as_lookup] = GraphNodes()
        if kmer not in graph[kmer_as_lookup].followingNodes:
            graph[kmer_as_lookup].followingNodes.append(kmer)

        #print(".." + kmer)
        #print(".." + kmer_as_key)
        #print(".." + kmer_as_lookup)
        #print_graph(graph)
        #print("--")        
    return graph


if __name__ == "__main__":
    start = time.process_time()

    if len(sys.argv) < 2:
        print("Usage: one param, filename, format:<int_kmer_len>\n<list_str_of_dna_kmers>.\nPositional switch -n will stop reading 1st line as kmer_len.")

    first_line_is_kmer_len = True
    if len(sys.argv) >= 3:
        if sys.argv[2] == "-n":
            first_line_is_kmer_len = False

    with open(sys.argv[1]) as f:
        kmers = [line.rstrip() for line in f]

    if not first_line_is_kmer_len:
        kmer_len = len(kmers[0])
    else:
        kmer_len = -1
        try:
            int(kmers[0])
            kmer_len = int(kmers.pop(0))
        except ValueError:
            kmer_len = len(kmers[0])

    graph = build_fragment_graph(kmers, kmer_len)
    print_graph(graph)

    end = time.process_time()
    #print("Time: {0}".format(end-start))
