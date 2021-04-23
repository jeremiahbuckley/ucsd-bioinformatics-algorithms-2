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
            if len(kvp[1].followingNodes) > 0:
                print("{0} -> {1}".format(kmer, ",".join(kvp[1].followingNodes)))

def build_fragment_graph(kmers):
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
    return graph


if __name__ == "__main__":
    start = time.process_time()

    if len(sys.argv) < 2:
        print("Usage: one param, filename, format is <list_str_of_dna_kmers>")

    with open(sys.argv[1]) as f:
        kmers = [line.rstrip() for line in f]

    graph = build_fragment_graph(kmers)
    print_graph(graph)

    end = time.process_time()
    print("Time: {0}".format(end-start))
