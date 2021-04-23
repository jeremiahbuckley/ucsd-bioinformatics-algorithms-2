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
    pass

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

def create_de_bruijn_graph(dna, kmer_len):
    graph = {}
    prev_kmer = ""
    for i in range(len(dna) - kmer_len + 2):
        kmer = dna[i:i+kmer_len-1]
        if kmer not in graph:
            graph[kmer] = []
        if len(prev_kmer) > 0: # and kmer not in graph[prev_kmer]:
            graph[prev_kmer].append(kmer)
        prev_kmer = kmer

    bad_keys = []
    for i in graph.keys():
        if len(graph[i]) == 0:
            print("got one")
            bad_keys.append(i)

    for i in range(len(bad_keys)):
        del graph[bad_keys[i]]
    return graph


if __name__ == "__main__":
    start = time.process_time()

    if len(sys.argv) < 2:
        print("Usage: one param, filename, format is <int_kmer_len>\n<str_dna>")

    with open(sys.argv[1]) as f:
        kmer_len = int(f.readline())
        dna = f.readline().rstrip()

    graph = create_de_bruijn_graph(dna, kmer_len)
    outlist = []
    for kvp in graph.items():
        outlist.append("{0} -> {1}".format(kvp[0], ",".join(kvp[1])))
    print(*sorted(outlist),sep="\n")

    #end = time.process_time()
    #print("Time: {0}".format(end-start))
