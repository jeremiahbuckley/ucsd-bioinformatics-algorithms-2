#! /bin/python3

import sys
import time


def split_dna_into_kmers(dna, kmer_len):
    kmers = []
    for i in range(len(dna)-kmer_len+1):
        kmer = dna[i:i+kmer_len]
        kmers.append(kmer)
    return sorted(kmers)


if __name__ == "__main__":
    starttime = time.process_time()

    if len(sys.argv) < 2:
       Print("Usage: param 1=filepath. File format: <int_kmer_len>\n<string_dna>")

    with open(sys.argv[1]) as fl:
        kmer_len = int(fl.readline())
        dna = fl.readline()
        dna = dna[0:len(dna)-1]

    kmers = split_dna_into_kmers(dna, kmer_len)
    print(*kmers, sep="\n")

    endtime = time.process_time()
    print("Time: {0}".format(endtime-starttime))


