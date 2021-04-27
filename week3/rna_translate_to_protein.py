#! /bin/python3

import sys
import time
import math

STOP_PROTEIN_STR = "STP"

def create_protein_lookup(rna_protein_table):
    proteins_by_rna = {}
    for line in rna_protein_table:
        space_idx = line.find(" ")
        if space_idx > -1:
            key = line[0:space_idx]
            value = line[space_idx+1:len(line)]
        else:
            key = line
            value = STOP_PROTEIN_STR
        proteins_by_rna[key] = value

    for kvp in proteins_by_rna.items():
        print("{0} {1}".format(kvp[0], kvp[1]))

    return proteins_by_rna    

def translate_rna_to_protein(rna_str, rna_protein_table):
    proteins_by_rna = create_protein_lookup(rna_protein_table)

    proteins_str = ""
    for i in range(math.ceil(len(rna_str)/3)):
        protein_codon = rna_str[i*3:i*3+3]
        if protein_codon in proteins_by_rna:
            if proteins_by_rna[protein_codon] != STOP_PROTEIN_STR:
                proteins_str += proteins_by_rna[protein_codon]
        else:
            raise ValueError("Could not find codon [{0}] in lookup.".format(protein_codon))
    return proteins_str


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: one param, filename, format:<str_rna>.") 

    else:
        start = time.process_time()

        with open(sys.argv[1]) as f:
            lines = [line.rstrip() for line in f]

        with open("./RNA_codon_table_1.txt") as f:
            rna_protein_table = [line.rstrip() for line in f]

        protein = translate_rna_to_protein(lines[0], rna_protein_table)
        print(protein)


        end = time.process_time()
        #print("Time: {0}".format(end-start))

