#! /bin/python3

import sys
import time
from lib.utils import utils

STOP_PROTEIN_STR = "STP"

def create_protein_rna_lookup(rna_protein_table):
    protein_by_rna = {}
    rna_by_protein = {}

    for line in rna_protein_table:
        space_idx = line.find(" ")
        if space_idx > -1:
            key = line[0:space_idx]
            value = line[space_idx+1:len(line)]
        else:
            key = line
            value = STOP_PROTEIN_STR
        protein_by_rna[key] = value
        if value not in rna_by_protein:
            rna_by_protein[value] = []
        rna_by_protein[value].append(key)

    #for kvp in protein_by_rna.items():
    #    print("{0} {1}".format(kvp[0], kvp[1]))
    #for kvp in rna_by_protein.items():
    #    print("{0} {1}".format(kvp[0], ",".join(kvp[1])))

    return protein_by_rna, rna_by_protein

def find_peptide_fragments_in_rna(rna, is_complement_strand, peptide_sequence, rna_protein_table):
    protein_by_rna, rna_by_protein = create_protein_rna_lookup(rna_protein_table)
    first_peptide_rna_list = rna_by_protein[peptide_sequence[0]][:]
    #print("<<" + ",".join(first_peptide_rna_list) + ">>")
    #print(rna)
    dna_fragments = []
    for i in range(len(rna)-3*len(peptide_sequence)+1):
       codon_idx = i
       potential_start = rna[codon_idx:codon_idx+3]
       if potential_start in first_peptide_rna_list:
          print("starting: " + str(i))
          print("x " + potential_start)
          fragment = ""
          peptide_idx = 0
          next_codon = potential_start
          next_peptide = peptide_sequence[peptide_idx]
          match_len = 0
          while next_codon in rna_by_protein[next_peptide] and match_len < len(peptide_sequence):
              fragment += next_codon
              match_len += 1
              #print("{0}{1}".format(" " * match_len, fragment))
              if codon_idx < len(rna)-6 and peptide_idx < (len(peptide_sequence)-1):
                  codon_idx += 3
                  peptide_idx += 1
                  next_codon = rna[codon_idx:codon_idx+3]
                  next_peptide = peptide_sequence[peptide_idx]
          if match_len == len(peptide_sequence):
              fragment = fragment.replace("U","T")
              if is_complement_strand:
                  fragment = utils.reverse_complement(fragment)
              dna_fragments.append(fragment)
          else:
              print("{0} {1}".format(match_len, fragment))

    return dna_fragments

def find_substrings_encoding_a_given_peptide_sequence(dna, peptide_sequence, rna_protein_table):
    if len(peptide_sequence) == 0 or len(dna) == 0:
        return []

    dna_complement = utils.reverse_complement(dna)
    #print(dna)
    #print(dna_complement)
    rna = dna.replace("T","U")
    rna_complement = dna_complement.replace("T","U")
    
    dna_fragments = []
    dna_frags = find_peptide_fragments_in_rna(rna, False, peptide_sequence, rna_protein_table)
    dna_fragments.extend(dna_frags)
    dna_frags = find_peptide_fragments_in_rna(rna_complement, True, peptide_sequence, rna_protein_table)
    dna_fragments.extend(dna_frags)
    

    #print("--")
    return dna_fragments

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: one param, filename, format:<str_dna>\n<str_peptide>.") 

    else:
        start = time.process_time()

        with open(sys.argv[1]) as f:
            lines = [line.rstrip() for line in f]

        with open("./RNA_codon_table_1.txt") as f:
            rna_protein_table = [line.rstrip() for line in f]

        dna = lines[0]
        for i in range(1,len(lines)-1):
            dna += lines[i]
        peptide_string = lines[len(lines)-1]

        dnas = find_substrings_encoding_a_given_peptide_sequence(dna, peptide_string, rna_protein_table)
        for dna in dnas:
            print(dna)


        end = time.process_time()
        #print("Time: {0}".format(end-start))

