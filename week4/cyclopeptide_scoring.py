#! /bin/python3

import sys
import time
import math

_amino_acid_by_mass_lookup_ = {}
_mass_by_amino_acid_lookup_ = {}

_amino_acid_masses_ = []
_all_available_amino_acids_ = []


def init_int_mass_tables(mass_table):
    for kvp in mass_table:
        print(kvp)
        key = kvp[0:kvp.index(" ")]
        value = int(kvp[kvp.index(" ")+1:len(kvp)])
        _mass_by_amino_acid_lookup_[key] = value
        _all_available_amino_acids_.append(key)
        if key != "I" and key != "K":
            _amino_acid_by_mass_lookup_[value] = key
            _amino_acid_masses_.append(value)
    print(_mass_by_amino_acid_lookup_)
    print(_amino_acid_by_mass_lookup_)

def calc_theoretic_spectrum(peptide_string, find_circular=True):
    masses = []
    cumulative_mass_peptide_string = peptide_string
    if find_circular:
        cumulative_mass_peptide_string += peptide_string #double so iterating through the loops is easy

    cumulative_mass_list = []
    cumulative_mass = 0
    cumulative_mass_list.append(cumulative_mass)
    for amino_acid in cumulative_mass_peptide_string:
        cumulative_mass += _mass_by_amino_acid_lookup_[amino_acid]
        cumulative_mass_list.append(cumulative_mass)

    masses.append(0)
    # intentionally do not add full-length substring in this loop
    for length in range(1, len(peptide_string)):
        upper_bound_offset = 0
        if not find_circular:
            upper_bound_offset = length
        for idx in range(1, len(peptide_string)+1-upper_bound_offset):
             new_mass = cumulative_mass_list[idx+length] - cumulative_mass_list[idx]
             masses.append(new_mass)
    masses.append(cumulative_mass_list[len(peptide_string)])
    return masses

def linear_theoretic_spec(peptide_string):
    return calc_theoretic_spectrum(peptide_string, False)

def circ_theoretic_spec(peptide_string):
    return calc_theoretic_spectrum(peptide_string, True)

def score_theoretic_spectrum_against_weights(peptide, experimental_weights, int_mass_table):
    init_int_mass_tables(int_mass_table)
    peptide_weights = circ_theoretic_spec(peptide)

    idx = 0
    score = 0
    while idx < len(experimental_weights):
        w = experimental_weights[idx]
        if w in experimental_weights and w in peptide_weights:
            experimental_weights.remove(w)
            peptide_weights.remove(w)
            score += 1
            # don't increment idx
        else:
            idx += 1

    return score

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: one param, filename, format:<str_peptide>\\n<list_theoretic_weights>.") 

    else:
        start = time.process_time()

        peptide = ""
        weights = []
        with open(sys.argv[1]) as f:
            peptide = f.readline().rstrip()
            weights_str = f.readline()
            weights = [int(i) for i in weights_str.split(" ")]

        with open("./integer_mass_table.txt") as f:
            integer_mass_table = [line.rstrip() for line in f]

        score = score_theoretic_spectrum_against_weights(peptide, weights, integer_mass_table)
        print(score)


        end = time.process_time()
        print("Time: {0}".format(end-start))

