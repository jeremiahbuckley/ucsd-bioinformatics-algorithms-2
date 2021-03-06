#! /bin/python3

import sys
import time
import itertools


_amino_acid_by_mass_lookup_ = {}
_mass_by_amino_acid_lookup_ = {}

_amino_acid_masses_ = []
_all_available_amino_acids_ = []

_debug_int_ = 0

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

def full_peptide_mass(peptide_string):
    mass = 0
    for amino_acid in peptide_string:
        mass += _mass_by_amino_acid_lookup_[amino_acid]
    return mass

def calc_theoretic_spectrum(peptide_string, find_circular=True):
    masses = []
    cumulative_mass_peptide_string = peptide_string
    if find_circular:
        cumulative_mass_peptide_string += peptide_string #double so iterating through the loops is easy

    global _debug_int_
    _debug_int_ += 1
    #if _debug_int_ % 100 == 0:
    #    print(cumulative_mass_peptide_string)
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
        for idx in range(0, len(peptide_string)+1-upper_bound_offset):
             new_mass = cumulative_mass_list[idx+length] - cumulative_mass_list[idx]
             masses.append(new_mass)
    masses.append(cumulative_mass_list[len(peptide_string)])
    return masses
'''
2
0 1 2 3 4 5
0 1 12 123 1234 12345
0 1 12 123 1234 12345 123451 1234512 12345123 123451234 1234512345

12-0 123-1  1234-12 12345-123 123451-1234 1234512-12345
2-0  3-1    4-2     5-3       6-4         7-5
0+len(substr) -to- len(str)+len(substr)
0             -to- len(str)

12-0 123-1 1234-12 12345-123
2-0  3-1   4-2     5-3
0+len(substr) -to- len(str)
0             -to- len(str)-len(substr)



3
123-0 1234-1 12345-12 123451-123 1234512-1234

4
1234-0 12345-1 123451-12 1234512-123 12345123-1234
'''
def linear_theoretic_spec(peptide_string):
    return calc_theoretic_spectrum(peptide_string, False)

def circ_theoretic_spec(peptide_string):
    return calc_theoretic_spectrum(peptide_string, True)

def score_spectrum_against_weights(peptide, spectrum, find_circular = True):
    peptide_weights = []
    if find_circular:
        peptide_weights = circ_theoretic_spec(peptide)
    else:
        peptide_weights = linear_theoretic_spec(peptide)

    experimental_weights = spectrum[:]
    idx = 0
    score = 0
    while idx < len(peptide_weights):
        w = peptide_weights[idx]
        if w in experimental_weights and w in peptide_weights:
            experimental_weights.remove(w)
            peptide_weights.remove(w)
            score += 1
            # don't increment idx
        else:
            idx += 1
   
    #score -= (len(peptide_weights) * .5) 
    if score >= 83:
        print("x {0} {1}".format(str(score), "-".join([str(i) for i in peptide_weights])))
    
    return score

def match_circ_spectrum(peptide, spectrum):
    return score_spectrum_against_weights(peptide, spectrum, True)

def match_linear_spectrum(peptide, spectrum):
    return score_spectrum_against_weights(peptide, spectrum, False)

def expand_search(candidate_peptides, spectrum):
    new_candidate_peptides = {}
    for base_peptide in candidate_peptides:
        for amino_acid in _all_available_amino_acids_:
            if amino_acid != "I" and amino_acid != "K":
                new_candidate_peptide = base_peptide + amino_acid
                #new_candidate_peptides[new_candidate_peptide] = match_linear_spectrum(new_candidate_peptide, spectrum)
                new_candidate_peptides[new_candidate_peptide] = match_circ_spectrum(new_candidate_peptide, spectrum)
    return new_candidate_peptides

def spectrums_match(peptide, spectrum):
    peptide_circ_theoretic_spec = circ_theoretic_spec(peptide)
    target_spec = spectrum[:]
    #print("x")
    #print(" ".join([str(s) for s in peptide_circ_theoretic_spec]))
    #print(" ".join([str(s) for s in target_spec]))
    matches = []
    for mass in peptide_circ_theoretic_spec:
        if mass in target_spec:
            matches.append(mass)
        else:
            return False

    for mass in matches:
        peptide_circ_theoretic_spec.remove(mass)
        if mass in target_spec:
            target_spec.remove(mass)
        else:
            return False # duplicates in peptide_circ_theoretic_spec but not in target_spec

    return len(peptide_circ_theoretic_spec) == 0 and \
           len(target_spec) == 0

def pep_consistent_with_spectrum(peptide, spectrum):
    lin_theor_spec = linear_theoretic_spec(peptide)
    target_spec = spectrum[:]
    remove_masses = []
    for m in lin_theor_spec:
        # if there's a mass in lin_theo that is not in target
        if m not in target_spec:
            return False
        remove_masses.append(m)
    for rm in remove_masses:
        # if there's a mass in lin_theo that shows up more times than in target
        if rm not in lin_theor_spec:
            return False
        lin_theor_spec.remove(rm)
    return True

def trim_candidate_peptides_list(candidate_peptides, n):
    trimmed_candidate_peptides = {}
    scores = [s for s in sorted(candidate_peptides.values(), reverse=True)]
    global _debug_int_
    if _debug_int_ % 10 == 0:
        print(",".join([str(i) for i in scores]))
    new_candidate_peptides = candidate_peptides
    if len(scores) > n:
        cutoff = scores[n]
        if _debug_int_ % 10 == 0:
            print(cutoff)
        new_candidate_peptides = {key:value for (key,value) in candidate_peptides.items() 
                                  if value >= cutoff}
    return new_candidate_peptides
    


def leaderboard_cyclopeptide_seqencing(spectrum, cutoff, integer_mass_table):
    init_int_mass_tables(integer_mass_table)

    candidate_peptides = {}
    candidate_peptides[""] = 0
    leader_peptides = []
    leader_peptide_score = 0
    parent_mass = spectrum[len(spectrum)-1]

    global _debug_int_
    _debug_int_ += 1

    while len(candidate_peptides) > 0:
        candidate_peptides = expand_search(candidate_peptides, spectrum)
        peps_to_remove = []
        for kvp in candidate_peptides.items():
            pep = kvp[0]
            pep_score = kvp[1]
            pep_full_mass = full_peptide_mass(pep)
            if pep_full_mass == parent_mass:
                if pep_score > leader_peptide_score:
                    print("new")
                    leader_peptides = []
                    leader_peptide_score = pep_score
                if pep_score == leader_peptide_score:
                    print("adding {0}".format(pep))
                    leader_peptides.append(pep)
                peps_to_remove.append(pep)
            elif pep_full_mass > parent_mass:
                peps_to_remove.append(pep)
        for i in peps_to_remove[::-1]:
            candidate_peptides.pop(i, None)
        candidate_peptides = trim_candidate_peptides_list(candidate_peptides, cutoff)

        if _debug_int_ % 10 == 0:
            print(" ".join(candidate_peptides))

    return leader_peptides                

def convert_amino_acid_strs_to_weights(sequences):
    aa_strs = []
    for s in sequences:
        print("-" + s + "-")
        weights_str = "-".join([str(_mass_by_amino_acid_lookup_[c]) for c in s])
        aa_strs.append(weights_str)
    return aa_strs

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: one param, filename, format:<int_cutoff>\\n<list_int_weights>.") 

    else:
        start = time.process_time()

        with open(sys.argv[1]) as f:
            cutoff = int(f.readline())
            spectrum_line = f.readline()

        spectrum = [int(i) for i in spectrum_line.split(" ")]

        with open("./integer_mass_table.txt") as f:
            integer_mass_table = [line.rstrip() for line in f]

        sequences = leaderboard_cyclopeptide_seqencing(spectrum, cutoff, integer_mass_table)
        #leader = convert_amino_acid_strs_to_weights([sequences])[0]
        #print(leader)
        seq_as_weights = convert_amino_acid_strs_to_weights(sequences)
        for s in seq_as_weights:
            print(s)

        end = time.process_time()
        print("Time: {0}".format(end-start))

