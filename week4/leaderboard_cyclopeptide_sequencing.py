#! /bin/python3

import sys
import time
import itertools
from typing import cast


'''
This ended up as an attempt to answer problem 1.1.10 Excercise Break - Run LeaderboardCyclopeptideSequencing on Spectrum(25) with N=1000.
You should find 38 linear peptides of max score 83 (corresponding to 15 different cyclic peptides).

I never got 38 peptides, I got 10. Also the max observed score was 84.
97-147-113-128-99-163-128-114-147-71-115
97-147-113-128-99-163-128-114-147-115-71
97-147-113-128-99-163-57-71-114-147-71-115
97-147-113-128-99-163-57-71-114-147-115-71
97-147-113-128-99-163-71-57-114-147-71-115
97-147-113-128-99-163-71-57-114-147-115-71
97-147-113-128-99-163-128-57-57-147-71-115
97-147-113-128-99-163-128-57-57-147-115-71
97-147-113-128-99-163-71-57-57-57-147-71-115
97-147-113-128-99-163-71-57-57-57-147-115-71

I could have written the scoring function wrong, still not 100% confident about those loops.

I also wonder about the search space described by expand/trim, it seems that other answers should have been found, the same peptide sequences but rotated around. But they weren't found because "PFLQVY" was profitable enough to eliminate them.
A competing sequence would have to contain (either linearly or ciclically) PLFQVY in order to stay in the running.

I think the answer "38 linear peptides" is a result of a particular search execution path, not part of the requirements. Mine doesn't take that path, and so doesn't come up with the same result.
'''

_amino_acid_by_mass_lookup_ = {}
_mass_by_amino_acid_lookup_ = {}

_amino_acid_masses_ = []
_all_available_amino_acids_ = []

_debug_int_ = 0

class PeptideInfo:

    def __init__(self, pep_str):
        self.peptide_string = pep_str
        self.linear_spec = sorted(linear_theoretic_spec(pep_str))
        self.circular_spec = sorted(circ_theoretic_spec(pep_str))
        self.full_mass = full_peptide_mass(pep_str)
        self.linear_score = -1
        self.circular_score = -1
        self.debug_circular_score = -1
        self.slice_idx = 100

    def score_against_experimental_spectrum(self, exp_spec_weights):


        self.linear_score = match_spec_weights(self.linear_spec, exp_spec_weights)
        self.circular_score = match_spec_weights(self.circular_spec, exp_spec_weights)
        self.debug_circular_score = match_spec_weights(self.circular_spec[0:self.slice_idx], exp_spec_weights[0:self.slice_idx])
        if self.peptide_string == 'PFLQVYAGGGFDA' or self.peptide_string == 'PADFNQYVQLF':
            print(self.peptide_string)
            print(self.linear_score)
            print(" ".join([str(w) for w in self.linear_spec]))
            print(self.circular_score)
            print(" ".join([str(w) for w in self.circular_spec]))
            print(" ".join([str(w) for w in exp_spec_weights]))
            linSpectrumDict = {}
            for w in self.linear_spec:
                linSpectrumDict[w] = linSpectrumDict.get(w, 0) + 1
            theoSpectrumDict = {}
            for w in self.circular_spec:
                theoSpectrumDict[w] = theoSpectrumDict.get(w, 0) + 1
            spectrumDict = {}
            for w in exp_spec_weights:
                spectrumDict[w] = spectrumDict.get(w, 0) + 1
            print("242... t:{0} s:{1}".format(str(theoSpectrumDict[242]), str(spectrumDict[242])))
            print(linSpectrumDict)
            print(theoSpectrumDict)
            print(spectrumDict)

    def printout(self):
        print("{0} {1}".format(self.peptide_string, str(self.full_mass)))
        print("L {0} {1}".format(str(self.linear_score), ",".join([str(i) for i in self.linear_spec])))
        print("C {0} {1}".format(str(self.circular_score), ",".join([str(i) for i in self.circular_spec])))
        print("D {0} {1}".format(str(self.debug_circular_score), ",".join([str(i) for i in self.circular_spec[0:self.slice_idx]])))

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
    cumulative_mass_peptide_string = peptide_string[:]
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
    if peptide_string == 'PFLQVYAGGGFDA' or peptide_string == 'PADFNQYVQLF':
        print(str(len(peptide_string)) + " " + " ".join([str(m) for m in cumulative_mass_list]))
    for length in range(1, len(peptide_string)):
        upper_bound = 0
        if find_circular:
            upper_bound = len(peptide_string)
        else:
            upper_bound = len(peptide_string) - length + 1
        for idx in range(0, upper_bound):
             new_mass = cumulative_mass_list[idx+length] - cumulative_mass_list[idx]
             masses.append(new_mass)
             if peptide_string == 'PFLQVYAGGGFDA':
                 if new_mass in [97, 244, 485, 584, 747, 818, 875, 932, 989, 1136, 1251]:
                     print("PFLQ.... {0} len: {1} idx: {2}  {3}-{4}".format(str(new_mass), str(length), str(idx), str(cumulative_mass_list[idx+length]), str(cumulative_mass_list[idx])))
    masses.append(cumulative_mass_list[len(peptide_string)])
    return masses

def linear_theoretic_spec(peptide_string):
    return calc_theoretic_spectrum(peptide_string, False)

def circ_theoretic_spec(peptide_string):
    return calc_theoretic_spectrum(peptide_string, True)

def match_spec_weights(pep_spec_weights, experimental_spec_weights):
    peptide_weights = pep_spec_weights[:]
    experimental_weights = experimental_spec_weights[:]
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
   
    return score


def score_spectrum_against_weights(peptide, spectrum, find_circular = True):
    peptide_weights = []
    if find_circular:
        peptide_weights = circ_theoretic_spec(peptide)
    else:
        peptide_weights = linear_theoretic_spec(peptide)

    return match_spec_weights(peptide_weights, spectrum)

def match_circ_spectrum(peptide, spectrum):
    return score_spectrum_against_weights(peptide, spectrum, True)

def match_linear_spectrum(peptide, spectrum):
    return score_spectrum_against_weights(peptide, spectrum, False)

def expand_search(candidate_peptides, spectrum):
    new_candidate_peptides = {}
    for base_peptide in candidate_peptides.values():
        for amino_acid in _all_available_amino_acids_:
            if amino_acid != "I" and amino_acid != "K":
                new_candidate_peptide_str = base_peptide.peptide_string + amino_acid
                pi = PeptideInfo(new_candidate_peptide_str)
                pi.score_against_experimental_spectrum(spectrum)
                new_candidate_peptides[new_candidate_peptide_str] = pi

                #new_candidate_peptide_str = amino_acid + base_peptide.peptide_string
                #pi = PeptideInfo(new_candidate_peptide_str)
                #pi.score_against_experimental_spectrum(spectrum)
                #new_candidate_peptides[new_candidate_peptide_str] = pi

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
    unsorted_scores = []
    for pi in candidate_peptides.values():
        unsorted_scores.append(pi.linear_score)
    scores = sorted(unsorted_scores, reverse=True)

    if __debug__:
        if len(scores) > n:
            print("cutoff: {0}".format(str(scores[n-1])))
        else:
            print("no cutoff")

    new_candidate_peptides = candidate_peptides
    if len(scores) > n:
        cutoff = scores[n-1]
        #if _debug_int_ % 10 == 0:
        #    print(cutoff)
        new_candidate_peptides = {key:value for (key,value) in candidate_peptides.items() 
                                  if value.linear_score >= cutoff}
    
    print("dropped: {0}".format(str(len(candidate_peptides) - len(new_candidate_peptides))))
    return new_candidate_peptides
    

def leaderboard_cyclopeptide_sequencing(spectrum, cutoff, integer_mass_table):
    init_int_mass_tables(integer_mass_table)

    candidate_peptides = {}
    candidate_peptides[""] = PeptideInfo("")
    leader_peptides = []
    leader_peptide_score = 0
    parent_mass = spectrum[len(spectrum)-1]

    global _debug_int_
    _debug_int_ += 1

    while len(candidate_peptides) > 0:
        candidate_peptides = expand_search(candidate_peptides, spectrum)
        peps_to_remove = []

        for pi in candidate_peptides.values():
            if pi.full_mass == parent_mass:
                if pi.circular_score > leader_peptide_score:
                    print("new high score: {0}".format(str(pi.circular_score)))
                    leader_peptides = []
                    leader_peptide_score = pi.circular_score
                if pi.circular_score >= leader_peptide_score:
                    print("adding {0}".format(pi.peptide_string))
                    print(" ".join([str(i) for i in sorted(pi.linear_spec)]))
                    leader_peptides.append(pi)
                #peps_to_remove.append(pi.peptide_string)
            elif pi.full_mass > parent_mass:
                peps_to_remove.append(pi.peptide_string)
        print("candidates: {0}, leaders: {1}, leader_score: {2}, removing: {3}".format(str(len(candidate_peptides)), str(len(leader_peptides)), str(leader_peptide_score), str(len(peps_to_remove))))
        for i in peps_to_remove:
            candidate_peptides.pop(i, None)
        candidate_peptides = trim_candidate_peptides_list(candidate_peptides, cutoff)

        #if _debug_int_ % 10 == 0:
        #    print(" ".join(candidate_peptides))

    return leader_peptides                

def convert_amino_acid_strs_to_weights(sequences):
    aa_strs = []
    for s in sequences:
        print("-" + s.peptide_string + "-")
        #weights_str = "-".join([str(c) for c in s.linear_spec])
        weights_str = "-".join([str(_mass_by_amino_acid_lookup_[c]) for c in s.peptide_string])
        aa_strs.append(weights_str)
    return aa_strs

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: one param, filename, format:\n<int_leaderboard_cutoff>\n<list_int_weights>.") 

    else:
        start = time.process_time()

        with open(sys.argv[1]) as f:
            cutoff = int(f.readline())
            spectrum_line = f.readline()

        spectrum = [int(i) for i in spectrum_line.split(" ")]

        with open("./integer_mass_table.txt") as f:
            integer_mass_table = [line.rstrip() for line in f]

        sequences = leaderboard_cyclopeptide_sequencing(spectrum, cutoff, integer_mass_table)
        #leader = convert_amino_acid_strs_to_weights([sequences])[0]
        #print(leader)
        seq_as_weights = convert_amino_acid_strs_to_weights(sequences)
        for s in seq_as_weights:
            print(s)

        end = time.process_time()
        print("Time: {0}".format(end-start))

