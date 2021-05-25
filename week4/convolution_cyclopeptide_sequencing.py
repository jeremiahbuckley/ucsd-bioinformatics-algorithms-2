#! /bin/python3

import sys
import time
import itertools


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

_debug_int_ = 0

_peptide_weights_set_ = []

class PeptideInfo:

    def __init__(self, pep_weights):
        self.peptide_weights = pep_weights[:]
        self.peptide_string = self.make_pep_str(self.peptide_weights)
        self.linear_spec = linear_theoretic_spec(self.peptide_weights)
        self.circular_spec = circ_theoretic_spec(self.peptide_weights)
        self.full_mass = full_peptide_mass(self.peptide_weights)

    def score_against_experimental_spectrum(self, exp_spec_weights):
        self.linear_score = match_spec_weights(self.linear_spec, exp_spec_weights)
        self.circular_score = match_spec_weights(self.circular_spec, exp_spec_weights)

    def make_pep_str(self, peptide_weights):
        return ("-".join([str(i) for i in peptide_weights]))


    def printout(self):
        print("{0} {1}".format(self.peptide_string, str(self.full_mass)))
        print("L {0} {1}".format(str(self.linear_score), "-".join([str(i) for i in self.linear_spec])))
        print("C {0} {1}".format(str(self.circular_score), "-".join([str(i) for i in self.circular_spec])))

def full_peptide_mass(peptide_weights):
    mass = 0
    for w in peptide_weights:
        mass += w
    return mass

def calc_theoretic_spectrum(peptide_weights, find_circular=False):
    masses = []
    #print(len(peptide_weights))
    cumulative_mass_peptide_weights = peptide_weights[:]
    if find_circular:
        cumulative_mass_peptide_weights += peptide_weights #double so iterating through the loops is easy

    global _debug_int_
    _debug_int_ += 1
    #if _debug_int_ % 100 == 0:
    #    print(cumulative_mass_peptide_weights)
    cumulative_mass_list = []
    cumulative_mass = 0
    cumulative_mass_list.append(cumulative_mass)
    for w in cumulative_mass_peptide_weights:
        cumulative_mass += w
        cumulative_mass_list.append(cumulative_mass)

    masses.append(0)
    # intentionally do not add full-length substring in this loop
    #print("cts: {0} {1}".format(str(len(peptide_weights)), "-".join([str(i) for i in peptide_weights])))
    for length in range(1, len(peptide_weights)):
        upper_bound_offset = 0
        if not find_circular:
            upper_bound_offset = length
        for idx in range(0, len(peptide_weights)+1-upper_bound_offset):
             new_mass = cumulative_mass_list[idx+length] - cumulative_mass_list[idx]
             masses.append(new_mass)
        #print("    cts-work: {0} {1}".format(length, "-".join([str(i) for i in masses])))
    masses.append(cumulative_mass_list[len(peptide_weights)])
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
def linear_theoretic_spec(peptide_weights):
    return calc_theoretic_spectrum(peptide_weights, False)

def circ_theoretic_spec(peptide_weights):
    return calc_theoretic_spectrum(peptide_weights, True)

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
   
    #score -= (len(peptide_weights) * .5) 
    #if score >= 83:
    #    print("x {0} {1}".format(str(score), "-".join([str(i) for i in peptide_weights])))
    
    return score


def score_spectrum_against_weights(peptide_weights, spectrum, find_circular = True):
    theo_spec = []
    if find_circular:
        theo_spec = circ_theoretic_spec(peptide_weights)
    else:
        theo_spec = linear_theoretic_spec(peptide_weights)

    return match_spec_weights(theo_spec, spectrum)

def match_circ_spectrum(peptide, spectrum):
    return score_spectrum_against_weights(peptide, spectrum, True)

def match_linear_spectrum(peptide, spectrum):
    return score_spectrum_against_weights(peptide, spectrum, False)

def expand_search(candidate_peptides, spectrum):
    new_candidate_peptides = {}
    for base_peptide in candidate_peptides.values():
        for weight in _peptide_weights_set_:
            new_candidate_peptide_weights = base_peptide.peptide_weights + [weight]
            pi = PeptideInfo(new_candidate_peptide_weights)
            pi.score_against_experimental_spectrum(spectrum)
            new_candidate_peptides[pi.peptide_string] = pi

            new_candidate_peptide_weights = [weight] + base_peptide.peptide_weights
            pi = PeptideInfo(new_candidate_peptide_weights)
            pi.score_against_experimental_spectrum(spectrum)
            new_candidate_peptides[pi.peptide_string] = pi

    return new_candidate_peptides

def trim_candidate_peptides_list(candidate_peptides, n):
    trimmed_candidate_peptides = {}
    unsorted_scores = []
    for pi in candidate_peptides.values():
        unsorted_scores.append(pi.linear_score)
    scores = sorted(unsorted_scores, reverse=True)
    #global _debug_int_
    #if _debug_int_ % 10 == 0:
    #    print(",".join([str(i) for i in scores]))
    new_candidate_peptides = candidate_peptides
    if len(scores) > n:
        cutoff = scores[n-1]
        #if _debug_int_ % 10 == 0:
        #    print(cutoff)
        new_candidate_peptides = {key:value for (key,value) in candidate_peptides.items() 
                                  if value.linear_score >= cutoff}
    return new_candidate_peptides
    

def spectral_convolution_likely_weights(spectrum):
    cvs = {}
    print(spectrum)
    for i in range(len(spectrum)-1):
        for j in range(i+1,len(spectrum)):
            w = spectrum[j] - spectrum[i]
            if w != 0:
                if w not in cvs:
                    cvs[w] = 0
                cvs[w] += 1

    return cvs


def remove_out_of_band_values(cvs):
    # delete values < 57
    delete_keys = []
    for k in cvs.keys():
        if k < 57 or k > 200:
            delete_keys.append(k)

    for k in delete_keys:
        cvs.pop(k)
    print(cvs)
    return cvs

def find_top_convoluted_weights(spectrum, cutoff):
    convoluted_likely_weights = spectral_convolution_likely_weights(spectrum)
    convoluted_likely_weights = remove_out_of_band_values(convoluted_likely_weights)
    scores = []
    for value in convoluted_likely_weights.values():
        scores.append(value)
    sorted_scores = sorted(scores, reverse=True)

    global _peptide_weights_set_    
    _peptide_weights_set_ = sorted(convoluted_likely_weights.keys())

    print("all w scores: {0}".format("_".join([str(i) for i in sorted_scores])))
    if cutoff < len(convoluted_likely_weights):
        cutoff_weight = sorted_scores[cutoff]

        pws = [key for (key, value) in convoluted_likely_weights.items() if value >= cutoff_weight]
        _peptide_weights_set_ = sorted(pws)

    print(",".join([str(i) for i in _peptide_weights_set_]))
    return _peptide_weights_set_ 

def leaderboard_cyclopeptide_seqencing(spectrum, cutoff, convoluted_weights_cutoff):
    find_top_convoluted_weights(spectrum, convoluted_weights_cutoff)

    candidate_peptides = {}
    candidate_peptides[""] = PeptideInfo([])
    leader_peptides = []
    leader_peptide_score = 0
    parent_mass = spectrum[len(spectrum)-1]

    global _debug_int_

    while len(candidate_peptides) > 0:
        _debug_int_ += 1

        candidate_peptides = expand_search(candidate_peptides, spectrum)
        peps_to_remove = []
        for pi in candidate_peptides.values():
            if pi.full_mass == parent_mass:
                if pi.circular_score > leader_peptide_score:
                    print("new")
                    leader_peptides = []
                    leader_peptide_score = pi.circular_score
                if pi.circular_score >= leader_peptide_score:
                    print("adding {0}".format(pi.peptide_string))
                    leader_peptides.append(pi)
                #peps_to_remove.append(pi.peptide_string)
            elif pi.full_mass > parent_mass:
                peps_to_remove.append(pi.peptide_string)
        for i in peps_to_remove[::-1]:
            candidate_peptides.pop(i, None)
        candidate_peptides = trim_candidate_peptides_list(candidate_peptides, cutoff)

        #if _debug_int_ % 10 == 0:
        #    print(" ".join(candidate_peptides))

    return leader_peptides                

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: one param, filename, format:\n<int_convoluted_weights_cutoff>\n<int_leaderboard_cutoff>\n<list_int_weights>.") 

    else:
        start = time.process_time()

        with open(sys.argv[1]) as f:
            convoluted_weights_cutoff = int(f.readline())
            leaderboard_cutoff = int(f.readline())
            spectrum_line = f.readline()

        spectrum = [int(i) for i in spectrum_line.split(" ")]

        sequences = leaderboard_cyclopeptide_seqencing(spectrum, leaderboard_cutoff, convoluted_weights_cutoff)
        #leader = convert_amino_acid_strs_to_weights([sequences])[0]
        #print(leader)
        for s in sequences:
            print(s.peptide_string)

        end = time.process_time()
        print("Time: {0}".format(end-start))

