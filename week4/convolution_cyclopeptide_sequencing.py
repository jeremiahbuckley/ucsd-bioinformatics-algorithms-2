#! /bin/python3

import sys
import time
import itertools

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
    cumulative_mass_peptide_weights = peptide_weights[:]
    if find_circular:
        cumulative_mass_peptide_weights += peptide_weights #double so iterating through the loops is easy

    cumulative_mass_list = []
    cumulative_mass = 0
    cumulative_mass_list.append(cumulative_mass)
    for w in cumulative_mass_peptide_weights:
        cumulative_mass += w
        cumulative_mass_list.append(cumulative_mass)

    masses.append(0)
    # intentionally do not add full-length substring in this loop
    for length in range(1, len(peptide_weights)):
        upper_bound = len(peptide_string)
        if not find_circular:
            upper_bound -= (length - 1)
        for idx in range(0, upper_bound):
             new_mass = cumulative_mass_list[idx+length] - cumulative_mass_list[idx]
             masses.append(new_mass)
    masses.append(cumulative_mass_list[len(peptide_weights)])
    return masses

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

    return new_candidate_peptides

def trim_candidate_peptides_list(candidate_peptides, n):
    trimmed_candidate_peptides = {}
    unsorted_scores = []
    for pi in candidate_peptides.values():
        unsorted_scores.append(pi.linear_score)
    scores = sorted(unsorted_scores, reverse=True)
    new_candidate_peptides = candidate_peptides
    if len(scores) > n:
        cutoff = scores[n-1]
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

    while len(candidate_peptides) > 0:
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
        for s in sequences:
            print(s.peptide_string)

        end = time.process_time()
        print("Time: {0}".format(end-start))

