
import pandas as pd
import math
import io

def reverse_complement(kmer):
    reversi_dict = {"A":"T","C":"G","G":"C","T":"A"}
    reverse=""
    for i in range(len(kmer)-1,-1,-1):
            reverse += reversi_dict[kmer[i]]
    return reverse


def find_neigh_recursive(kmer, hamming_dist, neighbors, used_idxs, check_rev_comp):
    for i in range(len(kmer)):
        if i in used_idxs:
            pass
        else:
            for nucleotide in ["A","C","G","T"]:
                new_kmer = kmer[0:i]+nucleotide+kmer[i+1:len(kmer)]
                # print(new_kmer)
                if hamming_dist == 1:
                    neighbors.add(new_kmer)
                else:
                    fresh_used_idxs = used_idxs[:]
                    fresh_used_idxs.append(i)
                    find_neigh_recursive(new_kmer,
                                         hamming_dist-1,
                                         neighbors,
                                         fresh_used_idxs,
                                         check_rev_comp)
                if (check_rev_comp):
                    reverse = reverse_complement(new_kmer)
                    # print(reverse)
                    if hamming_dist == 1:
                        neighbors.add(reverse)
                    else:
                        fresh_used_idxs = used_idxs[:]
                        fresh_used_idxs.append(i)
                        find_neigh_recursive(reverse,
                                             hamming_dist-1,
                                             neighbors,
                                             fresh_used_idxs,
                                             check_rev_comp)
    return neighbors

def find_neighbors(kmer, hamming_dist, check_rev_comp):
    neighbors = find_neigh_recursive(kmer,
                                     hamming_dist,
                                     set(),
                                     [],
                                     check_rev_comp)
    return neighbors

def calc_hamming_distance(snippet1, snippet2):
    hd = 0
    #print("{0}\n{1}".format(snippet1, snippet2))
    for i in range(len(snippet1)):
        if snippet1[i] != snippet2[i]:
            hd += 1
    #print(hd)
    return hd

def create_dataframe_from_list_of_strings(list_of_strings, has_eol=True, header=False):
    if has_eol:
        one_string = "".join(list_of_strings)
    else:
        one_string = "\n".join(list_of_strings)
    str_csv = io.StringIO(one_string)
    if header:
        df = pd.read_csv(str_csv, sep=" ")
    else: 
        df = pd.read_csv(str_csv, sep=" ", header=None)
    return df

# min_count defaults to 0
# sometimes you might want to have min_count = 1, or .1, a "Laplace correction"
def calc_basic_score_from_motifs_list(motifs, min_count=0):
    if len(motifs) < 1:
        return 0

    nuc_count = create_motifs_count_from_motifs(motifs, min_count)
    score = len(motifs[0]) * len(motifs) - sum([max(col) for col in zip(*nuc_count)])
    return score

# min_count defaults to 0
# sometimes you might want to have min_count = 1, or .1, a "Laplace correction"
def create_motifs_count_from_motifs(motifs, min_count=0):
    if len(motifs) < 1:
        raise ValueError("Need at least 1 motif for count.")

    nuc_count = [[min_count]*len(motifs[0]) for i in range(4)]

    for col, colvalues in enumerate(zip(*motifs)):
        for row, letter in enumerate('ACGT'):
            nuc_count[row][col] += colvalues.count(letter)

    return nuc_count

def create_profile_from_motifs(motifs, min_count=0):
    if len(motifs) < 1:
        raise ValueError("Need at least 1 motif for count.")

    nuc_count = create_motifs_count_from_motifs(motifs, min_count)

    divisor = len(motifs)+ (min_count*4)
    for row in range(4):
        for col in range(len(motifs[0])):
            nuc_count[row][col] /= divisor

    return nuc_count

def create_profile_from_count_by_nuc_dict(count_by_nuc_dict, num_motifs, min_count=1):
    motif_len=len(count_by_nuc_dict["A"])

    profile = {"A":[0]*motif_len,"C":[0]*motif_len,"G":[0]*motif_len,"T":[0]*motif_len}

    divisor = num_motifs+min_count*4 #if min_count=1, then add 4 to the divisor
    for kvp in count_by_nuc_dict.items():
        for i in range(motif_len):
            profile[kvp[0]][i] = kvp[1][i] / divisor

    return profile

def create_consensus_motif_from_profile(profile):
    motif = ""
    motif_len = len(profile["A"])
    for i in range(motif_len):
        max_ct = profile["A"][i]
        max_char = "A"
        for key in ["C","G","T"]:
            if profile[key][i] > max_ct:
                max_ct = profile[key][i]
                max_char = key
        motif+=max_char
       
    return motif

def calculate_entropy_score_of_motif(motif, profile):
    nuc_row={"A":0, "C":1, "G":2, "T":3}
    score = 0
    for i in range(len(motif)):
        row_idx = nuc_row[motif[i]]
        prob = profile[row_idx][i]
        score += -1*(prob*math.log2(prob))
    return score

def calculate_probability_of_motif(motif, profile):
    nuc_row={"A":0, "C":1, "G":2, "T":3}
    score = 1
    for i in range(len(motif)):
        row_idx = nuc_row[motif[i]]
        prob = profile[row_idx][i]
        score *= prob
    return score

