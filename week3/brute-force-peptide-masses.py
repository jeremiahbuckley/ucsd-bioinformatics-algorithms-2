#! /bin/python3

import sys
import time
import itertools
import array

_peptide_by_mass_lookup_ = {}
_mass_by_peptide_lookup_ = {}

_peptides_at_mass_idx_ = array.array('l', [0 for i in range(5000)])
_cache_p_res_at_m_idx_ = array.array('i', [0 for i in range(5000)])
_peptide_masses_ = []

def init_int_mass_tables(mass_table):
    for kvp in mass_table:
        key = kvp[0:kvp.index(" ")]
        value = int(kvp[kvp.index(" ")+1:len(kvp)])
        _mass_by_peptide_lookup_[key] = value
        if key != "I" and key != "K":
            _peptide_by_mass_lookup_[value] = key
            _peptides_at_mass_idx_[value] = 1
            _cache_p_res_at_m_idx_[value] = 1
            _peptide_masses_.append(value)
    print(_mass_by_peptide_lookup_)
    print(_peptide_by_mass_lookup_)

def brute_force_generate_peptide_combos(iterative_count, int_mass_table):
    init_int_mass_tables(int_mass_table)

    # G 57, *5 = 285
    # W 186 *5 = 930
    count_per_mass = {}
    for rep_ct in range(1,iterative_count+1):
        for pep_str in itertools.product('GASPVTCINDKEMHFRYW', repeat = rep_ct):
            m = 0
            for c in pep_str:
                m += _mass_by_peptide_lookup_[c]

            if m not in count_per_mass:
                count_per_mass[m] = []
            count_per_mass[m].append("".join(pep_str))

    write_result(count_per_mass)

def write_result(count_per_mass):
    with open("./count_per_mass.txt", "a") as f:
        for kvp in count_per_mass.items():
            if kvp[0] <= (57 * iterative_count): # G=57, if kvp[0] <= GGGGG then it has reached max-count at that weight
                #print(kvp[0])
                f.write("{0}    {1}\n".format(kvp[0],len(kvp[1])))
                if kvp[0] == 257 or kvp[0] == 273 or kvp[0] == 285 or len(kvp[1]) < 3:
                    print(str(kvp[0]) + " " + ",".join(kvp[1]))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: one param <int_mass>.") 

    else:
        sys.setrecursionlimit(12000)
        start = time.process_time()

        iterative_count = int(sys.argv[1])

        int_mass_table = []
        with open("./integer_mass_table.txt") as f:
            int_mass_table = [line.rstrip() for line in f]
        #print(int_mass_table)

        brute_force_generate_peptide_combos(iterative_count, int_mass_table)

        end = time.process_time()
        print("Time: {0}".format(end-start))
