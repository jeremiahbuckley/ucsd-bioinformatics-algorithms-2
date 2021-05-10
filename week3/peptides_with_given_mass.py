#! /bin/python3

import sys
import time
import array

_peptide_by_mass_lookup_ = {}


_peptide_str_count_at_mass_idx_ = array.array('l', [0 for i in range(5000)])
_cache_p_res_at_m_idx_ = array.array('i', [0 for i in range(5000)])
_peptide_strings_at_mass_idx_ = []
for i in range(5000):
    _peptide_strings_at_mass_idx_.append([])

_peptide_masses_ = []

def init_int_mass_tables(mass_table):
    mass_by_peptide = {}
    for kvp in mass_table:
        key = kvp[0:kvp.index(" ")]
        value = int(kvp[kvp.index(" ")+1:len(kvp)])
        mass_by_peptide[key] = value
        if key != "I" and key != "K":
            _peptide_by_mass_lookup_[value] = key
            _peptide_strings_at_mass_idx_[value].append(key)
            _peptide_str_count_at_mass_idx_[value] = 1
            _cache_p_res_at_m_idx_[value] = 1
            _peptide_masses_.append(value)
    #print(mass_by_peptide)
    #print(_peptide_by_mass_lookup_)

    #_peptide_str_count_at_mass_idx_[113] += 1 #I
    #_peptide_str_count_at_mass_idx_[128] += 1 #K
    _peptide_str_count_at_mass_idx_[114] += 1 #GG, GG (reversed)
    _peptide_str_count_at_mass_idx_[128] += 2 #GA,AG
    _peptide_str_count_at_mass_idx_[156] += 2 #GV, VG
    _peptide_str_count_at_mass_idx_[186] += 6 #EG,GE, SV, VS, AD, DA
    _peptide_strings_at_mass_idx_[114].append("GG")
    _peptide_strings_at_mass_idx_[128].append("AG")
    _peptide_strings_at_mass_idx_[128].append("GA")
    _peptide_strings_at_mass_idx_[156].append("GV")
    _peptide_strings_at_mass_idx_[156].append("VG")
    _peptide_strings_at_mass_idx_[186].append("EG")
    _peptide_strings_at_mass_idx_[186].append("GE")
    _peptide_strings_at_mass_idx_[186].append("SV")
    _peptide_strings_at_mass_idx_[186].append("VS")
    _peptide_strings_at_mass_idx_[186].append("AD")
    _peptide_strings_at_mass_idx_[186].append("DA")
    #print(_peptide_strings_at_mass_idx_)
    #print("z")
    #print(_peptide_str_count_at_mass_idx_)
    #print("z")
    return mass_by_peptide, _peptide_by_mass_lookup_

'''
idx: 0 1 2 3 4
mas: 1 2 3 4 5

cmi: 0 1 2 3  4  5
cmv: 0 1 3 6 10 15 
'''


def count_recursive_find_peptides_old (mass, known_results_cache, debug_indent):
    #_peptide_str_count_at_mass_idx_ = array.array('l', [0 for i in range(5000)])
    #_cache_p_res_at_m_idx_ = array.array('i', [0 for i in range(5000)])
    masses = []
    #print("{0}--> {1}.".format(" " * debug_indent, str(mass)))

    if _cache_p_res_at_m_idx_[mass] == 1:
        found_results = _peptide_str_count_at_mass_idx_[mass]
        #print("{0}cache: {1} {2}.".format(" " * debug_indent, str(mass), ",".join(known_results_cache[mass])))
        return found_results

    successful_finds = 0   
    for n_mass in _peptide_masses_:
        found_results = count_recursive_find_peptides(n_mass, known_results_cache, debug_indent+1)
        successful_finds += found_results

    if _cache_p_res_at_m_idx_[mass] == 1:
        raise RuntimeError("Cache miss. Cache should have been empty here. {0}.".format(mass))
    _peptide_str_count_at_mass_idx_[mass] = successful_finds
    _cache_p_res_at_m_idx_[mass] = 1
    return successful_finds


# track_peptide_strings is generally for debugging. on my system it will cause mem/cpu to spike a lot harder and cause failure when mass is somewhere  > 800
def count_recursive_find_peptides(mass, track_peptide_strings=False):
    recursion_stack = []
    recursion_stack.append([0, mass, 0, [""], []])
    recursion_stack.append([0, mass, 0, [""], []])
    while len(recursion_stack) > 1:
        stack_frame = recursion_stack.pop()
        previous_stack_frame_idx = len(recursion_stack)-1

        current_pep_idx = stack_frame[0]
        mass = stack_frame[1]
        result = stack_frame[2]
        pep_prefixes = stack_frame[3]
        pep_strs = stack_frame[4]

        if _cache_p_res_at_m_idx_[mass] == 1:
            result = _peptide_str_count_at_mass_idx_[mass]
            found_suffixes = _peptide_strings_at_mass_idx_[mass]

            if track_peptide_strings:
                for ps in pep_prefixes:
                    for ss in found_suffixes:
                        #print("w 1" + ps+ " 2" + ss)
                        recursion_stack[previous_stack_frame_idx][4].append(ps+ss)
            recursion_stack[previous_stack_frame_idx][2] += result
        elif mass < 57:
            recursion_stack[previous_stack_frame_idx][2] += 0
        elif current_pep_idx < len(_peptide_masses_):
            stack_frame[0] = current_pep_idx+1
            n_mass = mass-_peptide_masses_[current_pep_idx]
            prefix = _peptide_by_mass_lookup_[_peptide_masses_[current_pep_idx]]

            recursion_stack.append(stack_frame)
            recursion_stack.append([0, n_mass, 0, [prefix], []])
        else:
            _peptide_str_count_at_mass_idx_[mass] = result
            _cache_p_res_at_m_idx_[mass] = 1
            _peptide_strings_at_mass_idx_[mass] = pep_strs

            if track_peptide_strings:
                for ps in pep_prefixes:
                    for ss in pep_strs:
                        recursion_stack[previous_stack_frame_idx][4].append(ps+ss)
            recursion_stack[previous_stack_frame_idx][2] += result
    
    print("end: " + ",".join(recursion_stack[0][4]))        
    return recursion_stack[0][2]

def count_peptides_of_given_mass(mass):
    counts = count_recursive_find_peptides(mass)
 
    return counts

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: one param <int_mass>.") 

    else:
        start = time.process_time()

        mass = int(sys.argv[1])

        int_mass_table = []
        with open("./integer_mass_table.txt") as f:
            int_mass_table = [line.rstrip() for line in f]
        #print(int_mass_table)

        init_int_mass_tables(int_mass_table)
        masses = count_recursive_find_peptides(mass, mass < 800)
        #for mass in masses:
        #    print(mass)
        #print(len(masses))
        #print(len(masses) + 1) # to account for 0
        print(masses)
        print(masses + 1) # to account for 0

        end = time.process_time()
        print("Time: {0}".format(end-start))
