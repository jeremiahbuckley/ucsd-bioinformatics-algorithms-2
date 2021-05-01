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
    print(_peptide_strings_at_mass_idx_)
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

def recursive_find_peptides_old(mass, known_results_cache, debug_indent):
    masses = []
    #print("{0}--> {1}.".format(" " * debug_indent, str(mass)))

    if mass in known_results_cache:
        #print("{0}cache: {1} {2}.".format(" " * debug_indent, str(mass), ",".join(known_results_cache[mass])))
        found_results = known_results_cache[mass]
        return found_results, len(found_results) > 0

    recursive_success = False
    successful_finds = []   
    for kvp in _peptide_by_mass_lookup_.items():
        n_pep = kvp[1]
        n_mass = mass-kvp[0]

        if n_mass < 0:
            pass # success = False, but it's already False by default
        elif n_mass == 0:
            #print("{0}zero : {1} {2}".format(" " * debug_indent, n_pep))
            recursive_success = True
            successful_finds.extend([n_pep])
        else:
            found_peptide_strings,success = recursive_find_peptides(n_mass, known_results_cache, debug_indent+1)
            if success:
                #print("{0}recrs: '{1}' {2}".format(" " * debug_indent, n_pep, ",".join(found_peptide_strings)))
                n_pot_pep_str = []
                for recurs_str in found_peptide_strings:
                    n_pot_pep_str.append(n_pep + recurs_str)
                successful_finds += n_pot_pep_str
                #print("{0}s_finds: {1} {2}.".format(" " * debug_indent, ",".join(successful_finds), ",".join(n_pot_pep_str)))
                recursive_success = True

    #f_set = set()
    #for find in successful_finds:
    #    f_set.add("".join(sorted(find)))
    #print("{0}f_set : {1} {2}.".format(" " * debug_indent, str(mass), f_set))
    #unique_successful_finds = list(f_set)
    unique_successful_finds = successful_finds

    if mass in known_results_cache:
        raise RuntimeError("Cache miss. Cache should have been empty here. {0}.".format(mass))
    else:
        if recursive_success:
            #print("{0} write success cache: {1} {2}.".format(" " * debug_indent, str(mass), ",".join(unique_successful_finds)))
            known_results_cache[mass] = unique_successful_finds
        else:
            #print("{0} write failure cache: {1}".format(" " * debug_indent, str(mass) ))
            known_results_cache[mass] = []

    return unique_successful_finds, recursive_success


def count_recursive_find_peptides_old(mass, known_results_cache, debug_indent):
    masses = []
    #print("{0}--> {1}.".format(" " * debug_indent, str(mass)))

    if mass in known_results_cache:
        #print("{0}cache: {1} {2}.".format(" " * debug_indent, str(mass), ",".join(known_results_cache[mass])))
        found_results = known_results_cache[mass]
        return found_results

    successful_finds = 0   
    for kvp in _peptide_by_mass_lookup_.items():
        #n_pep = kvp[1]
        n_mass = mass-kvp[0]
        found_results = count_recursive_find_peptides(n_mass, known_results_cache, debug_indent+1)
        successful_finds += found_results

    if mass in known_results_cache:
        raise RuntimeError("Cache miss. Cache should have been empty here. {0}.".format(mass))
    known_results_cache[mass] = successful_finds

    return successful_finds

def count_peptides_of_given_mass_old(mass, int_mass_table):
    mass_by_peptide = init_int_mass_tables(int_mass_table)
    counts = 0
    known_results_cache = {}
    for kvp in _peptide_by_mass_lookup_.items():
        known_results_cache[kvp[0]] = 1
    #    print("{0} {1} {2}".format(kvp[1], kvp[0]*2, kvp[0]*3))
    #for kvp in _peptide_by_mass_lookup_.items():
    #    for kvp2 in _peptide_by_mass_lookup_.items():
    #        print("{0} {1}".format(kvp[0] + kvp2[0], kvp[1] + kvp2[1]))
    #known_results_cache[113] += 1 #I
    #known_results_cache[128] += 1 #K
    known_results_cache[114] += 2 #GG,GG(reversed)
    known_results_cache[128] += 2 #GA,AG
    known_results_cache[156] += 2 #GV, VG
    known_results_cache[186] += 6 #EG,GE, AD, DA, SV, VS
    print("z")
    print(known_results_cache)
    print("z")
    counts = recursive_find_peptides(mass, known_results_cache, [])
 
    return counts

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

def count_recursive_find_peptides(mass, debug_indent):
    recursion_stack = []
    recursion_stack.append([0, mass, 0, [""], []])
    recursion_stack.append([0, mass, 0, [""], []])
    while len(recursion_stack) > 1:
        stack_frame_idx = len(recursion_stack)-1
        previous_stack_frame_idx = stack_frame_idx-1
        sf = recursion_stack[stack_frame_idx]
        #if stack_frame_idx == 4 or stack_frame_idx ==6:
        #    print("{0} [{1} {2} {3} {4} {5}]".format(len(recursion_stack), sf[0], sf[1], sf[2], ",".join(sf[3]), ",".join(sf[4])))
        pep = sf[0]
        mass = sf[1]
        result = sf[2]
        pep_prefixes = sf[3]
        pep_strs = sf[4]
        #if len(recursion_stack) == 3 and (pep_prefixes[0] == "A" or pep_prefixes[0] == "W"):
        #    print("!!!! " + pep_prefixes[0])

        if _cache_p_res_at_m_idx_[mass] == 1:
            #print("x"+str(mass))
            #print(_peptide_strings_at_mass_idx_)
            found_results = _peptide_str_count_at_mass_idx_[mass]
            found_suffixes = _peptide_strings_at_mass_idx_[mass]
            for ps in recursion_stack[stack_frame_idx][3]:
                for ss in found_suffixes:
                    #print("w 1" + ps+ " 2" + ss)
                    recursion_stack[previous_stack_frame_idx][4].append(ps+ss)
            recursion_stack[previous_stack_frame_idx][2] += found_results
            recursion_stack.pop(stack_frame_idx)
        elif mass < 57:
            recursion_stack[previous_stack_frame_idx][2] += 0
            recursion_stack.pop(stack_frame_idx)
        elif pep < len(_peptide_masses_):
            n_mass = mass-_peptide_masses_[pep]
            recursion_stack[stack_frame_idx][0] = pep+1
            prefix = _peptide_by_mass_lookup_[_peptide_masses_[pep]]
            #print("z" + prefix)
            recursion_stack.append([0, n_mass, 0, [prefix], []])
        else:
            _peptide_str_count_at_mass_idx_[mass] = result
            _cache_p_res_at_m_idx_[mass] = 1
            #print("y" + ",".join(pep_strs))
            _peptide_strings_at_mass_idx_[mass] = pep_strs
            for ps in recursion_stack[stack_frame_idx][3]:
                for ss in pep_strs:
                    recursion_stack[previous_stack_frame_idx][4].append(ps+ss)
            recursion_stack[previous_stack_frame_idx][2] += result
            recursion_stack.pop(stack_frame_idx)
    
    print("end: " + ",".join(recursion_stack[0][4]))        
    return recursion_stack[0][2]

def count_peptides_of_given_mass(mass):
    #_peptide_str_count_at_mass_idx_ = array.array('l', [0 for i in range(5000)])
    #_cache_p_res_at_m_idx_ = array.array('i', [0 for i in range(5000)])
    counts = 0
    #    print("{0} {1} {2}".format(kvp[1], kvp[0]*2, kvp[0]*3))
    #for kvp in _peptide_by_mass_lookup_.items():
    #    for kvp2 in _peptide_by_mass_lookup_.items():
    #        print("{0} {1}".format(kvp[0] + kvp2[0], kvp[1] + kvp2[1]))
    counts = count_recursive_find_peptides(mass, 1)
 
    return counts

def peptides_of_given_mass(mass, int_mass_table):
    mass_by_peptide = init_int_mass_tables(int_mass_table)
    masses = []
    known_results_cache = {}
    for kvp in _peptide_by_mass_lookup_.items():
        known_results_cache[kvp[0]] = [kvp[1]]
    #    print("{0} {1} {2}".format(kvp[1], kvp[0]*2, kvp[0]*3))
    #for kvp in _peptide_by_mass_lookup_.items():
    #    for kvp2 in _peptide_by_mass_lookup_.items():
    #        if (kvp[0] +kvp2[0]) < 170 or ((kvp[0] + kvp2[0]) < 190 and (kvp[0] + kvp2[0] > 180)):
    #            print("{0} {1}".format(kvp[0] + kvp2[0], kvp[1] + kvp2[1]))
    #known_results_cache[113].append("I")
    #kinown_results_cache[128].append("K")
    known_results_cache[114].append("GG")
    #known_results_cache[114].append("GG")
    known_results_cache[128].append("GA")
    known_results_cache[128].append("AG")
    known_results_cache[156].append("GV")
    known_results_cache[156].append("VG")
    known_results_cache[186].append("AD")
    known_results_cache[186].append("DA")
    known_results_cache[186].append("EG")
    known_results_cache[186].append("GE")
    known_results_cache[186].append("VS")
    known_results_cache[186].append("SV")
    print("z")
    print(known_results_cache)
    print("z")
    masses,success = recursive_find_peptides(mass, known_results_cache, 1)
 
    return masses

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: one param <int_mass>.") 

    else:
        sys.setrecursionlimit(12000)
        start = time.process_time()

        mass = int(sys.argv[1])

        int_mass_table = []
        with open("./integer_mass_table.txt") as f:
            int_mass_table = [line.rstrip() for line in f]
        #print(int_mass_table)

        init_int_mass_tables(int_mass_table)
        masses = count_peptides_of_given_mass(mass)
        #for mass in masses:
        #    print(mass)
        #print(len(masses))
        #print(len(masses) + 1) # to account for 0
        print(masses)
        print(masses + 1) # to account for 0

        end = time.process_time()
        print("Time: {0}".format(end-start))
