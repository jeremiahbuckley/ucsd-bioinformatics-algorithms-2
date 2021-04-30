#! /bin/python3

import sys
import time
import array

_peptide_by_mass_lookup_ = {}

_peptides_at_mass_idx_ = array.array('l', [0 for i in range(5000)])
_cache_p_res_at_m_idx_ = array.array('i', [0 for i in range(5000)])
_peptide_masses_ = []

def init_int_mass_tables(mass_table):
    mass_by_peptide = {}
    for kvp in mass_table:
        key = kvp[0:kvp.index(" ")]
        value = int(kvp[kvp.index(" ")+1:len(kvp)])
        mass_by_peptide[key] = value
        if key != "I" and key != "K":
            _peptide_by_mass_lookup_[value] = key
            _peptides_at_mass_idx_[value] = 1
            _cache_p_res_at_m_idx_[value] = 1
            _peptide_masses_.append(value)
    print(mass_by_peptide)
    print(_peptide_by_mass_lookup_)
    return mass_by_peptide, _peptide_by_mass_lookup_

'''
idx: 0 1 2 3 4
mas: 1 2 3 4 5

cmi: 0 1 2 3  4  5
cmv: 0 1 3 6 10 15 
'''

def recursive_find_peptides(mass, known_results_cache, debug_indent):
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
    known_results_cache[113] += 1 #I
    known_results_cache[128] += 1 #K
    known_results_cache[114] += 2 #GG,GG(reversed)
    known_results_cache[128] += 2 #GA,AG
    known_results_cache[156] += 2 #GS,SG
    known_results_cache[186] += 2 #AD,DA
    known_results_cache[186] += 2 #EG,GE
    print("z")
    print(known_results_cache)
    print("z")
    counts = recursive_find_peptides(mass, known_results_cache, [])
 
    return counts

def count_recursive_find_peptides_old (mass, known_results_cache, debug_indent):
    #_peptides_at_mass_idx_ = array.array('l', [0 for i in range(5000)])
    #_cache_p_res_at_m_idx_ = array.array('i', [0 for i in range(5000)])
    masses = []
    #print("{0}--> {1}.".format(" " * debug_indent, str(mass)))

    if _cache_p_res_at_m_idx_[mass] == 1:
        found_results = _peptides_at_mass_idx_[mass]
        #print("{0}cache: {1} {2}.".format(" " * debug_indent, str(mass), ",".join(known_results_cache[mass])))
        return found_results

    successful_finds = 0   
    for n_mass in _peptide_masses_:
        found_results = count_recursive_find_peptides(n_mass, known_results_cache, debug_indent+1)
        successful_finds += found_results

    if _cache_p_res_at_m_idx_[mass] == 1:
        raise RuntimeError("Cache miss. Cache should have been empty here. {0}.".format(mass))
    _peptides_at_mass_idx_[mass] = successful_finds
    _cache_p_res_at_m_idx_[mass] = 1
    return successful_finds

def count_recursive_find_peptides(mass, known_results_cache, debug_indent):
    #_peptides_at_mass_idx_ = array.array('l', [0 for i in range(5000)])
    #_cache_p_res_at_m_idx_ = array.array('i', [0 for i in range(5000)])
    masses = []
    #print("{0}--> {1}.".format(" " * debug_indent, str(mass)))

    if _cache_p_res_at_m_idx_[mass] == 1:
        found_results = _peptides_at_mass_idx_[mass]
        #print("{0}cache: {1} {2}.".format(" " * debug_indent, str(mass), ",".join(known_results_cache[mass])))
        return found_results

    successful_finds = 0   
    for n_mass in _peptide_masses_:
        found_results = count_recursive_find_peptides(n_mass, known_results_cache, debug_indent+1)
        successful_finds += found_results

    if _cache_p_res_at_m_idx_[mass] == 1:
        raise RuntimeError("Cache miss. Cache should have been empty here. {0}.".format(mass))
    _peptides_at_mass_idx_[mass] = successful_finds
    _cache_p_res_at_m_idx_[mass] = 1




    recursion_stack.append([0, mass, -1])
    recursion_stack.append([0, mass, -1])
    while len(recursion_stack) > 1:
        mass = recursion_stack[len(recursion_stack)-1][1]

        if mass in known_results_cache:
            found_results = known_results_cache[mass]
            recursion_stack[len(recursion_stack)-2][2] += found_results
            recursion_stack.pop(len(recursion_stack)-1]

        succesful_finds = 0
        for kvp in _peptide_by_mass_lookup_.items():
            n_mass = mass-kvp[0]
            recursion_stack.append([kvp[0],n_mass,-1])

        
            

    return successful_finds

def count_peptides_of_given_mass(mass, int_mass_table):
    #_peptides_at_mass_idx_ = array.array('l', [0 for i in range(5000)])
    #_cache_p_res_at_m_idx_ = array.array('i', [0 for i in range(5000)])
    mass_by_peptide = init_int_mass_tables(int_mass_table)
    counts = 0
    #    print("{0} {1} {2}".format(kvp[1], kvp[0]*2, kvp[0]*3))
    #for kvp in _peptide_by_mass_lookup_.items():
    #    for kvp2 in _peptide_by_mass_lookup_.items():
    #        print("{0} {1}".format(kvp[0] + kvp2[0], kvp[1] + kvp2[1]))
    _peptides_at_mass_idx_[113] += 1 #I
    _peptides_at_mass_idx_[128] += 1 #K
    _peptides_at_mass_idx_[114] += 2 #GG, GG (reversed)
    _peptides_at_mass_idx_[128] += 2 #GA,AG
    _peptides_at_mass_idx_[156] += 2 #GS,SG
    _peptides_at_mass_idx_[186] += 2 #AD,DA
    _peptides_at_mass_idx_[186] += 2 #EG,GE
    print("z")
    print(_peptides_at_mass_idx_)
    print("z")
    counts = recursive_find_peptides(mass, known_results_cache, 1)
 
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
    #        print("{0} {1}".format(kvp[0] + kvp2[0], kvp[1] + kvp2[1]))
    #known_results_cache[113].append("I")
    #kinown_results_cache[128].append("K")
    known_results_cache[114].append("GG")
    known_results_cache[128].append("GA")
    known_results_cache[156].append("GS")
    known_results_cache[186].append("AD")
    known_results_cache[186].append("EG")
    print("z")
    print(known_results_cache)
    print("z")
    masses,success = recursive_find_peptides(mass, known_results_cache, 1)
 
    return masses

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

        masses = peptides_of_given_mass(mass, int_mass_table)
        #for mass in masses:
        #    print(mass)
        print(len(masses))


        end = time.process_time()
        print("Time: {0}".format(end-start))
