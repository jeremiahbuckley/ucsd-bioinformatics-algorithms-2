#! /bin/python3

import sys
import time

def init_int_mass_table(mass_table):
    int_mass_table = {}
    for kvp in mass_table:
        key = kvp[0:kvp.index(" ")]
        value = int(kvp[kvp.index(" ")+1:len(kvp)])
        int_mass_table[key] = value
    #print(int_mass_table)
    return int_mass_table

'''
idx: 0 1 2 3 4
mas: 1 2 3 4 5

cmi: 0 1 2 3  4  5
cmv: 0 1 3 6 10 15 
'''

def theoretical_spectrum(peptide_string, int_mass_table):
    mass_table = init_int_mass_table(int_mass_table)
    masses = []
    cummulative_mass_peptide_string = peptide_string + peptide_string #double so iterating through the loops is easy
    cummulative_mass_list = []
    cummulative_mass = 0
    cummulative_mass_list.append(cummulative_mass)
    for char in cummulative_mass_peptide_string:
        cummulative_mass += mass_table[char]
        cummulative_mass_list.append(cummulative_mass)

    masses.append(0)
    for length in range(1, len(peptide_string)):
        for idx in range(1, len(peptide_string)+1):
             new_mass = cummulative_mass_list[idx+length-1] - cummulative_mass_list[idx-1]
             masses.append(new_mass)
    masses.append(cummulative_mass_list[len(peptide_string)])

    return masses
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: one param, filename, format:<str_peptide>.") 

    else:
        start = time.process_time()

        peptide_string = ""
        with open(sys.argv[1]) as f:
            peptide_string = f.readline()
        peptide_string = peptide_string[0:len(peptide_string)-1]


        int_mass_table = []
        with open("./integer_mass_table.txt") as f:
            int_mass_table = [line.rstrip() for line in f]
        #print(int_mass_table)

        masses = theoretical_spectrum(peptide_string, int_mass_table)
        masses.sort()
        for mass in masses:
            print(str(mass))


        end = time.process_time()
        #print("Time: {0}".format(end-start))
