#! /bin/python3

import sys
import peptides_with_given_mass as pep 

import time
import os

def mfunction(garbage):
    return 10

def test_count_peptides_of_given_mass():
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, "count_per_mass.txt")
    with open(filename) as f:
        lines = f.readlines()

    int_mass_table = []
    with open("./integer_mass_table.txt") as f:
        int_mass_table = [line.rstrip() for line in f]

    test_values = []
    for l in lines:
        in_val = int(l[0:l.index(" ")])
        out_val = int(l[l.index(" ")+4:len(l)])
        test_values.append([in_val, out_val])

    pep.init_int_mass_tables(int_mass_table)
    i = 0
    for t in test_values:
        result = pep.count_peptides_of_given_mass(t[0])
        print("{0} {1} {2}".format(t[0], t[1], result))
        assert t[1] == result
        i += 1
    print("{0} tests passed.".format(i))

if __name__ == "__main__":
    start = time.process_time()

    masses = test_count_peptides_of_given_mass()

    end = time.process_time()
