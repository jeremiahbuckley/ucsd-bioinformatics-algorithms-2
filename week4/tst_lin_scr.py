#! /bin/python3

import sys
import leaderboard_cyclopeptide_sequencing as lcs

def test_spec(pep_str, spectrum):
    pi = lcs.PeptideInfo(pep_str)
    #pi.printout()
    pi.score_against_experimental_spectrum(spectrum)
    pi.printout()


def test_pep_against_known_score(pep_str, spectrum, score):
    pi = lcs.PeptideInfo(pep_str)
    pi.score_against_experimental_spectrum(spectrum)
    print("{0} {1}".format(str(pi.linear_score), pi.linear_spec))
    assert pi.linear_score == score
    print("success")

if __name__ == "__main__":

    
    with open(sys.argv[1]) as f:
        cutoff = int(f.readline())
        spectrum_line = f.readline()

    spectrum = [int(i) for i in spectrum_line.split(" ")]

    with open("./integer_mass_table.txt") as f:
        integer_mass_table = [line.rstrip() for line in f]

    lcs.init_int_mass_tables(integer_mass_table)

    test_pep_against_known_score("PW", spectrum, 4)

    test_spec("PADFNQYVQLF", spectrum)
    test_spec("PFLQVYAGGGFDA", spectrum)



#    test_spec([97,71,115,147,114,128,163,99,128,113,147])
#    test_spec([97,147,113,128,99,163,71,57,57,57,147,115,71])
