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
    pi.printout()
    assert pi.linear_score == score
    print("success")

def test_full_leaderboard(spectrum, cutoff, integer_mass_table):
    expected_values = ['71-115-147-114-128-163-99-128-113-147-97', \
'97-71-115-147-114-128-163-99-128-113-147', \
'97-147-113-128-99-163-128-114-147-71-115', \
'97-147-113-128-99-163-128-114-147-115-71', \
'99-163-128-114-147-71-115-97-147-113-128', \
'99-163-128-114-147-115-71-97-147-113-128', \
'113-128-99-163-128-114-147-71-115-97-147', \
'113-128-99-163-128-114-147-115-71-97-147', \
'114-128-163-99-128-113-147-97-71-115-147', \
'114-128-163-99-128-113-147-97-115-71-147', \
'115-71-147-114-128-163-99-128-113-147-97', \
'115-147-114-128-163-99-128-113-147-97-71', \
'128-99-163-128-114-147-71-115-97-147-113', \
'128-99-163-128-114-147-115-71-97-147-113', \
'128-114-147-115-71-97-147-113-128-99-163', \
'147-113-128-99-163-128-114-147-71-115-97', \
'147-113-128-99-163-128-114-147-115-71-97', \
'147-114-128-163-99-128-113-147-97-71-115', \
'147-114-128-163-99-128-113-147-97-115-71', \
'163-128-114-147-115-71-97-147-113-128-99', \
'57-57-128-163-99-128-113-147-97-71-115-147', \
'97-147-113-128-99-163-57-71-114-147-71-115', \
'97-147-113-128-99-163-57-71-114-147-115-71', \
'97-147-113-128-99-163-71-57-114-147-71-115', \
'97-147-113-128-99-163-71-57-114-147-115-71', \
'97-147-113-128-99-163-128-57-57-147-71-115', \
'97-147-113-128-99-163-128-57-57-147-115-71', \
'99-163-128-114-147-71-115-97-147-113-57-71', \
'99-163-128-114-147-71-115-97-147-113-71-57', \
'99-163-128-114-147-115-71-97-147-113-57-71', \
'99-163-128-114-147-115-71-97-147-113-71-57', \
'147-113-128-99-163-128-57-57-147-115-71-97', \
'147-114-128-163-99-57-71-113-147-97-71-115', \
'147-114-128-163-99-57-71-113-147-97-115-71', \
'147-114-128-163-99-71-57-113-147-97-71-115', \
'147-114-128-163-99-71-57-113-147-97-115-71', \
'97-147-113-128-99-163-71-57-57-57-147-71-115', \
'97-147-113-128-99-163-71-57-57-57-147-115-71']

    sequences = lcs.leaderboard_cyclopeptide_sequencing(spectrum, cutoff, integer_mass_table)
    seq_as_weights = lcs.convert_amino_acid_strs_to_weights(sequences)

    assert len(seq_as_weights) == len(expected_values)

    for ev in expected_values:
        assert ev in seq_as_weights

    print("success - full leaderboard test")

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

    test_full_leaderboard(spectrum, cutoff, integer_mass_table)

#    test_spec([97,71,115,147,114,128,163,99,128,113,147])
#    test_spec([97,147,113,128,99,163,71,57,57,57,147,115,71])
