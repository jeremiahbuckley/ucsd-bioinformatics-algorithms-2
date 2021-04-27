#! /bin/python3

import utils
import functools
import time

def test_reverse_complement():
    out = utils.reverse_complement("ACGT")
    assert out == "ACGT"
    out = utils.reverse_complement("")
    assert len(out) == 0
    out = utils.reverse_complement("AAAGGG")
    assert out == "CCCTTT"

def test_create_profile_from_motifs():
    motifs = ["ACCGG", "ACCTA", "TCCGG", "AGCTT"]
    profile=utils.create_profile_from_motifs(motifs)
    #print(profile)
    assert len(profile) == 4
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[0],[0.75, 0.0, 0.0, 0.0, 0.25]),
                            True)
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[1],[0.0, 0.75, 1.0, 0.0, 0.0]),
                            True)
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[2],[0.0, 0.25, 0.0, 0.5, 0.5]),
                            True)
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[3],[0.25, 0.0, 0.0, 0.5, 0.25]),
                            True)
    profile=utils.create_profile_from_motifs(motifs, 1)
    #print(profile)
    assert len(profile) == 4
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[0],[0.5,0.125,0.125,0.125,0.25]),
                            True)
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[1],[0.125,0.5,0.625,0.125,0.125]),
                            True)
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[2],[0.125,0.25,0.125,0.375,0.375]),
                            True)
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[3],[0.25,0.125,0.125,0.375,0.25]),
                            True)

    motifs = ["TCGGGGGTTTTT","CCGGTGACTTAC","ACGGGGATTTTC","TTGGGGACTTTT","AAGGGGACTTCC","TTGGGGACTTCC","TCGGGGATTCAT","TCGGGGATTCCT","TAGGGGAACTAC","TCGGGTATAACC"]
    profile=utils.create_profile_from_motifs(motifs)
    #print(profile)
    assert len(profile) == 4
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[0],[0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0]),
                            True)
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[1],[0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6]),
                            True)
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[2],[0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0]),
                            True)
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[3],[0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]),
                            True)
    profile=utils.create_profile_from_motifs(motifs,1)
    #print(profile)
    assert len(profile) == 4
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[0],[0.21428571428571427, 0.21428571428571427, 0.07142857142857142, 0.07142857142857142, 0.07142857142857142, 0.07142857142857142, 0.7142857142857143, 0.14285714285714285, 0.14285714285714285, 0.14285714285714285, 0.2857142857142857, 0.07142857142857142]),
                            True)
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[1],[0.14285714285714285, 0.5, 0.07142857142857142, 0.07142857142857142, 0.07142857142857142, 0.07142857142857142, 0.07142857142857142, 0.35714285714285715, 0.14285714285714285, 0.21428571428571427, 0.35714285714285715, 0.5]),
                            True)
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[2],[0.07142857142857142, 0.07142857142857142, 0.7857142857142857, 0.7857142857142857, 0.7142857142857143, 0.7142857142857143, 0.14285714285714285, 0.07142857142857142, 0.07142857142857142, 0.07142857142857142, 0.07142857142857142, 0.07142857142857142]),
                            True)
    assert functools.reduce(lambda n1, n2: n1 and n2, 
                            map(lambda out,test: out==test, profile[3],[0.5714285714285714, 0.21428571428571427, 0.07142857142857142, 0.07142857142857142, 0.14285714285714285, 0.14285714285714285, 0.07142857142857142, 0.42857142857142855, 0.6428571428571429, 0.5714285714285714, 0.2857142857142857, 0.35714285714285715]),
                            True)
def test_calc_basic_score_from_motifs_list():
    motifs = ["ACCGG", "ACCTA", "TCCGG", "AGCTT"]
    score=utils.calc_basic_score_from_motifs_list(motifs)
    #print(score)
    assert score == 6
    score=utils.calc_basic_score_from_motifs_list(motifs, 1)
    #print(score)
    assert score == 1

    motifs = ["TCGGGGGTTTTT","CCGGTGACTTAC","ACGGGGATTTTC","TTGGGGACTTTT","AAGGGGACTTCC","TTGGGGACTTCC","TCGGGGATTCAT","TCGGGGATTCCT","TAGGGGAACTAC","TCGGGTATAACC"]
    score=utils.calc_basic_score_from_motifs_list(motifs)
    #print(score)
    assert score == 30
    score=utils.calc_basic_score_from_motifs_list(motifs, 1)
    #print(score)
    assert score == 18


if __name__ == "__main__":
    start = time.process_time()

    test_reverse_complement()

    test_create_profile_from_motifs()
    
    test_calc_basic_score_from_motifs_list()

    end = time.process_time()
    print("Time: {0}".format(end-start))
