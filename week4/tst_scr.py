#! /bin/python3

import sys
import convolution_cyclopeptide_sequencing as cps
#calc_theoretic_spectrum(spectrum, circular=False)

def generate_test_spectrum(val, increment):
    expected_answer = [0]

    #print("val = {0}".format(str(val)))
    # intentionally do not include max-length
    for length in range(1, val):
        #print("length = {0}".format(str(length)))
        for idx in range(1, val-length+2):
            #print("idx = {0}".format(str(idx)))
            v = 0
            for s in range(idx, idx+length):
                v += s*increment
            #print("new_val = {0}".format(str(v)))
            expected_answer.append(v)
        print("   e_a {0}".format("-".join([str(i) for i in expected_answer])))
    v = 0
    for s in range(val+1):
        v += s*increment
    print("e_a max = {0}".format(str(v)))
    expected_answer.append(v)

    print(" ".join([str(i) for i in expected_answer]))
    return sorted(expected_answer)

def generate_test_peptide_weights(val, increment):
    w = []
    for i in range(1, val+1):
        w.append(i*increment)
    return w

def test_scoring(val, increment=1):
    print("test-scoring {0}".format(str(val)))
    expected_result = generate_test_spectrum(val, increment)
    peptide_weights = generate_test_peptide_weights(val, increment)

    program_result_unsorted = cps.calc_theoretic_spectrum(peptide_weights)
    program_result = sorted(program_result_unsorted)

    #assert len(program_result) == len(expected_result)
    print("lens: {0} {1}".format(str(len(program_result)), str(len(expected_result))))
    print("pr: {0}".format("-".join([str(i) for i in program_result])))
    print("er: {0}".format("-".join([str(i) for i in expected_result])))

    for i in range(len(program_result)):
        assert program_result[i] == expected_result[i]
    print("success {0}".format(str(val)))

def test_convoluted_values():
    test_spectrum = [4,12,13,23,25,35,37,40,44,48,49]
    expected_result_list_u = [8,9,1,19,11,10,21,13,12,2,31,23,22,12,10,33,25,24,14,12,2,36,28,27,17,15,5,3,40,32,31,21,19,9,7,4,44,36,35,25,23,13,11,8,4,45,37,36,26,24,14,12,9,5,1]
    expected_result_list = sorted(expected_result_list_u)
    expected_result = {}
    for e in expected_result_list:
        if e not in expected_result:
            expected_result[e] = 0
        expected_result[e] += 1


    program_weights = cps.spectral_convolution_likely_weights(test_spectrum)

    print("pw: {0} er: {1}".format(str(len(expected_result)), str(len(program_weights))))
    assert len(expected_result) == len(program_weights)

    for kvp in program_weights.items():
        assert kvp[0] in expected_result
        assert expected_result[kvp[0]] == program_weights[kvp[0]]

    print("success")


def test_convoluted_values_2():
    test_spectrum = [0,113,114,128,129,227,242,242,257,355,356,370,371,484]
    expected_result_list = [113]
    expected_result_list.extend([114,1])
    expected_result_list.extend([128,15,14])
    expected_result_list.extend([129,16,15,1])
    expected_result_list.extend([227,114,113,99,98])
    expected_result_list.extend([242,129,128,114,113,15])
    expected_result_list.extend([242,129,128,114,113,15])
    expected_result_list.extend([257,144,143,129,128,30,15,15])
    expected_result_list.extend([355,242,241,227,226,128,113,113,98])
    expected_result_list.extend([356,243,242,228,227,129,114,114,99,1])
    expected_result_list.extend([370,257,256,242,241,143,128,128,113,15,14])
    expected_result_list.extend([371,258,257,243,242,144,129,129,114,16,15,1])
    expected_result_list.extend([484,371,370,356,355,257,242,242,227,129,128,114,113])
    expected_result = {}
    for e in expected_result_list:
        if e not in expected_result:
            expected_result[e] = 0
        expected_result[e] += 1

    program_weights = cps.spectral_convolution_likely_weights(test_spectrum)

    print("pw: {0} er: {1}".format(str(len(expected_result)), str(len(program_weights))))
    #assert len(expected_result) == len(program_weights)

    print(expected_result)
    print(program_weights)

    for kvp in program_weights.items():
        print("{0}:{1}".format(str(kvp[0]), str(kvp[1])))
        print("er: {0}".format(str(expected_result[kvp[0]])))
        assert kvp[0] in expected_result
        assert expected_result[kvp[0]] == program_weights[kvp[0]]

    print("success")
    
if __name__ == "__main__":
    test_scoring(1)
    test_scoring(3)
    test_scoring(5)
    test_scoring(10)
    test_scoring(3, 3)
    test_scoring(10, 3)

    test_convoluted_values()
    test_convoluted_values_2()
