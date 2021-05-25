#! /bin/python3

import sys


if __name__ == '__main__':
    with open("./a.out") as fl:
        r_str = fl.readline()

    with open("./spectral_convolution_2_expected_result.txt") as fl:
        e_str = fl.readline()

    r_vals_1 = [int(i) for i in r_str.split(" ")]
    e_vals_1 = [int(i) for i in e_str.split(" ")]
    r_vals = sorted(r_vals_1)
    e_vals = sorted(e_vals_1)

    print("{0} {1}".format(len(r_vals), len(e_vals)))
    #assert len(r_vals) == len(e_vals)

    diff = len(r_vals) - len(e_vals)
    for i in range(len(r_vals)-1, -1+diff, -1):
        print("{0} {1}".format(r_vals[i], e_vals[i-diff]))

    for i in range(diff-1, -1, -1):
        print("{0} x".format(r_vals[i]))

    for i in range(len(r_vals)):
        assert r_vals[i] == e_vals[i]

    
