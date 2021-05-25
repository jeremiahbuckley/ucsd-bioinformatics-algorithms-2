#! /bin/python3

import sys
import time
import itertools

def spectral_convolution_likely_weights(spectrum):
    cvs = []

    for i in range(len(spectrum)-1):
        for j in range(i+1,len(spectrum)):
            w = spectrum[j] - spectrum[i]
            if w != 0:
                cvs.append(w)

    return cvs


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: one param, filename, format:<list_int_weights>.") 

    else:
        start = time.process_time()

        with open(sys.argv[1]) as f:
            spectrum_line = f.readline()

        spectrum = [int(i) for i in spectrum_line.split(" ")]


        weights = spectral_convolution_likely_weights(sorted(spectrum))
        print(" ".join([str(w) for w in sorted(weights)]))

        end = time.process_time()
        print("Time: {0}".format(end-start))

