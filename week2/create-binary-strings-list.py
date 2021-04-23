#! /bin/python3

import sys
import time
import itertools

if __name__ == "__main__":
    start = time.process_time()
    if len(sys.argv) < 2:
        print("Need a string length (int)")

    for i in itertools.product(range(2), repeat=int(sys.argv[1])-1):
        print("".join([str(istr) for istr in i]))
    end = time.process_time()
    #print("Time: {0}".format(end-start))
