#! /bin/python3

import sys
import time

def rprint(ilevel, str):
    #print("    "*ilevel + "-" + str)
    return

def build_dna_recursive(dna, append_kmer, kmers, kmer_len, ilevel):
    ilevel+=1
    dnas = []
    success = False
    rprint(ilevel, "{0} {1} {2}".format(kmer_len, dna, append_kmer))
    rprint(ilevel, " ".join(kmers))
    rprint(ilevel, "{0} =? {1}".format(dna[len(dna)-kmer_len+1:len(dna)], append_kmer[0:kmer_len-1]))
    if dna[len(dna) - kmer_len+1:len(dna)] == append_kmer[0:kmer_len-1]:
        dna += append_kmer[kmer_len-1:kmer_len]
        if len(kmers) == 0:
            dnas.append(dna)
            rprint(ilevel, "-reached bottom")
            return True, dnas
        else:
            for i in range(len(kmers)):
                d = kmers.pop(0)
                next_level_success, dna_list = build_dna_recursive(dna, d, kmers[:], kmer_len, ilevel)
                if next_level_success:
                    dnas.extend(dna_list)
                    rprint(ilevel, "--match win {0}".format(d))
                    success = True
                else:
                    kmers.append(d) # put it back on the list so other iterations can use it
            rprint(ilevel, "-try return {0} {1}".format(success, dnas))
            return success, dnas
    else:
        rprint(ilevel, "nope {0}".format(dnas))
        return False, dnas

def reconstruct_dna_from_kmer_list(kmers):
    dna = kmers[0]
    for i in range(1, len(kmers)):
        dna += kmers[i][len(kmers[i])-1:len(kmers[i])]
    return dna

def reconstruct_dna_from_kmer_list_orig(kmers):
    d = kmers.pop(0)
    dna_list = []
    print(*kmers, sep=",")
    for i in range(len(kmers)):
        append_dna = kmers.pop(0)
        success, dnas = build_dna_recursive(d, append_dna, kmers[:], len(kmers[0]), 0)
        if success:
            dna_list.extend(dnas)
        else:
            kmers.append(append_dna) # put it back on the list so other iterations can use it

    print("end..")
    print(*dna_list,sep="! ")
    ret_val = ""
    if len(dna_list) > 1:
        ret_val =  dna_list[0]
    return ret_val

if __name__ == "__main__":
    start = time.process_time()

    if len(sys.argv) < 2:
        print("Usage: need one argument, file path. File format <list_str_kmers>.")
    else:
        kmer_list = []
        with open(sys.argv[1]) as fl:
            for line in fl:
                kmer_list.append(line.rstrip())
            # can also do kmer_list = [line.rstrip() for line in fl]

        full_dna = reconstruct_dna_from_kmer_list(kmer_list)

        print(full_dna)

    end = time.process_time()
    print("Time: {0}".format(end-start))
