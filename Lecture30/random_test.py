import math
import random
import time

poor_seed = 123                         # maintain seed for poor generator
def poor():                             # example from Numerical Recipes
    global poor_seed
    seed = poor_seed
    IM = 6075 ; IA = 106 ; IC = 1283
    seed = (IA * seed + IC) % IM        # linear congruential algorithm
    poor_seed = seed
    return seed / float(IM)

park_miller_seed = 1                    # maintain seed for Park-Miller
def park_miller():                      # Park-Miller generator
    global park_miller_seed
    seed = park_miller_seed
    IA = 16807 ; IC = 2147483647 ; IQ = 127773 ; IR = 2836
    h = seed // IQ
    l = seed % IQ;
    seed = IA * l - IR * h
    if seed <= 0:
        seed += IC
    park_miller_seed = seed
    return seed / float(IC)

def test(generator, name, remark):
    print " " + name + "() " + remark

    print " checking period ..."
    start = generator()
    max_steps = 1000000
    for i in range(max_steps):
        if generator() == start:
            print " repeats after " + repr(i) + " steps"
            break
    if i+1 >= max_steps:
        print " period larger than " + repr(max_steps)

    start_time = time.clock()
    bins = 10000
    print " binning " + repr(max_steps) + " tries in " + repr(bins) + " bins"
    bin = [ 0.0 ] * bins
    for i in range(max_steps):
        b = int(generator() * bins)
        if b >= 0 and b < bins:
            bin[b] += 1
    chisqr = 0.0
    expect = max_steps / float(bins)
    for b in range(bins):
        diff = bin[b] - expect
        chisqr += diff**2
    chisqr /= expect
    cpu_time = time.clock() - start_time
    print " chi-squar/d.o.f. = " + repr(chisqr/(bins-1.0))
    print " CPU time = " + repr(cpu_time) + " seconds"

    file_name = name + ".data"
    print " writing eyeball test data to " + file_name
    file = open(file_name, "w")
    for i in range(10000):
        file.write(repr(generator()) + "\t" + repr(generator()) + "\n")
    file.close()
    print

print " Testing random number generators()"
print " ----------------------------------"
test(poor, "poor", "from Numerical Recipes period 6,075")
test(park_miller, "park_miller", "Park-Miller multiplicative generator")
test(random.random, "random.random", "from Python standard library")
