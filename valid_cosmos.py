#!/usr/bin/env python2.7

# Written by Dominique Boutigny with modifications by Michael Wood-Vasey

import os.path
import sys

import matplotlib.pylab as plt
import numpy as np

from check_astrometry import main, loadAndMatchData

def defaultData(repo):
    # List of visits to be considered
    visits = [176837, 176846, 176850]

    # Reference visit (the other visits will be compared to this one)
    ref = 176837

    # List of CCD to be considered (source catalogs will be concateneted)
    ccd = [10] #, 12, 14, 18]
    filter = 'z'
    
    # Reference values for the median astrometric scatter and the number of matches
    good_mag_limit = 21
    medianRef = 25
    matchRef = 5600

    return visits, ref, ccd, filter, good_mag_limit, medianRef, matchRef

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("""Usage: valid_cosmos repo
where repo is the path to a repository containing the output of processCcd
""")
        sys.exit(1)

    repo = sys.argv[1]
    if not os.path.isdir(repo):
        print("Could not find repo %r" % (repo,))
        sys.exit(1)

    visits, ref, ccd, filter, good_mag_limit, medianRef, matchRef = defaultData(repo)
    main(repo, visits, ref, ccd, filter, good_mag_limit, medianRef, matchRef)
