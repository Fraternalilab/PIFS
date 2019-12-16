import csv
from pifsDraw import *

for pdb in ['A3A', 'A3H']:
    for dna in ['5keg', '5sww', '5td5', '6bux']:
        fn = pdb + '_graft' + dna
        print fn
        loadNetworkOnPyMOL('output/' + fn + '_network.txt', 'pdb/' + fn + '.pdb', \
            'output/' + fn + '.pse', cutoff=6.0)
