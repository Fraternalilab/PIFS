'''
heuristic to run through all A3-DNA graft structures and detect entanglement:
- DNA entangled in loop1
- DNA entangled in loop3
considering both DNA chain and individual DNA base.
'''
import csv
import glob
import os
import subprocess
import sys
from pifsCalcInteractions import *

# PDB definitions
pdb_dna = 'A3_DNA_Definition'
pdbs = dict()
with open(pdb_dna, 'r') as instream:
    read = csv.reader(instream, delimiter="\t")
    for row in read:
        if row[0] != 'pdb':
            pdbs[row[0]] = [int(row[1]), int(row[2])]

domain = sys.argv[1]
print domain
for pdb in pdbs.keys():
    print pdb
    tarfile = domain + '_graft' + pdb + '.tar.gz'
    comparisonFile = domain + '_graft' + pdb + '.summary'
    print tarfile
    subprocess.call('mkdir ' + domain + '_graft' + pdb, shell=True )
    subprocess.call('tar -xzvf ' + tarfile + ' -C ' + domain + '_graft' + pdb , shell=True)
    pdbfiles = glob.glob(domain + '_graft' + pdb + '/' + domain + '_graft' + pdb + '_dna_*.pdb')
    orig = domain + '_graft' + pdb + '/' + domain + '_graft' + pdb + '.pdb'
    name_struct = domain + '_graft' + pdb 
    find_contacts(orig, name_struct, pdbs[pdb], 'A', 'B', distCutoff=6.0)
    for pdbfile in pdbfiles:
        print pdbfile + ' ...'
        name_struct = pdbfile.split('/')[1][:-4]
        find_contacts(pdbfile, name_struct, pdbs[pdb], 'A', 'B', distCutoff=6.0)
    subprocess.call('rm -rf ' + domain + '_graft' + pdb + '/', shell=True)
