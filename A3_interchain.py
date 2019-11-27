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
    tarfile_restrained = domain + '_graft' + pdb + '_restrained.tar.gz'
    comparisonFile = domain + '_graft' + pdb + '.summary'
    comparisonFile_restrained = domain + '_graft' + pdb + '_restrained.summary'
    # first the non-restrained ones
    print tarfile
    subprocess.call('tar -xzvf ' + tarfile, shell=True)
    pdbfiles = glob.glob(domain + '_graft' + pdb + '_dna_*.pdb')
    orig = domain + '_graft' + pdb + '.pdb'
    find_contacts(orig, orig[:-4], pdbs[pdb], 'A', 'B', distCutoff=6.0)
    for pdbfile in pdbfiles:
        print pdbfile + ' ...'
        find_contacts(pdbfile, pdbfile[:-4], pdbs[pdb], 'A', 'B', distCutoff=6.0)
        # compareHbonds(orig[:-4] + '.hbonds', pdbfile[:-4] + '.hbonds', comparisonFile)
        # comparePiStack(orig[:-4] + '.pistack', pdbfile[:-4] + '.pistack', comparisonFile)
    subprocess.call('rm ' + domain + '_graft' + pdb + '*.pdb', shell=True)
    # then the restrained ones
    print tarfile_restrained
    subprocess.call('tar -xzvf ' + tarfile_restrained, shell=True)
    pdbfiles = glob.glob(domain + '_graft' + pdb + '_dna_.pdb')
    orig = domain + '_graft' + pdb + '.pdb'
    find_contacts(orig, orig[:-4], pdbs[pdb], 'A', 'B', distCutoff=6.0)
    for pdbfile in pdbfiles:
        print pdbfile + ' ...'
        find_contacts(pdbfile, pdbfile[:-4], pdbs[pdb], 'A', 'B', distCutoff=6.0)
        # compareHbonds(orig[:-4] + '.hbonds', pdbfile[:-4] + '.hbonds', comparisonFile_restrained)
        # comparePiStack(orig[:-4] + '.pistack', pdbfile[:-4] + '.pistack', comparisonFile_restrained)
    subprocess.call('rm ' + domain + '_graft' + pdb + '*.pdb', shell=True)
