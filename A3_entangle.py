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
from pifsEntanglement import *

# A3 loop definitions 
A3_loops = 'A3_loopDefinition'
A3s = dict()
with open(A3_loops, 'r') as instream:
    read = csv.reader(instream, delimiter="\t")
    for row in read:
        if row[0] != 'domain':
            loop1 = row[1].split('-')
            loop1 = [int(loop1[0]), int(loop1[1])]
            loop3 = row[2].split('-')
            loop3 = [int(loop3[0]), int(loop3[1])]
            A3s[row[0]] = {'loop1': loop1, 'loop3': loop3}

# PDB definitions
pdb_dna = 'A3_DNA_Definition'
pdbs = dict()
with open(pdb_dna, 'r') as instream:
    read = csv.reader(instream, delimiter="\t")
    for row in read:
        if row[0] != 'pdb':
            pdbs[row[0]] = [int(row[1]), int(row[2])]


domain = sys.argv[1]
definitions = A3s[domain]

print domain
os.chdir(domain + '/')
for pdb in pdbs.keys():
    print pdb
    tarfile = 'pdb/' + domain + '_graft' + pdb + '.tar.gz'
    chainfile = 'out/' + domain + '_graft' + pdb + '_chainEntangle.txt'
    basefile = 'out/' + domain + '_graft' + pdb + '_baseEntangle.txt'
    print tarfile
    subprocess.call('tar -xzvf ' + tarfile + ' -C ' +  domain + '_graft' + pdb , shell=True)
    pdbfiles = glob.glob(domain + '_graft' + pdb + '/' + domain + '*.pdb')
    for pdbfile in pdbfiles:
        print pdbfile + ' ...'
        name_struct = pdbfile.split('/')[1]
        chain_entangle(pdbfile, 'loop1_' + name_struct, A3s[domain]['loop1'], pdbs[pdb], 'A', 'B', chainfile)
        chain_entangle(pdbfile, 'loop3_' + name_struct, A3s[domain]['loop3'], pdbs[pdb], 'A', 'B', chainfile)
        base_entangle(pdbfile, 'loop1_' + name_struct, A3s[domain]['loop1'], pdbs[pdb], 'A', 'B', basefile)
        base_entangle(pdbfile, 'loop3_' + name_struct, A3s[domain]['loop3'], pdbs[pdb], 'A', 'B', basefile)
    subprocess.call('rm -rf ' + domain + '_graft' + pdb + '/', shell=True)
