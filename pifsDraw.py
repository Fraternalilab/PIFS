'''
draw C-alpha network in PyMOL
'''
import csv
import sys
sys.path.append('/home/josephn/pymol-open-source-v2.1.0/lib/python/')
import os
from pymol import cmd

def drawCaNetwork(interactions, pdb, outfile, cutoff=6.0):
    '''
    draw C-alpha network
    '''
    # read the structure
    cmd.load(pdb, 'structure')
    # generate a list of each Calpha and select each Calpha
    C_alpha = set()
    for pair in interactions:
        C_alpha.add(pair[0])
        C_alpha.add(pair[1])
    C_alpha = list(C_alpha)    
    for res in C_alpha:
        # get chain, resn and resi
        resn, chain, resi = readResInfo(res)
        cmd.select('ca_sele', '/structure//' + chain + '/' + resi + '/' + 'CA', merge=1)
    cmd.create('ca', 'ca_sele')
    cmd.create('struct', 'structure')
    cmd.delete('structure')
    #cmd.delete('ca_sele')
    # configure the appearance
    cmd.show_as('spheres', 'ca')
    cmd.alter('/ca////', 'vdw=0.2')
    cmd.bg_color(color="white")
    cmd.color('deepsalmon', '/ca////')
    cmd.set('label_color', 'black')
    cmd.set('label_font_id', 11)
    cmd.set('label_position', (1,1,1))
    cmd.label('/ca////', 'resn+resi' )
    # generate the network
    for pair in interactions:
        atom1 = readResInfo(pair[0])
        atom2 = readResInfo(pair[1])
        cmd.distance('dist', '/ca//'+ atom1[1] + '/' + atom1[2], \
            '/ca//' + atom2[1] + '/' + atom2[2], cutoff=cutoff)
    cmd.color('blue', 'dist')
    cmd.save(outfile)
    cmd.delete('ca')
    cmd.delete('ca_sele')
    cmd.delete('dist')
    cmd.delete('structure')
    return None

def readResInfo(res):
    '''
    break down res strings ('TYR-A:101' --> ['TYR', 'A', '101'])
    '''
    res_split = res.split(':')
    resi=res_split[1]
    res_split = res_split[0].split('-')
    chain=res_split[1]
    resn=res_split[0]
    return [resn, chain, resi]

def replaceProteinChainID(res, protein_chain = None):
    '''
    replace protein chain ID in 'res' by that indicated in 'protein_chain'
    '''
    if protein_chain is None:
        return res
    res = readResInfo(res)
    return ':'.join(['-'.join([res[0], protein_chain]), res[2]])

def loadNetworkOnPyMOL(network_file, pdb, outfile, cutoff=6.0, protein_chain = None):
    '''
    read the network file and generate a PyMOL session
    @param protein_chain: if provided override what is indicated in network_file.
    '''
    interactions = []
    with open(network_file, 'r') as instream:
        read = csv.reader(instream, delimiter="\t")
        for row in read:
            if row[0] != 'res_a':
                atom1 = replaceProteinChainID(row[1], protein_chain)
                atom2 = replaceProteinChainID(row[3], protein_chain)
                if ([atom1, atom2] not in interactions) and ([atom2, atom1] not in interactions):
                    interactions.append([atom1, atom2])
    drawCaNetwork(interactions, pdb, outfile, cutoff)
    return None
