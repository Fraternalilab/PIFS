'''
compute intramolecular interactions (h-hond, pi-pi stacking in protein-nucleic acid grafts
'''
from Bio.PDB import *
import csv
import itertools
import sys
import glob
import os
from sympy import Point3D, Line3D, Plane
import math
import numpy as np
from shapely.geometry import Point, Polygon

class PDBcoord:
    '''
    deals with reading and extracting atom coords
    '''
    def __init__(self):
        self.pdb = None

    def readPDB(self, pdb, name):
        structParser = PDBParser()
        self.pdb = structParser.get_structure(name, pdb)

    def getResName(self, res, chain):
        models = self.pdb.get_list()
        if len(models) > 0:
            pdb1 = self.pdb[0] # just take the first model for simplicity
            if pdb1.has_id(chain):
                pdbchain = pdb1[chain]
                if pdbchain.has_id(int(res)):
                    # get res in pdbchain
                    resEntry = pdbchain[int(res)]
                    return resEntry.get_resname().strip()
                else:
                    return "getAtomsCoord: pdb does not have residue %s" % int(res)
            else:
                return "getAtomsCoord: chain not in given pdb. Exit."
        else:
            return "getAtomsCoord: pdb file has no models in it. Exit."

    def getAtomsCoord(self, res, chain, atomType):
        '''
        Get atom coordinates as array (float, size=3)
        '''
        models = self.pdb.get_list()
        if len(models) > 0:
            pdb1 = self.pdb[0] # just take the first model for simplicity
            if pdb1.has_id(chain):
                pdbchain = pdb1[chain]
                if pdbchain.has_id(int(res)):
                    # get res in pdbchain
                    resEntry = pdbchain[int(res)]
                    # fetch C-alphas and append
                    if resEntry.has_id(atomType):
                        return resEntry[atomType].get_coord().tolist()
                    else:
                        return "getAtomsCoord: res %s does not have the atom type %s" % (int(res), atomType)
                else:
                    return "getAtomsCoord: pdb does not have residue %s" % int(res)
            else:
                return "getAtomsCoord: chain not in given pdb. Exit."
        else:
            return "getAtomsCoord: pdb file has no models in it. Exit."

    def getAtom(self, res, atomType):
        '''
        Get atom object
        '''
        if res.has_id(atomType):
            return res[atomType]
        else:
            return "getAtom: %s atom type is not present for this residue." % atomType

def calcAtomDist(n1, n2, res=False):
    '''
    calculate distance between atoms n1 and n2. if res=True then calculate distance between the 
    center-of-mass of atom groups (residues) n1 and n2, and output the maximum of this list.
    '''
    if res is False:
        return n1 - n2
    else:
        # calculate center of mass of n1 and of n2
        atom1 = [atom.get_coord().tolist() for atom in Selection.unfold_entities(n1,'A')]
        atom2 = [atom.get_coord().tolist() for atom in Selection.unfold_entities(n2,'A')]
        com1 = [ np.mean([i[0] for i in atom1]), np.mean([i[1] for i in atom1]), np.mean([i[2] for i in atom1]) ]
        com2 = [ np.mean([i[0] for i in atom2]), np.mean([i[1] for i in atom2]), np.mean([i[2] for i in atom2]) ]
        return ( (com1[0] - com2[0]) ** (2) + (com1[1] - com2[1]) ** (2) + (com1[2] - com2[2]) ** (2) ) ** (0.5)

def stacking_angle(plane1, plane2):
    '''
    test whether the nucleic acid is embedded/wrapped/entangled in a specified loop
    a position (base) of the nucleic acid molecule is specified as the vector and
    part of the loop is specified as a plane.
    Amounts to calculating the angle between the vector and the plane
    '''
    if any(type(x) is str for x in plane1) or any(type(x) is str for x in plane2):
        return None
    plane01 = Plane(Point3D(plane1[0], evaluate=False), Point3D(plane1[1], evaluate=False), \
        Point3D(plane1[2], evaluate=False))
    plane02 = Plane(Point3D(plane2[0], evaluate=False), Point3D(plane2[1], evaluate=False), \
        Point3D(plane2[2], evaluate=False))
    angle_rad = plane01.angle_between(plane02)
    return angle_rad

def classify_stacking(angle):
    '''
    classify stacking as 'face-to-face' or 'edge-to-face'
    see DOI: 10.1039/C7MD00381A 
    '''
    if angle is None:
        return None
    angle = math.degrees(angle)
    if angle <= 30.0 or angle >= 150.0:
        return 'face-to-face'
    elif 60.0 <= angle <= 120.0:
        return 'edge-to-face'
    else:
        return None

def definePlane_res(pdbfile, pdbname, res):
    '''
    return Plane object define by atoms, depending on the res type
    '''
    pdb = PDBcoord()
    pdb.readPDB(pdbfile, pdbname)
    restype = res.get_resname()

    if restype in [' DT', ' DC']:
        atomCoords = [res[atomType].get_coord().tolist() for atomType in ['N1', 'N3', 'C5']]
    elif restype in [' DA', ' DG']:
        atomCoords = [res[atomType].get_coord().tolist() for atomType in ['N9', 'C2', 'C6']]
    elif restype == 'HIS':
        atomCoords = [res[atomType].get_coord().tolist() for atomType in ['CB', 'CE1', 'NE2']]
    elif restype == 'TYR':
        atomCoords = [res[atomType].get_coord().tolist() for atomType in ['CB', 'CE1', 'CE2']]
    elif restype == 'TRP':
        atomCoords = [res[atomType].get_coord().tolist() for atomType in ['CB', 'CZ2', 'CZ3']]
    elif restype == 'PHE':
        atomCoords = [res[atomType].get_coord().tolist() for atomType in ['CB', 'CE1', 'CE2']]
    else:
        return None

    if any(type(x) is str for x in atomCoords):
        return None
    elif Point3D.are_collinear(Point3D(atomCoords[0]), Point3D(atomCoords[1]), Point3D(atomCoords[2])) is False: 
        return atomCoords

    return None

def inlist(test, ls):
    '''
    test whether any element in 'test' is in list
    '''
    o = [i for i in test if i in ls]
    if len(o) > 0:
        return o
    else:
        return None

"""
def hbond_isdonor(pdbfile, pdbname, res, chain):
    '''
    test whether a residue can be a hbond donor based on atom type it has
    '''
    pdb = PDBcoord()
    pdb.readPDB(pdbfile, pdbname)
    atoms = [atom.name for atom in Selection.unfold_entities(res, 'A')]
    atomTypes = list()
    if 'N' in atoms and res.get_resname() != 'PRO':
        atomTypes.append('N')
    elif inlist(['O', 'OE1', 'OE2', 'NE1', 'NE2', 'ND1', 'OD', 'OG', 'NZ'], atoms) is not None:
        atomTypes += inlist(['O', 'OE1', 'OE2', 'NE1', 'NE2', 'ND1', 'OD', 'OG', 'NZ'], atoms) 
    elif inlist(['N1', 'N3', 'N7'], atoms) is not None and res.get_resname() == ' DA':
        atomTypes += inlist(['N1', 'N3', 'N7'], atoms) 
    elif inlist(['N3', 'N7', 'O6'], atoms) is not None and res.get_resname() == ' DG':
        atomTypes += inlist(['N3', 'N7', 'O6'], atoms) 
    elif inlist(['N3', 'O2'], atoms) is not None and res.get_resname() == ' DC':
        atomTypes += inlist(['N3', 'O2'], atoms) 
    elif inlist(['O2', 'O4'], atoms) is not None and res.get_resname() == ' DT':
        atomTypes += inlist(['O2', 'O4'], atoms) 
    if inlist(['OP1', 'OP2', 'O1P', 'O2P', "O1'", "O2'", "O3'", "O4'", "O5'"], atoms) is not None and res.get_resname() in [' DT', ' DA', ' DC', ' DG']:
        atomTypes += inlist(['OP1', 'OP2', 'O1P', 'O2P', "O1'", "O2'", "O3'", "O4'", "O5'"], atoms) 
    if len(atomTypes) > 0:
        return {chain + ':' + res.get_resname().strip() + str(res.get_id()[1]) + ':' + atomType: pdb.getAtom(res, atomType) for atomType in atomTypes}
    else:
        return False
    return None

def hbond_isacceptor(pdbfile, pdbname, res, chain):
    '''
    test whether a residue can be a hbond acceptor based on atom type it has
    '''
    pdb = PDBcoord()
    pdb.readPDB(pdbfile, pdbname)
    atoms = [atom.name for atom in Selection.unfold_entities(res, 'A')]
    atomTypes = list()
    if inlist(['NH1', 'NH2', 'OH'], atoms) is not None:
        atomTypes += inlist(['NH1', 'NH2', 'OH'], atoms)
    elif 'N6' in atoms and res.get_resname() == ' DA':
        atomTypes.append('N6') 
    elif 'N2' in atoms and res.get_resname() == ' DG':
        atomTypes.append('N2') 
    elif 'N4' in atoms and res.get_resname() == ' DC':
        atomTypes.append('N4')
    elif 'N3' in atoms and res.get_resname() == ' DT':
        atomTypes.append('N3')
    if len(atomTypes) > 0:
        return {chain + ':' + res.get_resname().strip() + str(res.get_id()[1]) + ':' + atomType: pdb.getAtom(res, atomType) for atomType in atomTypes}
    else:
        return False
    return None

def calcHBond(donor, acceptor, distcutoff=5.5):
    '''
    calculate hydrogen bonds between donor and acceptor
    '''
    hbonds = {}
    i = 1
    for donorKey, donorAtom in donor.iteritems():
        for acceptorKey, acceptorAtom in acceptor.iteritems():
            if donorAtom.get_parent() != acceptorAtom.get_parent() and donorKey[0] != acceptorKey[0]: 
                dist = calcAtomDist(donorAtom, acceptorAtom)
                if dist <= distcutoff:
                    hbonds[str(i)] = {'donor': donorKey, 'acceptor': acceptorKey, 'dist': dist}
                    i += 1
    return hbonds
"""
def NeighbourBetweenChains(structure, name, chainToDefine, chainToSearch, residToDefine, \
    atomToDefine = 'CA', atomToSearch = "CA", distCutoff = 6.0):
    '''
    Extract neighbours of a given list of residues in "resid", based on distances below a cut-off specified in "distCutoff"
    '''
    if isinstance(distCutoff, float) is False:
        return "distCutoff must be a float. Exit."
    #filepath = structure.split('/')[:-1] + name + '.pdb'
    if glob.glob(structure):
        structParser = PDBParser()
        pdb = structParser.get_structure(name, structure)
        models = pdb.get_list()
        if len(models) > 0:
            pdb1 = pdb[0] # just take the first model for simplicity
            if pdb1.has_id(chainToSearch) and pdb1.has_id(chainToDefine):
                pdbchain = pdb1[chainToSearch]
                all_search = [atom for atom in Selection.unfold_entities(pdbchain,'A') if atom.name == atomToSearch]
                c_define = list()
                neighbors = list()
                for res in residToDefine:
                    if pdb1[chainToDefine].has_id(int(res)):
                        # get res in pdbchain
                        resEntry = pdb1[chainToDefine][int(res)]
                        # put itself into the neighbor list
                        neighbors.append(resEntry)
                        # fetch C-alphas and append
                        if resEntry.has_id(atomToDefine):
                            c = resEntry[atomToDefine]
                            c_define.append(c)
                        else:
                            print "res %s not in given chain. Skip." % int(res)
                            pass
                if len(c_define) > 0:
                    print "Fetched a total of %s %s." % (len(c_define), atomToDefine)
                    search = NeighborSearch(all_search)
                    for c in c_define:
                         neighbors = neighbors + search.search(c.get_coord(), distCutoff, level="R")  # Fetch on residue level
                         neighbors = list(set(neighbors)) # unique entries
                         # print "A total of %s residues to consider." % len(neighbors)
                else:
                    print "Error in extracting c-alphas from structure. Exit."
                return neighbors
            else:
                return "chain not in given pdb. Exit."
        else:
           return "pdb file has no models in it. Exit."
    else:
        return "File not found. Exit."

def find_contacts(pdb, name, resrangeNA, chainProt, chainNA, distCutoff=6.0):
    '''
    the pipeline to find h-bonds and pi-stacking and write output files
    '''
    neighbours = NeighbourBetweenChains(pdb, name, chainNA, chainProt, range(resrangeNA[0], resrangeNA[1] + 1), \
        atomToDefine = "C3'", atomToSearch = 'CA', distCutoff = 8.0)
    """
    # first caclulate h-bonds
    print 'calculating h-bonds ...'
    donors = {}
    acceptors = {}
    for neighbour in neighbours:
        isdonor = hbond_isdonor(pdb, name, neighbour, neighbour.get_parent().get_id())
        isacceptor = hbond_isacceptor(pdb, name, neighbour, neighbour.get_parent().get_id())
        if isdonor is not None and isdonor is not False:
            donors.update(isdonor)
        if isacceptor is not None and isacceptor is not False:
            acceptors.update(isacceptor)
    hbonds = calcHBond(donors, acceptors, distcutoff=distCutoff)
    """
    # then calculate pi-pi stacking
    print 'calculating stacking ...'
    pi_stacks = {}
    x = 1
    n_pair = itertools.combinations(neighbours, 2)
    for i, j in n_pair:
        dist = calcAtomDist(i, j, res=True)
        if dist <= distCutoff and i.get_parent() != j.get_parent():
            plane_i = definePlane_res(pdb, name, i)
            plane_j = definePlane_res(pdb, name, j)
            if plane_i is not None and plane_j is not None:
                angle = stacking_angle(plane_i, plane_j)
                classification = classify_stacking(angle)
                res_i = i.get_parent().get_id() + ':' + i.get_resname().strip() + str(i.get_id()[1])
                res_j = j.get_parent().get_id() + ':' + j.get_resname().strip() + str(j.get_id()[1])
                pi_stacks[x] = {'res_i': res_i, 'res_j': res_j, \
                    'angle': math.degrees(angle), 'class': classification}
                x += 1
    """
    with open(name + '.hbonds', 'w') as outstream:
        wr = csv.writer(outstream, delimiter="\t")
        wr.writerow(['pdb', 'donor', 'acceptor', 'distance'])
        for key, hbond in hbonds.iteritems():
            wr.writerow([name, hbond['donor'], hbond['acceptor'], str(hbond['dist'])])
    """
    with open(name + '.pistack', 'w') as outstream:
        wr = csv.writer(outstream, delimiter="\t")
        wr.writerow(['pdb', 'res_1', 'res_2', 'angle', 'classification'])
        for key, pistack in pi_stacks.iteritems():
            wr.writerow([name, pistack['res_i'], pistack['res_j'], str(pistack['angle']), str(pistack['class'])])
    
    return None
"""
def compareHbonds(hb1, hb2, compareSummaryFile):
    '''
    compare two hbonds list, output the number and details of bonds found in both hb1 and hb2,
    and those found only in hb1 but not hb2
    write to compareSummaryFile
    '''
    hb1_bonds = []
    hb2_bonds = []
    with open(hb1, 'r') as instream:
        read = csv.reader(instream, delimiter="\t")
        for row in read:
            if row[0] != 'pdb':
                hb1_bonds.append(row[1] + '/' + row[2])
    with open(hb2, 'r') as instream:
        read = csv.reader(instream, delimiter="\t")
        for row in read:
            if row[0] != 'pdb':
                hb2_bonds.append(row[1] + '/' + row[2])
    overlap = [ bond for bond in list(set(hb1_bonds + hb2_bonds)) if (bond in hb1_bonds and bond in hb2_bonds) ]
    unique = [ bond for bond in list(set(hb1_bonds + hb2_bonds)) if (bond in hb1_bonds and bond not in hb2_bonds) ]
    if os.path.isfile(compareSummaryFile) is False:
        with open(compareSummaryFile, 'a') as outstream:
            wr = csv.writer(outstream, delimiter="\t")
            wr.writerow(['pdb1', 'pdb2', 'type', 'overlap_n', 'overlap', 'unique_in_pdb1_n', 'unique_in_pdb1'])
    with open(compareSummaryFile, 'a') as outstream:
        wr = csv.writer(outstream, delimiter="\t")
        wr.writerow([hb1[:-7], hb2[:-7], 'H-bond', str(len(overlap)), ';'.join(overlap), str(len(unique)), ';'.join(unique)])

def comparePiStack(ps1, ps2, compareSummaryFile):
    '''
    compare two pistack list, output the number and details of pi-stacks found in both ps1 and ps2,
    and those found only in ps1 but not ps2
    write to compareSummaryFile
    '''
    ps1_bonds = []
    ps2_bonds = []
    with open(ps1, 'r') as instream:
        read = csv.reader(instream, delimiter="\t")
        for row in read:
            if row[0] != 'pdb':
                res1, res2 = sorted([row[1], row[2]])
                ps1_bonds.append(res1 + '/' + res2 + '(' + row[4] + ')')
    with open(ps2, 'r') as instream:
        read = csv.reader(instream, delimiter="\t")
        for row in read:
            if row[0] != 'pdb':
                res1, res2 = sorted([row[1], row[2]])
                ps2_bonds.append(res1 + '/' + res2 + '(' + row[4] + ')')
    overlap = [ bond for bond in list(set(ps1_bonds + ps2_bonds)) if (bond in ps1_bonds and bond in ps2_bonds) ]
    unique = [ bond for bond in list(set(ps1_bonds + ps2_bonds)) if (bond in ps1_bonds and bond not in ps2_bonds) ]
    if os.path.isfile(compareSummaryFile) is False:
        with open(compareSummaryFile, 'a') as outstream:
            wr = csv.writer(outstream, delimiter="\t")
            wr.writerow(['pdb1', 'pdb2', 'type', 'overlap_n', 'overlap', 'unique_in_pdb1_n', 'unique_in_pdb1'])
    with open(compareSummaryFile, 'a') as outstream:
        wr = csv.writer(outstream, delimiter="\t")
        wr.writerow([ps1[:-8], ps2[:-8], 'pi-stack', str(len(overlap)), ';'.join(overlap), str(len(unique)), ';'.join(unique)])
"""
