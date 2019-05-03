'''
draw C-alpha network in PyMol and compute vectors and planes
'''
from Bio.PDB import *
import csv
import itertools
import sys
sys.path.append('/home/josephn/pymol-open-source-v2.1.0/lib/python/')
import os
import cmd
from sympy import Point3D, Line3D, Plane
import math
import numpy as np
from shapely.geometry import Point, Polygon

def drawCaNetwork(C_alphas, pdb, outfile, cutoff=6.0):
    '''
    draw C-alpha network
    '''
    # read the structure
    cmd.load(pdb, 'structure')
    # select each Calpha
    for res in C_alphas:
        # get chain, resn and resi
        res_split = res.split(':')
        resi=res_split[1]
        res_split = res_split.split('-')
        chain=res_split[1]
        resn=res_split[0]
        cmd.select('ca_sele', '/structure//' + chain + '/' + resi + '/' + 'CA', merge=1)
    cmd.create('ca', 'ca_sele')
    cmd.delete('structure')
    cmd.delete('ca_sele')
    # configure the appearance
    cmd.show_as('spheres', 'ca')
    cmd.alter('/ca////', vdw=3.0)
    cmd.bg_color(color="white")
    cmd.color('/ca////', 'deepsalmon')
    cmd.set('label_color', 'black')
    cmd.set('label_font_id', 9)
    cmd.label('/ca////', '"%s%s" %(resn, resi)' )
    # generate the network
    cmd.distance('dist', '/ca////', cutoff=cutoff)
    cmd.save(outfile, 'all')
    cmd.quit()
    return None

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

def definePlane_Prot(pdbfile, pdbname, resrange, chain='A', atomType = 'CA'):
    '''
    return Plane object defined by atoms.
    All possible combinations of triple C-alphas.
    '''
    pdb = PDBcoord()
    pdb.readPDB(pdbfile, pdbname)
    resis = range(resrange[0], resrange[1]+1)
    triples = list(itertools.combinations(resis, 3))
    planes = dict()
    for triple in triples:
        atomCoords = [pdb.getAtomsCoord(i, chain, atomType) for i in triple]
        planes[str(triple)] = atomCoords
    return planes    

def defineVector_NA(pdbfile, pdbname, resrange, chain='B'):
    '''
    return line3D (aka vector) object defined by the nucleic acid molecule.
    defined as vector from C5' (i.e. at the (deoxy)ribose) (at pos n-1) to C3' (at (deoxy)ribose of pos n)
    '''
    pdb = PDBcoord()
    pdb.readPDB(pdbfile, pdbname)
    resis = range(resrange[0], resrange[1]+1)
    lines = dict()
    for n in range(resrange[0], resrange[1]):
        line_start = pdb.getAtomsCoord(n, chain, "C5'")
        line_end = pdb.getAtomsCoord(n+1, chain, "C3'")
        lines['-'.join([str(n), str(n+1)])] = [line_start, line_end]
    return lines

def defineVector_base(pdbfile, pdbname, base, chain='B'):
    '''
    return line3D (aka vector) object defined by the indicated base.
    For purine (A/G), defined as the vector from N9 to C6
    For pyrimidine (U/T/C), defined as the vector from N1 to C4
    '''
    pdb = PDBcoord()
    pdb.readPDB(pdbfile, pdbname)
    baseType = pdb.getResName(base, chain)
    if baseType in ['A', 'G', 'DA', 'DG']:
        line_start =  pdb.getAtomsCoord(base, chain, "N9")
        line_end =  pdb.getAtomsCoord(base, chain, "C6")
    elif baseType in ['U', 'T', 'C', 'DT', 'DC']:
        line_start =  pdb.getAtomsCoord(base, chain, "N1")
        line_end =  pdb.getAtomsCoord(base, chain, "C4")
    out = dict()
    out[str(base)] = [line_start, line_end]
    return out 

def knot_angle(vector, plane):
    '''
    test whether the nucleic acid is embedded/wrapped/knotted in a specified loop
    a position (base) of the nucleic acid molecule is specified as the vector and
    part of the loop is specified as a plane.
    Amounts to calculating the angle between the vector and the plane
    '''
    vector = Line3D(Point3D(vector[0], evaluate=False), Point3D(vector[1], evaluate=False))
    plane = Plane(Point3D(plane[0], evaluate=False), Point3D(plane[1], evaluate=False), \
        Point3D(plane[2], evaluate=False))
    angle_rad = plane.angle_between(vector)
    return angle_rad

def knot_contains(vector, plane):
    '''
    test whether the nucleic acid is embedded/wrapped/knotted in a specified loop
    a position (base) of the nucleic acid molecule is specified as the vector and
    part of the loop is specified as a plane.
    this functions project both points which define the vector onto the plane and test whether the polygon which defines the plane contain both of these points
    '''
    # project both points defining the vector (v1 and v2) to the plane
    plane3d = Plane(Point3D(plane[0]), Point3D(plane[1]), \
        Point3D(plane[2]))
    vector = Line3D(Point3D(vector[0]), Point3D(vector[1]))
    v_p = plane3d.intersection(vector) # intersection point
    P_t = Polygon([[plane[0][0], plane[0][1]], [plane[1][0], plane[1][1]], [plane[2][0], plane[2][1]]])
    return P_t.contains(Point(v_p.x, v_p.y))

def angleToScore(radian):
    '''
    convert angle A (in radian) to 1-cosA which ranges from 0 to 1 in -pi/2 <= A <- pi/2.
    The higher this value, the 
    '''
    return 1 - math.cos(radian)

def chain_knot(pdb, name, resrangeProt, resrangeNA, chainProt, chainNA, outname):
    '''
    the pipeline to test for chain-knots (i.e. knots involving the nucleic-acid chain (rather than individual base) and write output file
    '''
    planes = definePlane_Prot(pdb, name, resrangeProt, chain=chainProt, atomType='CA')
    vectors = defineVector_NA(pdb, name, resrangeNA, chain=chainNA)
    with open(outname, 'w') as outstream:
        wr = csv.writer(outstream, delimiter="\t")
        wr.writerow(['chain_pos', 'plane_pos', 'state', 'angle', 'score'])
    for keyNA, vector in vectors.iteritems():
        for keyProt, plane in planes.iteritems():
            state = knot_contains(vector, plane)
            angle = knot_angle(vector, plane)
            score = angleToScore(knot_angle(vector, plane))
            wr.writerow([keyNA, keyProt, str(state), str(angle), str(score)])
    return None        

def base_knot(pdb, name, resrangeProt, resrangeNA, chainProt, chainNA, outname):
    '''
    the pipeline to test for base-knots (i.e. knots involving individual bases) and write output file
    '''
    planes = definePlane_Prot(pdb, name, resrangeProt, chain=chainProt, atomType='CA')
    vectors = [defineVector_base(pdbfile, pdbname, n, chain=chainNA) for n in range(resrangeNA[0], resrange[1]+1)]
    with open(outname, 'w') as outstream:
        wr = csv.writer(outstream, delimiter="\t")
        wr.writerow(['base_pos', 'plane_pos', 'state', 'angle', 'score'])
    for keyNA, vector in vectors.iteritems():
        for keyProt, plane in planes.iteritems():
            state = knot_contains(vector, plane)
            angle = knot_angle(vector, plane)
            score = angleToScore(knot_angle(vector, plane))
            wr.writerow([keyNA, keyProt, str(state), str(angle), str(score)])
    return None
