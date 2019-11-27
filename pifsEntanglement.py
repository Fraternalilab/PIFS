'''
compute vectors and planes to detect entanglements in protein-nucleic acid grafts
'''
from Bio.PDB import *
import csv
import itertools
import sys
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

def definePlane_Prot(pdbfile, pdbname, resrange, chain='A', atomType = 'CA'):
    '''
    return Plane object defined by atoms.
    All possible combinations of triple C-alphas.
    '''
    pdb = PDBcoord()
    pdb.readPDB(pdbfile, pdbname)
    resis = range(resrange[0], resrange[1])
    # triples = list(itertools.combinations(resis, 3))
    triples = [(resrange[0], i, resrange[1]) for i in resis]
    planes = dict()
    for triple in triples:
        atomCoords = [pdb.getAtomsCoord(i, chain, atomType) for i in triple]
        if any(type(x) is str for x in atomCoords):
            planes[str(triple)] = None
        elif Point3D.are_collinear(Point3D(atomCoords[0]), Point3D(atomCoords[1]), Point3D(atomCoords[2])) is False: 
            planes[str(triple)] = atomCoords
    return planes    

def defineVector_NA(pdbfile, pdbname, resrange, chain='B'):
    '''
    return line3D (aka vector) object defined by the nucleic acid molecule.
    defined as vector from C5' (i.e. at the (deoxy)ribose) (at pos n-1) to C3' (at (deoxy)ribose of pos n)
    '''
    pdb = PDBcoord()
    pdb.readPDB(pdbfile, pdbname)
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
        line_start =  pdb.getAtomsCoord(base, chain, "C3'")
        line_end =  pdb.getAtomsCoord(base, chain, "C6")
    elif baseType in ['U', 'T', 'C', 'DT', 'DC']:
        line_start =  pdb.getAtomsCoord(base, chain, "C3'")
        line_end =  pdb.getAtomsCoord(base, chain, "C4")
    else:
        line_start = 'None'
        line_end = 'None'
    return base, [line_start, line_end]

def entangle_angle(vector, plane):
    '''
    test whether the nucleic acid is embedded/wrapped/entangled in a specified loop
    a position (base) of the nucleic acid molecule is specified as the vector and
    part of the loop is specified as a plane.
    Amounts to calculating the angle between the vector and the plane
    '''
    if any(type(x) is str for x in vector) or any(type(x) is str for x in plane):
        return None
    vector = Line3D(Point3D(vector[0], evaluate=False), Point3D(vector[1], evaluate=False))
    plane = Plane(Point3D(plane[0], evaluate=False), Point3D(plane[1], evaluate=False), \
        Point3D(plane[2], evaluate=False))
    angle_rad = plane.angle_between(vector)
    return angle_rad

def point_relative_to_plane(plane, point):
    '''
    calculate dot product between v (the normal of the plane), and
    (p - p0) where p0 is the projection of a point p on the plane.
    '''
    p = point
    p0 = [x for x in plane.projection(p)]
    v = [x for x in plane.normal_vector]
    v_p = [p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]]
    dot_product = np.dot(v, v_p)
    return dot_product

def barycentric(a, b, c, p):
    '''
    calculate barycentric coordinates of a point p based on 
    points a, b, c which defines a triangle on a plane.
    solution is based on Cramer's rule solving a system of 
    2 x 2 linear equations following the definition of 
    barycentric coordinate:

    _[ P = uA + vB + wC
     [ u + v + w = 1

    which rearranges and simplifies to solving for v and w : 

    _[ (P - A) . (B - A) = v[(B - A) . (B - A)] + w [(C - A) . (B - A)]
     [ (P - A) . (C - A) = v[(B - A) . (C - A)] + w [(C - A) . (C - A)]
    (where '.' signifies a dot product)

    and then applying Cramer's rule.
    '''
    v0 = b - a
    v1 = c - a
    v2 = p - a
    A = np.array([[np.dot(v0, v0), np.dot(v1, v0)], [np.dot(v0, v1), np.dot(v1, v1)]], dtype='float')
    B = np.array([np.dot(v2, v0), np.dot(v2, v1)], dtype='float')
    v, w = np.linalg.solve(A, B)
    # d00 = float(np.dot(v0, v0))
    # d01 = float(np.dot(v0, v1))
    # d11 = float(np.dot(v1, v1))
    # d20 = float(np.dot(v2, v0))
    # d21 = float(np.dot(v2, v1))
    # denom = d00 * d11 - d01 * 01
    # v = (d11 * d20 - d01 * d21) / denom
    # w = (d00 * d21 - d01 * d20) / denom
    u = 1.0 - v - w
    return u, v, w

def entangle_contains(vector, plane):
    '''
    test whether the nucleic acid is embedded/wrapped/entangled in a specified loop
    a position (base) of the nucleic acid molecule is specified as the vector and
    part of the loop is specified as a plane.
    this functions project both points which define the vector onto the plane and test whether 
    the intersection point of the vector and the plane is within the triangle which defines the plane
    '''
    # project both points defining the vector (v1 and v2) to the plane
    if any(type(x) is str for x in vector) or any(type(x) is str for x in plane):
        return None
    plane3d = Plane(Point3D(plane[0]), Point3D(plane[1]), Point3D(plane[2]))
    vector3d = Line3D(Point3D(vector[0]), Point3D(vector[1]))
    v_p = plane3d.intersection(vector3d)[0] # intersection point
    # test whether each of v0, v1 which define the vector is below / above the plane. 
    # Carry forward if one below and the other above.
    dp_v0 = point_relative_to_plane(plane3d, Point3D(vector[0]))
    dp_v1 = point_relative_to_plane(plane3d, Point3D(vector[1]))
    if (dp_v0 < 0 and dp_v1 > 0) or (dp_v0 > 0 and dp_v1 < 0):
        # calculate barycentric coordinate of v_p relative to the triangle.
        barycoord = barycentric(np.array(plane[0]), np.array(plane[1]), \
            np.array(plane[2]), np.array([x for x in v_p]))
        if all([( 0 <= i <= 1 ) for i in barycoord]):
            return True
        else:
            return False
    else:
        return False

def angleToScore(radian):
    '''
    convert angle A (in radian) to 1-cosA which ranges from 0 to 1 in -pi/2 <= A <- pi/2.
    The higher this value, the 
    '''
    if radian is None:
        return None
    return 1 - math.cos(radian)

def chain_entangle(pdb, name, resrangeProt, resrangeNA, chainProt, chainNA, outname):
    '''
    the pipeline to test for chain-entanglements (i.e. entanglements involving the nucleic-acid chain (rather than individual base) and write output file
    '''
    planes = definePlane_Prot(pdb, name, resrangeProt, chain=chainProt, atomType='CA')
    vectors = defineVector_NA(pdb, name, resrangeNA, chain=chainNA)
    if os.path.isfile(outname) is False:
        with open(outname, 'w') as outstream:
            wr = csv.writer(outstream, delimiter="\t")
            wr.writerow(['pdb', 'chain_pos', 'plane_pos', 'state', 'angle', 'score'])
    for keyNA, vector in vectors.iteritems():
        for keyProt, plane in planes.iteritems():
            if plane is None:
                state = None
                angle = None
                score = None
            else:
                state = entangle_contains(vector, plane)
                angle = entangle_angle(vector, plane)
                score = angleToScore(entangle_angle(vector, plane))
            with open(outname, 'a') as outstream:
                wr = csv.writer(outstream, delimiter="\t")
                if angle is None:
                    angle = str(None)
                else:
                    angle = str(math.degrees(angle))
                wr.writerow([name, keyNA, keyProt, str(state), angle, str(score)])
    return None        

def base_entangle(pdb, name, resrangeProt, resrangeNA, chainProt, chainNA, outname):
    '''
    the pipeline to test for base-entanglements (i.e. entanglements involving individual bases) and write output file
    '''
    planes = definePlane_Prot(pdb, name, resrangeProt, chain=chainProt, atomType='CA')
    vectors = [defineVector_base(pdb, name, n, chain=chainNA) for n in range(resrangeNA[0], resrangeNA[1]+1)]
    if os.path.isfile(outname) is False:
        with open(outname, 'w') as outstream:
            wr = csv.writer(outstream, delimiter="\t")
            wr.writerow(['pdb', 'base_pos', 'plane_pos', 'state', 'angle', 'score'])
    for (base, vector) in vectors:    
        for keyProt, plane in planes.iteritems():
            if plane is None:
                state = None
                angle = None
                score = None
            else:
                state = entangle_contains(vector, plane)
                angle = entangle_angle(vector, plane)
                score = angleToScore(entangle_angle(vector, plane))
            with open(outname, 'a') as outstream:
                wr = csv.writer(outstream, delimiter="\t")
                if angle is None:
                    angle = str(None)
                else:
                    angle = str(math.degrees(angle))
                wr.writerow([name, base, keyProt, str(state), angle, str(score)])
    return None
