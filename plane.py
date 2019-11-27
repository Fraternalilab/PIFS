'''
Described at PyMOL wiki:
https://pymolwiki.org/index.php/Plane_Wizard

Authors : Troels Schwarz-Linnet
Date    : Dec 2016
Modified: From previous contributors. 

Modified by Joseph Ng July 2019 to plot only triangles instead of quadrilaterals (the 'plane' function)

'''
import sys
sys.path.append('/home/josephn/pymol-open-source-v2.1.0/lib/python/')
import pymol
from pymol import cmd, cgo, CmdException
from pymol.wizard import Wizard
from chempy import cpv
from pymol.cgo import COLOR, SPHERE, CYLINDER, BEGIN, TRIANGLE_STRIP, NORMAL, VERTEX, END, ALPHA

def makePrimitive(cgo, name):
    az = cmd.get('auto_zoom', quiet=1)
    cmd.set('auto_zoom', 0, quiet=1)
    cmd.load_cgo(cgo, name)
    cmd.set('auto_zoom', az, quiet=1)

def point(p):
    x, y, z = p
    return [COLOR, 1, 1, 1, SPHERE, float(x), float(y), float(z), 0.5]


def line(p1, p2):
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    return [CYLINDER, float(x1), float(y1), float(z1), float(x2), float(y2), float(z2), 0.25, 1, 1, 1, 1, 1, 1]

def triangle(corner1, corner2, corner3, normal, settings):
    planeObj = []
    planeObj.extend(point(corner1))
    planeObj.extend(point(corner2))
    planeObj.extend(point(corner3))
    planeObj.extend(line(corner1, corner2))
    planeObj.extend(line(corner2, corner3))
    planeObj.extend(line(corner3, corner1))
    # planeObj.extend(line(corner4, corner1))

    # Make settings
    if 'ALPHA' in settings:
        planeObj.extend([ALPHA, settings['ALPHA']])

    if 'COLOR' in settings:
        planeObj.extend([COLOR, settings['COLOR'][0], settings['COLOR'][1], settings['COLOR'][2]])
    else:
        planeObj.extend([COLOR, 0.8, 0.8, 0.8]) # greyish

    planeObj.extend([BEGIN, TRIANGLE_STRIP])
    planeObj.append(NORMAL)

    if 'INVERT' in settings:
        if settings['INVERT']==True:
            planeObj.extend(cpv.negate(normal))
        else:
            planeObj.extend(normal)
    else:
        planeObj.extend(normal)

    for corner in [corner1, corner2, corner3, corner1]:
        planeObj.append(VERTEX)
        planeObj.extend(corner)
    planeObj.append(END)
    return planeObj


def planeFromPoints(p1, p2, p3, vm1=None, vm2=None, center=True, settings={}):
    v1 = cpv.sub(p1, p2)
    v2 = cpv.sub(p3, p2)
    normal = cpv.cross_product(v1, v2)

    if 'translate' in settings:
        vtran = cpv.scale(cpv.normalize(normal), settings['translate'])
        p1_t = cpv.sub(p1, vtran)
        p2_t = cpv.sub(p2, vtran)
        p3_t = cpv.sub(p3, vtran)
        print("New coordinates are:")
        print_info("New", p1_t, p2_t, p3_t)
        print("New coordinates are for normalized plane:")
        v1_t = cpv.normalize(cpv.sub(p1_t, p2_t))
        v2_t = cpv.normalize(cpv.sub(p3_t, p2_t))
        normal_t = cpv.normalize(cpv.cross_product(v1_t, v2_t))
        v2_t = cpv.normalize(cpv.cross_product(normal_t, v1_t))
        p1_t2 = cpv.add(v1_t, p2_t)
        p3_t2 = cpv.add(v2_t, p2_t)
        print_info("Newnormal", p1_t2, p2_t, p3_t2)

    if vm1!=None:
        v1 = cpv.scale(cpv.normalize(v1), vm1)
    if vm2!=None:
        v2 = cpv.scale(cpv.normalize(v2), vm2)

    centrum = p2
    if center:
        corner1 = p1
        corner2 = p2
        corner3 = p3
        #corner1 = cpv.add(cpv.add(centrum, v1), v2)
        #corner2 = cpv.sub(cpv.add(centrum, v1), v2)
        #corner3 = cpv.sub(cpv.sub(centrum, v1), v2)
        #corner4 = cpv.add(cpv.sub(centrum, v1), v2)
    else:
        corner1 = p1
        corner2 = p2
        corner3 = p3
        #corner1 = cpv.add(cpv.add(centrum, v1), v2)
        #corner2 = cpv.add(centrum, v1)
        #corner3 = centrum
        #corner4 = cpv.add(centrum, v2)

    return triangle(corner1, corner2, corner3, normal, settings)


def print_info(name, coor1, coor2, coor3):
    cs1 = (map(float, [ '%.2f' % elem for elem in coor1 ]) )
    cs2 = (map(float, [ '%.2f' % elem for elem in coor2 ]) )
    cs3 = (map(float, [ '%.2f' % elem for elem in coor3 ]) )
    print("You can also use the function calls with these coordinates")
    print("plane.make_plane_points(name='%s', l1=%s, l2=%s, l3=%s)"%(name, cs1, cs2, cs3))


def make_plane(name,a1='(pk1)',a2='(pk2)',a3='(pk3)', vm1=None, vm2=None, center=True, makepseudo=True, settings={}):
    """
    DESCRIPTION
    Create a CGO plane from three atomic coordinates

    USAGE
    make_plane name, a1, a2, a3

    where each atom is a standard PyMOL selection (defaults to pk1,pk2 and pk3)
    """
    # get coordinates from atom selections
    coor1 = cmd.get_model(a1).get_coord_list()[0]
    coor2 = cmd.get_model(a2).get_coord_list()[0]
    coor3 = cmd.get_model(a3).get_coord_list()[0]

    # Help with alternative
    print_info(name, coor1, coor2, coor3)

    # Get the plane
    line(coor1, coor2)
    
    plane = planeFromPoints(p1=coor1, p2=coor2, p3=coor3, vm1=vm1, vm2=vm2, center=center, settings=settings)
    makePrimitive(plane, name)
    #cmd.show("cgo", "plane*")

    if makepseudo:
        cmd.pseudoatom("%s_%s"%(name, "l1"), color=settings['COLOR'], pos=coor1)
        cmd.pseudoatom("%s_%s"%(name, "l2"), color=settings['COLOR'], pos=coor2)
        cmd.pseudoatom("%s_%s"%(name, "l3"), color=settings['COLOR'], pos=coor3)
        cmd.show_as('spheres', "%s_%s"%(name, "l1"))
        cmd.show_as('spheres', "%s_%s"%(name, "l2"))
        cmd.show_as('spheres', "%s_%s"%(name, "l3"))

# Extend function to be called inside pymol
cmd.extend("make_plane", make_plane)

def make_plane_points(name,l1=None,l2=None,l3=None, vm1=None, vm2=None, center=True, makepseudo=True, settings={}):
    """
    DESCRIPTION
    Create a CGO plane from three atomic coordinates

    USAGE
    make_plane name, l1, l2, l3

    where each xys is a list with floats of x,y,z coordinates
    """
    if l1==None or l2==None or l3==None:
        print("Please provide a list of xyz floats for each 3 positions")
        return
    if type(l1) is not list or type(l2) is not list or type(l3) is not list:
        print(type(l1),type(l2),type(l3))
        print("Please provide 3 list of xyz floats for each 3 positions")
        return

    plane = planeFromPoints(p1=l1, p2=l2, p3=l3, vm1=vm1, vm2=vm2, center=center, settings=settings)
    makePrimitive(plane, name)

    if makepseudo:
        cmd.pseudoatom("%s_%s"%(name, "l1"), color=settings['COLOR'], pos=l1)
        cmd.pseudoatom("%s_%s"%(name, "l2"), color=settings['COLOR'], pos=l2)
        cmd.pseudoatom("%s_%s"%(name, "l3"), color=settings['COLOR'], pos=l3)
        cmd.show_as('spheres', "%s_%s"%(name, "l1"))
        cmd.show_as('spheres', "%s_%s"%(name, "l2"))
        cmd.show_as('spheres', "%s_%s"%(name, "l3"))
    
    cmd.delete("%s_%s"%(name, "l1"))
    cmd.delete("%s_%s"%(name, "l2"))
    cmd.delete("%s_%s"%(name, "l3"))
    return None  

# Extend function to be called inside pymol
cmd.extend("make_plane_points", make_plane_points)

'''
cgo_arrow function

http://pymolwiki.org/index.php/cgo_arrow

(c) 2013 Thomas Holder, Schrodinger Inc.
Modified by Joseph Ng July 2019

License: BSD-2-Clause
'''

def cgo_arrow(name='', atom1='pk1', atom2='pk2', radius=0.5, gap=0.0, hlength=-1, hradius=-1,
              color1 = [0.0, 0.75, 1.0], color2 = [0.0, 0.75, 1.0]):
    '''
DESCRIPTION
    Create a CGO arrow between two picked atoms.
ARGUMENTS
    atom1 = string: single atom selection or list of 3 floats {default: pk1}
    atom2 = string: single atom selection or list of 3 floats {default: pk2}
    radius = float: arrow radius {default: 0.5}
    gap = float: gap between arrow tips and the two atoms {default: 0.0}
    hlength = float: length of head
    hradius = float: radius of head
    color = string: one or two color names {default: blue red}
    color1, color2 = list of 3 floats indicating color {default: [0.0, 0.75, 1.0] for skyblue}
    name = string: name of CGO object
    '''
    from chempy import cpv

    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)

    #try:
    #    color1, color2 = color.split()
    #except:
    #    color1 = color2 = color
    #color1 = list(cmd.get_color_tuple(color1))
    #color2 = list(cmd.get_color_tuple(color2))

    def get_coord(v):
        if not isinstance(v, str):
            return v
        if v.startswith('['):
            return cmd.safe_list_eval(v)
        return cmd.get_atom_coords(v)

    xyz1 = get_coord(atom1)
    xyz2 = get_coord(atom2)
    normal = cpv.normalize(cpv.sub(xyz1, xyz2))

    if hlength < 0:
        hlength = radius * 3.0
    if hradius < 0:
        hradius = hlength * 0.6

    if gap:
        diff = cpv.scale(normal, gap)
        xyz1 = cpv.sub(xyz1, diff)
        xyz2 = cpv.add(xyz2, diff)

    xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

    obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
          [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
          [1.0, 0.0]

    if not name:
        name = cmd.get_unused_name('arrow')

    cmd.load_cgo(obj, name)

cmd.extend('cgo_arrow', cgo_arrow)
