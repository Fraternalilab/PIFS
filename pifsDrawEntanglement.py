'''
function to draw entanglement; extension to PyMOL
uses plane.py from PyMOLWiki
source it in pymol: do `run pifsDrawEntanglement.py` befor using
'''
import plane

def entangle(selection, p1, p2, p3, n1, chain_prot='A', chain_NA='B', n2 = None, plane_color=[1.0, 0.5, 0.0], vector_color=[0.0, 0.75, 1.0]):
    '''
DESCRIPTION

    Depict an 'entanglement' with a given list of residues which define a plane, and a single (or pair of) nucleic acid base(s) defining a vector.

USAGE

    entangle [ selection , p1, p2, p3, n1, chain_prot , chain_NA , [, n2 [, plane_color [, vector_color]]]]

ARGUMENTS

    selection = str: name of object

    p1 = int: residue number of 1st residue on protein to define plane

    p2 = int: residue number of 2nd residue on protein to define plane

    p3 = int: residue number of 3rd residue on protein to define plane

    n1 = int: residue number of 1st base on nucleic acid chain to define vector. If only n1 is indicated but not n2,
    the vector is defined based on this base only.

    n2 = int: residue number of 2nd base on nucleic acid chain to define vector. 
    The 5' to 3' from n1 to n2 defines the vector.

    chain_prot = str: chain ID for protein in selection (default: 'A')

    chain_NA = str: chain ID for nucleic acid in selection (default: 'B')

    plane_color = list of 3 floats: color of the plane. Defined as [R, G, B] where each of R, G and B are in the range [0, 1]. (Default: orange [1.0, 0.5, 0.0])
    
    vector_color = list of 3 floats: color of the vector. Defined as [R, G, B] where each of R, G and B are in the range [0, 1]. (Default: skyblue [0.0, 0.75, 1.0]) If only a single-nucleic acid base is selected this will not be depicted with a line.

    '''
    cmd.select('p1', '/' + selection + '//' + chain_prot + '/' + p1 + '/CA')
    cmd.select('p2', '/' + selection + '//' + chain_prot + '/' + p2 + '/CA')
    cmd.select('p3', '/' + selection + '//' + chain_prot + '/' + p3 + '/CA')
    p1_com = cmd.centerofmass('p1')
    p2_com = cmd.centerofmass('p2')
    p3_com = cmd.centerofmass('p3')
    plane.make_plane_points(chain_prot + '-' + ','.join([str(i) for i in [p1, p2, p3]]), p1_com, p2_com, p3_com, settings = {'ALPHA':0.4, 'COLOR':plane_color})
    cmd.delete('p1')
    cmd.delete('p2')
    cmd.delete('p3')

    if n2 is not None:
        # generate a line
        cmd.select('l1', '/' + selection + '//' + chain_NA + '/' + n1)
        cmd.select('l2', '/' + selection + '//' + chain_NA + '/' + n2)
        l1_com = cmd.centerofmass('l1')
        l2_com = cmd.centerofmass('l2')
        plane.cgo_arrow(chain_NA + '-' + ','.join([str(i) for i in [n1, n2]]), l1_com, l2_com, color1=vector_color, color2=vector_color )
        cmd.delete('l1')
        cmd.delete('l2')

    cmd.orient( chain_prot + '-' + ','.join([str(i) for i in [p1, p2, p3]]) )
    #cmd.zoom('/' + selection + '//' + chain_prot)

cmd.extend('entangle', entangle)


