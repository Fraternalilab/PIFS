import pifsDraw
import numpy as np

#planes = pifsDraw.definePlane_Prot('../pifs_cases/grafts/analysis/packbind_A3D-N_graft5keg_dna_0074_0001.pdb', 'A3D-N', [58, 78], chain='A', atomType='CA')
#vectors = pifsDraw.defineVector_NA('../pifs_cases/grafts/analysis/packbind_A3D-N_graft5keg_dna_0074_0001.pdb', 'A3D-N', [1, 4], chain='B')
planes = pifsDraw.definePlane_Prot('../pifs_cases/grafts/analysis/packbind_A3A_graft5keg_dna_0008_0001.pdb', 'A3A', [58, 70], chain='A', atomType='CA')
vectors = pifsDraw.defineVector_NA('../pifs_cases/grafts/analysis/packbind_A3A_graft5keg_dna_0008_0001.pdb', 'A3A', [1, 4], chain='B')

print pifsDraw.defineVector_base('../pifs_cases/grafts/analysis/packbind_A3A_graft5keg_dna_0008_0001.pdb', 'A3A', 2, 'DT', chain='B')

#print pifsDraw.embedding_contains(vectors['1-2'], planes['(64, 68, 70)'])
#print pifsDraw.embedding_contains(vectors['2-3'], planes['(64, 68, 70)'])
#print pifsDraw.embedding_contains(vectors['3-4'], planes['(64, 68, 70)'])
#print pifsDraw.embedding_contains(vectors['2-3'], planes['(58, 61, 64)'])
#print '------'
#print pifsDraw.embedding_angle(vectors['1-2'], planes['(64, 68, 70)'])
#print pifsDraw.embedding_angle(vectors['2-3'], planes['(64, 68, 70)'])
#print pifsDraw.embedding_angle(vectors['3-4'], planes['(64, 68, 70)'])
#print pifsDraw.embedding_angle(vectors['2-3'], planes['(58, 61, 64)'])

