'''
invoking PIFS
'''

import argparse
import csv
import itertools
import os
import pifsUtils
import subprocess

if __name__=="__main__":
	parser = argparse.ArgumentParser(description="PIFS: Protein Interface Fingerprinting using spectral decomposition. This script deals with extracting residue networks from protein structure. It invokes POPSCOMP for interface calculation and uses biopython functionalities to manipulate structures and calculate Calpha-Calpha distances.")
	parser.add_argument('--structList', help = 'a list of structures and chains to analyse.', type = str, required = True)
	parser.add_argument('--inDir', help = 'directory where input PDB files (or their gzipped versions) are.', type = str, required = True)
	parser.add_argument('--networkOutDir', help = 'directory to network file output', type = str, required = True)
	parser.add_argument('--popscompOutDir', help = 'directory to POPSCOMP outputs', type = str, required = True)
	parser.add_argument('--pifsDir', help = 'directory to PIFS network feature summaries and census', type = str, required = True)
	parser.add_argument('--pifsTag', help = 'name of this dataset, to be noted in the filenames of network feature summaries.', type = str, required = True)
	parser.add_argument('--popsDir', help = 'directory where the POPS executable is.', default = "/usr/local/bin", type = str)
	parser.add_argument('--popscompDir', help = 'directory where the POPSCOMP executable is.', default = "~/Documents/POPSCOMP/POPSCOMP", type = str)
	parser.add_argument('--errorFile', help = 'filename for log file of errors. Note: NO ABSOLUTE PATH since It will go in outDir.', type = str, default = "error.log")
	parser.add_argument('--distCutoff', help = 'distance cut-off to extract residue neighbours.', default = 6.0, type = float)
	parser.add_argument('--atomType', help = 'atom type based on which distances between residues are calculated.', default = "CA", type = str)
	parser.add_argument('--gzipped', help = 'whether input PDBs are gzipped.', type = bool, default = True)
	args = parser.parse_args()

	structDict = {}
	if os.path.isfile(args.structList):
		with open(args.structList, 'r') as instream:
			read = csv.reader(instream, delimiter="\t")
			for row in read:
				print row[0]
				protein = row[1].split(',')
				other = row[2].split(',')
				if len(protein) + len(other) <= 20:
					structDict[row[0]] = {"protein": protein, "other": other}

	print "To analyse a total of  " + str(len(structDict.keys())) + " structures."

	if args.gzipped  is True:
		inputExt = "pdb.gz"
	else:
		inputExt = "pdb"

	print "Input structure files are of format: " + inputExt + ' .'

	print "Now run POPSCOMP."

	for pdbcode, chains in structDict.items():
		pdbdir = pdbcode[1:3]
		if os.path.isfile(os.path.join(args.pifsDir, pdbcode + '_TopolFeat.pifs')) is False:
			run_popscomp = pifsUtils.runPOPSCOMP(args.popsDir, args.popscompDir, \
				pdbcode, args.inDir + pdbdir + '/', args.popscompOutDir + pdbdir + '/', \
				args.errorFile, inputExt = inputExt)
			if run_popscomp == 'Done' and \
				os.path.isfile(args.popscompOutDir + pdbdir + '/popscomp.' + pdbcode + '.pdb.all.diff'):
				# get interface residues
				allNeighbours = list()
				for chains in list(itertools.product(chains["protein"], chains["other"])):
					protein_chain = chains[0]
					chain1, chain2 = sorted(chains)
					interface = pifsUtils.extractInteractRes(pdbcode, chain1, chain2, args.popscompOutDir + pdbdir + '/')
					if interface == []:
					# sometimes the chain pair in POPSCOMP is not sorted, in this case
					# try swapping and see if this helps!
						interface = pifsUtils.extractInteractRes(pdbcode, chain2, chain1, args.popscompOutDir + pdbdir + '/')
					if interface != []:
						print "Extract neighbour and calculate network ..."
						if interface[0].chain_a == protein_chain:
							resid = list([interaction.chain_a_res for interaction in interface][0])
						else:
							resid = list([interaction.chain_b_res for interaction in interface][0])
						# extract neigbours
						neighbours = pifsUtils.extractNeighbours(args.inDir + pdbdir + '/' + pdbcode + '.' + inputExt, \
							protein_chain, resid, distCutoff = args.distCutoff)
						if isinstance(neighbours, str) is True:
							print neighbours # So the error will get printed out
						else:
							allNeighbours = allNeighbours + neighbours
						# calculate distances
						distanceDict = pifsUtils.calcDistance(neighbours, \
							atomType = args.atomType)
						# write output network flat-file
						sifFile = args.networkOutDir + pdbcode + '-' + chain1 + 'vs' \
							+ chain2 + '_network.txt'
						if pifsUtils.writeSif(distanceDict, sifFile) is "Done":
							print pdbcode + ':' + chain1 + 'vs' + chain2 + ' residue network extracted.'
				# calculate network based on entire structure
				print "Now calculate a network on all interfaces in the structure ..."
				allDistances = pifsUtils.calcDistance(list(set(allNeighbours)), \
					atomType = args.atomType)
				allSifFile = args.networkOutDir + pdbcode + '_network.txt'
				if pifsUtils.writeSif(allDistances, allSifFile, directNetwork = True, distCutoff = args.distCutoff) is "Done":
					print pdbcode + ' residue network extracted.'
					# extract features
			                pifsUtils.extractOrganicFeatures(allSifFile, \
			                        os.path.join(args.pifsDir, pdbcode + '_TopolFeat.pifs'))
			
				os.system(''.join(['rm ', pdbcode, '.pdb']))
	# perform a census of features
#	pifsUtils.census(args.pifsDir, args.pifsTag)
	print 'Exit'
