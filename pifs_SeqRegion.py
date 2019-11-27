'''
invoking PIFS with indication of sequence region of interest 
(instead of computing interfaces and fetching residue network local to the interfaces)
'''

import argparse
import csv
import itertools
import os
import pifsUtils
import subprocess

if __name__=="__main__":
	parser = argparse.ArgumentParser(description="PIFS: Protein Interface Fingerprinting with Subgraphs. This script deals with extracting residue networks from protein structure. It takes input for what sequence regions are of interest for each structure and uses biopython functionalities to manipulate structures and calculate Calpha-Calpha distances.")
	parser.add_argument('--structList', help = 'a list of structures and residue ranges to analyse.', type = str, required = True)
        parser.add_argument('--inDir', help = 'directory where input PDB files (or their gzipped versions) are.', type = str, required = True)
	parser.add_argument('--networkOutDir', help = 'directory to network file output', type = str, required = True)
	# parser.add_argument('--pifsDir', help = 'directory to PIFS network feature summaries and census', type = str, required = True)
	# parser.add_argument('--pifsTag', help = 'name of this dataset, to be noted in the filenames of network feature summaries.', type = str, required = True)
	parser.add_argument('--errorFile', help = 'filename for log file of errors. Note: NO ABSOLUTE PATH since It will go in outDir.', type = str, default = "error.log")
	parser.add_argument('--distCutoff', help = 'distance cut-off to extract residue neighbours.', default = 6.0, type = float)
	parser.add_argument('--atomType', help = 'atom type based on which distances between residues are calculated.', default = "CA", type = str)
	parser.add_argument('--gzipped', help = 'whether input PDBs are gzipped.', type = str, default = True)
	args = parser.parse_args()

	structDict = {}
	if os.path.isfile(args.structList):
		with open(args.structList, 'r') as instream:
			read = csv.reader(instream, delimiter="\t")
			for row in read:
				print row[0]
				protein = row[1].split(',')
				region = row[2].split(',')
				structDict[row[0]] = {"file": os.path.join(args.inDir, row[0]), "protein": protein, "range": region}

	print "To analyse a total of  " + str(len(structDict.keys())) + " structures."

	if args.gzipped != "True":
		inputExt = "pdb"
	else:
		inputExt = "pdb.gz"

	print "Input structure files are of format: " + inputExt + ' .'

	for pdbcode, chains in structDict.items():
		# pdbdir = pdbcode[1:3]
		if os.path.isfile(os.path.join(args.pifsDir, pdbcode + '_TopolFeat.pifs')) is False:
			allNeighbours = list()
			protein_chain = chains['protein'][0]
			resid = range(int(chains['range'][0]), int(chains['range'][1]) + 1)
			# extract neigbours
			neighbours = pifsUtils.extractNeighbours(chains['file'], \
				protein_chain, resid, distCutoff = float(args.distCutoff))
			if isinstance(neighbours, str) is True:
				print neighbours # So the error will get printed out
			else:
				allNeighbours = allNeighbours + neighbours
				# calculate distances
				distanceDict = pifsUtils.calcDistance(neighbours, \
					atomType = args.atomType)
				# write output network flat-file
				if protein_chain == ' ':
					protein_chain = ''
				#sifFile = args.networkOutDir + pdbcode + '-' + protein_chain + '_network.txt'
				#if pifsUtils.writeSif(distanceDict, sifFile) is "Done":
				#	print pdbcode + ':' + protein_chain + ' residue network extracted.'
			# calculate network based on entire structure
			print "Now calculate a network on all interfaces in the structure ..."
			allDistances = pifsUtils.calcDistance(list(set(allNeighbours)), \
				atomType = args.atomType)
			allSifFile = args.networkOutDir + args.pifsTag + "_" + pdbcode + '_network.txt'
			if pifsUtils.writeSif(allDistances, allSifFile, directNetwork = True, distCutoff = float(args.distCutoff)) is "Done":
				print pdbcode + ' residue network extracted.'
				# extract features
			        # pifsUtils.extractOrganicFeatures(allSifFile, \
			        #         os.path.join(args.pifsDir, pdbcode + '_TopolFeat.pifs'))
			# os.system(''.join(['rm ', pdbcode, '.pdb']))
	# perform a census of features
#	pifsUtils.census(args.pifsDir, args.pifsTag)
	print 'Exit'
