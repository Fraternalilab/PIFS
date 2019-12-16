"""
PIFS: Protein Interface Fingerprinting using spectral decomposition
"""
from Bio.PDB import *
import collections
import csv
import glob
import gzip
import numpy as np
import networkx as nx
import os
import re
import sys

#list of amino acids
aa = ['ALA', 'GLY', 'ILE', 'LEU', 'PRO', 'VAL', 'PHE', 'TRP', 'TYR', \
	  'ASP', 'GLU', 'ARG', 'HIS', 'LYS', 'SER', 'THR', \
	  'CYS', 'MET', 'ASN', 'GLN']

class Interaction:
	'''
	Class to store interface residues
	From Anna Laddach
	'''
	def __init__(self, pdb, chain_a, chain_b):
		self.pdb = pdb
		self.chain_a = chain_a
		self.chain_b = chain_b
		self.chain_a_res = set()
		self.chain_b_res = set()

	def add_residue(self, residue, chain):
		if chain == self.chain_a:
			self.chain_a_res.add(residue)
		elif chain == self.chain_b:
			self.chain_b_res.add(residue)
		else:
			print 'incorrect chain specified ' + chain + ' ' + self.pdb

def extractInteractRes(structName, chain1, chain2, \
	popscompOutDir, interaction_cutoff = 15.0):
	"""
	read POPSCOMP output and get Interface Res
	:param structName: name of the structure; as used in runPOPSCOMP, string
	:param chain1: chain ID for the first chain involved in the interface to probe, char
	:param chain2: chain ID for the second chain involved in the interface to probe, char
	:param popscompOutDir: directory where POPSCOMP outputs are found, string
	:param interaction_cutoff: cutoff in Angstrom-squared, the change in SASA upon complex formation greater than which are considered as interacting interface, float. Default = 15.0
	Modified from Anna Laddach
	"""
	pat = re.compile('\*\*\* Complex popscomp\.pair\.(.+)\.pdb-' \
		+ chain1 + '.NtoC:' + chain2 + '.NtoC(.*?)TOTAL DELTA', re.M|re.S)
	interactions = []
	infile = popscompOutDir + '/popscomp.' + structName + '.pdb.all.diff'
	with open(infile) as f:
		text = f.read()
		for match in pat.finditer(text):
			interactions.append(Interaction(match.group(1),chain1, \
				chain2))
			for line in match.group(0).split('\n'):
				info = line.split()
				if len(info)== 7 and info[0] != 'Resid' and float(info[5]) > interaction_cutoff:
					interactions[-1].add_residue(info[2],info[1])
	return interactions


def runPOPSCOMP(popsDir, popscompDir, struct, indir, outdir, error_file, \
	inputExt = "pdb.gz"):
	'''
	invoke POPSCOMP over all pdb/pdb.gz files in a given directory (indir)
	:param popsDir: directory where the pops executable is/linked to, string
	:param popscompDir: directory where popscomp_dom.pl is, string.
	:param struct:
	:param indir: directory where the structure file(s) is, string. The files are in either *.pdb or *.pdb.gz formats.
	:param outdir: directory of the destination of POPS/POPSCOMP output, stsring.
	:param error_file: filename of error outputs, string. Note the error file will be output to outdir.
	:param inputExpt: indicate whether files are in *.pdb or *.pdb.gz format. Default *.pdb.gz. A string.
	Modified from Anna Laddach
	'''
	if inputExt == 'pdb.gz':
		stripBack = -7
		zipflag = "--zipped"
	elif inputExt == 'pdb':
		stripBack = -4
		zipflag = ""
	else:
		return "inputExt must be either pdb.gz or pdb. Exit."

	with open(error_file, 'w') as error:
                for structure in glob.glob(os.path.join(indir, struct + '.' + inputExt)):
			#Only run POPSCOMP if output does not already exist.
#			if not glob.glob(outdir + '/popscomp*'  + struct + '*'):
			os.system('perl ' + popscompDir \
				+ '/popscomp_dom.pl --popsPath ' + popsDir \
				+ '/pops --absPath -u ' +  structure \
				+' --residueOut --totalOut --singleOut --pairwiseOut --diffOut ' \
				+ zipflag + ' --outputDir ' + outdir)

			#If output still doesn't exist record this in error file.
			if not glob.glob(outdir + '/popscomp*'  + struct + '*'):
				error.write('Error, no POPSCOMP output for ' + struct + '\n')

			#Cleanup unwanted POPSCOMP output.
			try:
				os.remove('sigma.out')
				os.remove('faillist')
				os.remove('stderr')
				os.remove('popsb.out')
				os.remove(name + '.pdb')
			except:
				pass

	return "Done"

def extractNeighbours(structure, chain, resid, \
	atomType = "CA", distCutoff = 6.0):
	'''
	Extract neighbours of a given list of residues in "resid", based on distances below a cut-off specified in "distCutoff"
	:param structure: filepath of the structure, string
	:param chain: chain, character
	:param resid: list of residues between which neighbours are extracted, list of integers
	:param atomType: type of atoms based on which to extract neighbours. a string. By default neighbours are extracted based on distances between C-alphas of the given residues
	:param distCutoff: distance cut-off in Angstrom subject to atomType, below which the pair of residues are considered neighbours. Default = 10.0
	'''

	if isinstance(distCutoff, float) is False:
		return "distCutoff must be a float. Exit."
	#filepath = structure.split('/')[:-1] + name + '.pdb'
	if glob.glob(structure):
		if "pdb.gz" in structure:
			stripBack = -7
			name = structure.split('/')[-1][:stripBack]
			structParser = PDBParser()
			pdb = structParser.get_structure(name, gzip.open(structure))
		elif "pdb" in structure:
			stripBack = -4
			name = structure.split('/')[-1][:stripBack]
			structParser = PDBParser()
			pdb = structParser.get_structure(name, structure)
		else:
			return "inputExt must be either pdb.gz or pdb. Exit."

		models = pdb.get_list()
		if len(models) > 0:
			pdb1 = pdb[0] # just take the first model for simplicity
			if pdb1.has_id(chain):
				pdbchain = pdb1[chain]
                                pdb_A = pdb1['A']
				all_ca = [atom for atom in Selection.unfold_entities(pdb_A,'A') if atom.name == atomType]
				calphas = list()
				neighbors = list()
				for res in resid:
					if pdbchain.has_id(int(res)):
						# get res in pdbchain
						resEntry = pdbchain[int(res)]
						# put itself into the neighbor list
						neighbors.append(resEntry)
						# fetch C-alphas and append
						if resEntry.has_id(atomType):
							ca = resEntry[atomType]
							calphas.append(ca)
					else:
						print "res %s not in given chain. Skip." % int(res)
						pass

				if len(calphas) > 0:
					print "Fetched a total of %s %s." % (len(calphas), atomType)
					search = NeighborSearch(all_ca)
					for calpha in calphas:
						neighbors = neighbors + \
							search.search(calpha.get_coord(), \
							distCutoff, level="R")  # Fetch on residue level
					neighbors = list(set(neighbors)) # unique entries
					print "A total of %s residues to consider." % len(neighbors)
				else:
					print "Error in extracting c-alphas from structure. Exit."

				return neighbors

			else:
				return "chain not in given pdb. Exit."
		else:
			return "pdb file has no models in it. Exit."

	else:
		return "File not found. Exit."

def calcDistance(resid, atomType = "CA"):
	'''
	Calculate all pairwise distances between residues given in "resid".
	:param resid: list of residues between which distances are calculated, array of residue objects output from extractNeighbours.
	:param atomType: type of atoms based on which to calculate the distance. a string of length = 1. By default distances are calculated between C-alphas of the given residues
	'''
	# calculate pairwise distance
	calphas = [atom for atom in Selection.unfold_entities(resid,'A') if atom.name == atomType]
	allDistances = list()
	allDistances_ids = list()
	for calpha in calphas:
		d = list()
		res1 = calpha.get_parent().get_full_id()
		res1name = calpha.get_parent().get_resname()
		res1num = res1[3][1]
		res1_id = ''.join([str(item) for item in res1[3]])
		res1 = re.sub('\s', '',':'.join([res1[2], res1_id]))
		for c in calphas:
			res2 = c.get_parent().get_full_id()
			res2name = c.get_parent().get_resname()
			res2num = res2[3][1]
			res2_id = ''.join([str(item) for item in res2[3]])
			res2 = re.sub('\s', '', ':'.join([res2[2], res2_id]))
			d.append({"distance": calpha - c,\
				"res1": res1name, "res_1_id": res1name + '-' + res1,\
				"res2": res2name, "res_2_id": res2name + '-' + res2, \
				"res1num": res1num, "res2num": res2num, \
				"uniqueID": '||'.join(sorted([res1name + '-' + res1, \
					res2name + '-' + res2]))})
		for entry in d:
			# keep track so that duplicated edges are filtered away
			# to filter away possible self-cycles
			if entry["uniqueID"] not in allDistances_ids \
				and entry['distance'] > 0 \
				and abs(entry['res1num'] - entry['res2num']) > 1:
				allDistances_ids.append(entry["uniqueID"])
				allDistances.append(entry)
	return allDistances

def writeSif(dictList, sif, nodeAKey = "res_1_id", nodeBKey = "res_2_id",
	nodeAattrKey = "res1", nodeBattrKey = "res2", weightKey = "distance", \
	sifSep = '\t', weightTransform = "inverse", directNetwork = True, distCutoff = 6.0):
	'''
	Write network to a text file.
	:param dictList: a list of dictionaries containing pairwise distances to be
	included in the network
	:param sif: filepath for the output text file.
	:param nodeAKey: key for the node A in each dict in dictList, a string.
	By default this is an unique identifier of the form "<residue>-<chain>:<resID>".
	:param nodeBKey: key for the node B in each dict in dictList, a string.
	:param nodeAattrKey: key for annotation of node A in each dict in dictList, a string.
	By default the annotation is the residue.
	:param nodeBattrKey: key for annotation of node B in each dict in dictList, a string.
	:param weightKey: key for the weights of the edge for each pairwise interaction, key in
	each dict in dictList, a string.
	:param sifSep: separator in the output file. Default = tab ("\t").
	:param weightTransform: should the weight values be transformed? Default "inverse"
	i.e. take the inverse of the values pointed by "weightKey".
	:param directNetwork: should the output network contains only direct edges
	(i.e. edges that fall within distanceCutoff). If false, network will contain every pairwise
	edges possible between all combinations of nodes. Default: True
	:param distCutoff: distCutoff based on which to exclude edges with distance
	beyond distCutoff. Only considered if directNetwork is True.
	TO DO: more methods of transformation (weightTransform)?
	'''
	if isinstance(sif, str) is False:
		return "sif must be a string. Exit."
	with open(sif, 'w') as outstream:
		wr = csv.writer(outstream, delimiter = sifSep)
		wr.writerow(["res_a", "res_a_id", "res_b", "res_b_id", \
			"distance", "weight"])
		for item in dictList:
			if isinstance(item,dict) is False:
				return "dictList contains items which are not dict. Exit."
			else:
				if nodeAKey not in item.keys():
					return "'%s' was given as nodeAKey but it was not found in dict keys. Exit." % nodeAKey
				elif nodeBKey not in item.keys():
					return "'%s' was given as nodeBKey but it was not found in dict keys. Exit." % nodeBKey
				elif nodeAattrKey not in item.keys():
					return "'%s' was given as nodeAattrKey but it was not found in dict keys. Exit." % nodeAattrKey
				elif nodeBattrKey not in item.keys():
					return "'%s' was given as nodeBattrKey but it was not found in dict keys. Exit." % nodeBattrKey
				elif weightKey not in item.keys():
					return "'%s' was given as weightKey but it was not found in dict keys. Exit." % weightKey
				else:
					if weightTransform == "inverse":
						weight = 1 / item[weightKey]
					if directNetwork is True:
						if item[weightKey] < distCutoff:
							wr.writerow([item[nodeAattrKey], item[nodeAKey], \
								item[nodeBattrKey], item[nodeBKey], \
								item[weightKey], weight])
					else:
							wr.writerow([item[nodeAattrKey], item[nodeAKey], \
								item[nodeBattrKey], item[nodeBKey], \
								item[weightKey], weight])
	return "Done"
"""
class OrganicFeatures:
	'''
	class to handle extraction of organic features of network
	'''
	def __init__(self):
		self.graph = nx.Graph()
		self.cycles = list()
		self.cliques = list()
		self.bridges = list()
		self.articulation_points = list()
		self.simple_path_fragment = list()

	def addGraph(self, networkFP):
		with open(networkFP, 'r') as instream:
			read = csv.reader(instream, delimiter='\t')
			for row in read:
				if row[0] != 'res_a':
					if row[1] not in list(self.graph.nodes()):
						self.graph.add_node(row[1], attr_dict={'res': row[0]})
					if row[3] not in list(self.graph.nodes()):
						self.graph.add_node(row[3], attr_dict={'res': row[2]})
					self.graph.add_edge(row[1], row[3], attr_dict={'dist': row[4]})
		self.simple_path_fragment = list(self.graph.edges())

	def getNode(self, node):
		return self.graph.node[node]['attr_dict']['res']

	def getEdge(self, node1, node2):
		edge = self.graph.edges[node1, node2]
		return [node1, node2, edge['attr_dict']['dist']]

	def findCycles(self, maxLength = 6):
		cycles = nx.cycle_basis(self.graph)
		for cycle in list(cycles):
			if len(cycle) <= maxLength:
				cycle.append(cycle[0])
				self.cycles.append(cycle)
				for i in range(len(cycle) - 1):
					try:
						self.simple_path_fragment.remove((cycle[i], cycle[i+1]))
					except ValueError:
						pass

	def findCliques(self, minSize = 4):
		cliques = nx.find_cliques(self.graph)
		if nx.graph_clique_number(self.graph, cliques=cliques) <4:
			pass
		else:
			for clique in list(cliques):
				if len(clique) >= 4:
					self.cliques.append(clique)
					for node1, node2 in list(itertools.combinations(clique, 2)):
						self.simple_path_fragment.remove((node1, node2))

	def findBridges(self):
		bridges = nx.bridges(self.graph)
		for bridge in list(bridges):
			self.bridges.append(bridge)
			self.simple_path_fragment.remove((bridge[0], bridge[1]))

	def findArticulationPoints(self):
		arPt = nx.articulation_points(self.graph)
		for point in list(arPt):
			self.articulation_points.append(point)

	def getSimplePathFragment(self):
		out = list()
		for edge in self.simple_path_fragment:
			out.append(self.getEdge(edge[0], edge[1]))
		return out

	def getCycles(self):
		cycles = list()
		for cycle in self.cycles:
			out = list()
			for i in range(len(cycle) - 1):
				out.append(self.getEdge(cycle[i], cycle[i+1]))
			cycles.append(out)
		return cycles

	def getCliques(self):
		cliques = list()
		for clique in self.cliques:
			out = list()
			for node1, node2 in list(itertools.combinations(clique, 2)):
				out.append(self.getEdge(node1, node2))
			cliques.append(out)
		return cliques

	def getBridges(self):
		bridges = list()
		for bridge in self.bridges:
			bridges.append(self.getEdge(bridge[0], bridge[1]))
		return bridges

	def getArticulationPoints(self):
		articulation_points = list()
		for point in self.articulation_points:
			articulation_points.append(point)
		return articulation_points


def extractOrganicFeatures(networkFP, outputFile):
	'''
	extract organic topological features from a network:
	(1) cliques
	(2) cycles (lengths 3, 4, 5, 6)
	(3) bridges
	(4) articulation points
	(5) the remaining edges - "simple path fragments"
	OUTPUT:
	(1), (2) and (5): extract both the residue members in the unit, and distances between the members (for all the edges)
	i.e. for cliques, a list of length (n(V) \choose 2) pairwise distances, and for cycles, a list of p
	airwise distance of lenth = cycle length
	(3): list of edges and the distances associated with the getEdges
	(4): list of vertices
    :param networkFP: filepath pointing to the network file
	:param outputFile: filepath pointing to the output
	'''
	OF = OrganicFeatures()
	if os.path.isfile(networkFP):
		OF.addGraph(networkFP)
		OF.findCycles()
		OF.findCliques()
		OF.findBridges()
		OF.findArticulationPoints()
		cycles = OF.getCycles()
		cliques = OF.getCliques()
		bridges = OF.getBridges()
		articulationPoints = OF.getArticulationPoints()
		simple_path_fragment = OF.getSimplePathFragment()

		with open(outputFile, 'w') as outstream:
			outstream.write('Network: ' + networkFP + '\n\n')
			outstream.write('===== Simple path fragments =====\n')
			outstream.write('res1\tres2\tres1_fullid\tres2_fullid\tdistance\n')
			for edge in simple_path_fragment:
				outstream.write(OF.getNode(edge[0]) + '\t' + OF.getNode(edge[1]) \
					+ '\t' + edge[0] + '\t' + edge[1] + '\t' + edge[2] + '\n')
			outstream.write('\n*** END ***')
			outstream.write('\n\n===== Cycles =====\n\n')
			outstream.write('res\tres_fullid\tdist\n')
			for l in range(len(cycles)):
				members = OF.cycles[l][:-1]
				fullid = ';'.join(members)
				members = [OF.getNode(member) for member in members]
				distances = list()
				for edge in cycles[l]:
					distances.append(edge[0] + '-' + edge[1] + ':' + edge[2])
				outstream.write((';').join(members) + '\t' + fullid  + '\t' + \
					';'.join(distances) + '\n')
			outstream.write('\n*** END ***')
			outstream.write('\n\n===== Cliques =====\n\n')
			outstream.write('res\tres_fullid\tdist\n')
			for c in range(len(cliques)):
				members = OF.cliques[l]
				fullid = ';'.join(members)
				members = [OF.getNode(member) for member in members]
				distances = list()
				for edge in clique:
					distances.append(edge[0] + '-' + edge[1] + ':' + edge[2])
				outstream.write((';').join(members) + '\t' + fullid + '\t' + \
					';'.join(distances) + '\n')
			#for bridgees, output not the resType but the full res ID
			outstream.write('\n*** END ***')
			outstream.write('\n\n===== Bridges =====\n\n')
			outstream.write('res1\tres2\tres1_fullid\tres2_fullid\tdistance\n')
			for edge in bridges:
				outstream.write(OF.getNode(edge[0]) + '\t' + OF.getNode(edge[1]) + \
					'\t' + edge[0] + '\t' + edge[1] + '\t' + edge[2] + '\n')
			#for ar. points, output not the resType but the full res ID
			outstream.write('\n*** END ***')
			outstream.write('\n\n===== Articulation points =====\n\n')
			outstream.write('res\tres_fullid\n')
			for vertex in articulationPoints:
				outstream.write(OF.getNode(vertex) + '\t' + vertex + '\n')
			outstream.write('\n*** END ***')

	else:
		print '"networkFP" does not point to an existing file. Exit.'

def census(directory, tag):
	'''
	census of topological OrganicFeatures
	:param directory: the directory where all pifs output will be parsed and counted
	:param tag: a tag indicating the name of the dataset/cohort; to be included in the output filenames
	'''
	# save in a file all possible configurations seen for each type of features
	# collapse the same pattern from different chains into 1

	# list files
	files = [f for f in os.listdir(directory) if re.match('[\w]{4}_TopolFeat.pifs', f)]
	cycles = list()
	cliques = list()
	fragments = list()
	bridges = list()
	articulation_points = list()

	for file in files:
		with open(os.path.join(directory, file), 'r') as instream:
			pifsOutput = instream.read()
			# start with fragments
			fragment_pat = re.compile('===== Simple path fragments =====(.*?)\*\*\* END \*\*\*', re.M|re.S)
			for match in fragment_pat.finditer(pifsOutput):
				fragment_set = set()
				for line in match.group(0).split('\n'):
					info = line.split('\t')
					if len(info) == 5 and info[0] != 'res1' and info[0] != ' CA' and info[1] != ' CA':
						pattern = (info[2], info[3])
						pattern = (re.sub("-[\w]?:", "-X:", x) for x in pattern)
						if pattern not in fragment_set:
							fragment_set.update(pattern)
							fragment = (info[0], info[1])
							fragment = sorted((re.sub("MSE", "MET", f) for f in fragment))
							fragments.append(fragment[0] + '-' + fragment[1])
			# cycles
			cycle_pat = re.compile('===== Cycles =====(.*?)\*\*\* END \*\*\*', re.M|re.S)
			for match in cycle_pat.finditer(pifsOutput):
				cycle_set = set()
				for line in match.group(0).split('\n'):
					info = line.split('\t')
					if len(info) == 3 and info[0] != 'res' and ' CA' not in info[0]:
						pattern = tuple(info[1].split(';'))
						pattern = (re.sub("-[\w]?:", "-X:", x) for x in pattern)
						if pattern not in cycle_set:
							cycle_set.update(pattern)
							cycle = info[0].split(';')
							cycle = sorted([re.sub("MSE", "MET", l) for l in cycle])
							cycles.append(';'.join(cycle))
			# Cliques
			clique_pat = re.compile('===== Cliques =====(.*?)\*\*\* END \*\*\*', re.M|re.S)
			for match in clique_pat.finditer(pifsOutput):
				clique_set = set()
				for line in match.group(0).split('\n'):
					info = line.split('\t')
					if len(info) == 3 and info[0] != 'res' and ' CA' not in info[0]:
						pattern = tuple(info[1].split(';'))
						pattern = (re.sub("-[\w]?:", "-X:", x) for x in pattern)
						if pattern not in clique_set:
							clique_set.update(pattern)
							clique = info[0].split(';')
							clique = sorted([re.sub('MSE', 'MET', c) for c in clique])
							cliques.append(';'.join(clique))
			# Articulation Points
			articulationPoint_pat = re.compile('===== Articulation points =====(.*?)\*\*\* END \*\*\*', re.M|re.S)
			for match in articulationPoint_pat.finditer(pifsOutput):
				articulationPoint_set = set()
				for line in match.group(0).split('\n'):
					info = line.split('\t')
					if len(info) == 2 and info[0] != 'res' and info[0] != ' CA':
						pattern = info[1]
						pattern = re.sub("-[\w]?:", "-X:", pattern)
						if pattern not in articulationPoint_set:
							articulationPoint_set.update(pattern)
							articulation_points.append(re.sub('MSE', 'MET', info[0]))
			# Bridges
			bridge_pat = re.compile('===== Bridges =====(.*?)\*\*\* END \*\*\*', re.M|re.S)
			for match in bridge_pat.finditer(pifsOutput):
				bridge_set = set()
				for line in match.group(0).split('\n'):
					info = line.split('\t')
					if len(info) == 5 and info[0] != 'res1' and info[0] != ' CA' and info[1] != ' CA':
						pattern = tuple([info[2], info[3]])
						pattern = (re.sub("-[\w]?:", "-X:", x) for x in pattern)
						if pattern not in bridge_set:
							bridge_set.update(pattern)
							bridges.append(re.sub('MSE', 'MET', info[0]) + '-' + re.sub('MSE', 'MET', info[1]))

	# Now do the counting
	cycle_count = collections.Counter(cycles)
	bridge_count = collections.Counter(bridges)
	articulationPoints_count = collections.Counter(articulation_points)
	clique_count = collections.Counter(cliques)
	fragment_count = collections.Counter(fragments)
	with open(os.path.join(directory, tag + '_cycles_census.txt'), 'w') as cycle_wr:
		cycle_wr.write('#PIFS: dataset "' + tag + '": CYCLE CENSUS\n')
		for key, value in cycle_count.items():
			cycle_wr.write(key + '\t' + str(value) + '\n')

	with open(os.path.join(directory, tag + '_bridges_census.txt'), 'w') as bridge_wr:
		bridge_wr.write('#PIFS: dataset "' + tag + '": BRIDGE CENSUS\n')
		for key, value in bridge_count.items():
			bridge_wr.write(key + '\t' + str(value) + '\n')

	with open(os.path.join(directory, tag + '_cliques_census.txt'), 'w') as clique_wr:
		clique_wr.write('#PIFS: dataset "' + tag + '": CLIQUE CENSUS\n')
		for key, value in clique_count.items():
			clique_wr.write(key + '\t' + str(value) + '\n')

	with open(os.path.join(directory, tag + '_articulationPoints_census.txt'), 'w') as articulationPoint_wr:
		articulationPoint_wr.write('#PIFS: dataset "' + tag + '": ARTICULATION POINT CENSUS\n')
		for key, value in articulationPoints_count.items():
			articulationPoint_wr.write(key + '\t' + str(value) + '\n')

	with open(os.path.join(directory, tag + '_fragments_census.txt'), 'w') as fragment_wr:
		fragment_wr.write('#PIFS: dataset "' + tag + '": SIMPLE PATH FRAGMENT CENSUS\n')
		for key, value in fragment_count.items():
			fragment_wr.write(key + '\t' + str(value) + '\n')

	return "CENSUS DONE for " + tag
"""
