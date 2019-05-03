"""
utilities to curate list of proteins with structures of both on its own and in
complex with other proteins, and
extract and compare interface and exposed residues calculated using POPSCOMP
"""

import csv
import os.path
import re
import requests


URL = 'https://www.ebi.ac.uk/pdbe/api/'

def getChainTypes(pdb, type1 = 'polypeptide(L)', type2 = 'polydeoxyribonucleotide'):
	"""
	* USE PDBe API
	give chains of type specified by type1 and type2 for the  supplied pdb code
	:param pdb: list of 4-letter pdb code.
	:param type1: chain type 1
	:param type2: chain type 2
	return a dict of <pdb>: {type1: [], type2: []} 
	"""
	if isinstance(pdb, list) is True:
		postdata = requests.post(URL + 'pdb/entry/molecules/', data = ','.join(pdb), verify = False)
		if postdata.status_code == 200: # test if API is working
			output = postdata.json()
			out = dict()
			for pdbcode in output.keys():
				type1chain = list()
				type2chain = list()
				for entry in output[pdbcode]:
					if type1 in entry['molecule_type']:
						for chain in entry['in_chains']:
							type1chain.append(chain)
					if type2 in entry['molecule_type']:
						for chain in entry['in_chains']:
							type2chain.append(chain)
				out[pdbcode] = {type1: ','.join(type1chain), type2: ','.join(type2chain)}
			return out
		else:
			return 'Error: from HTTP requests - %s' % postdata.reason
	else:
		return 'Error: "pdb" is not a list. Exit.'

if __name__ == '__main__':
	type1='polypeptide(L)'
	type2='polyribonucleotide'
	pdbcodes = list()
	with open('protein-rna-complexes_20181106.csv', 'r') as instream:
		read = csv.reader(instream)
		for row in read:
			if row[0] != 'pdb_id':
				pdbcodes.append(row[0])

	chunk_size = 300
	nested_list = [pdbcodes[i:i + chunk_size] for i in xrange(0, len(pdbcodes), chunk_size)] # nested list of size 300
	for l in nested_list:
		print l
		output = getChainTypes(l, type1, type2)
		with open('protein-rna-complexes_annotated.txt', 'a') as outstream:
			wr = csv.writer(outstream, delimiter="\t")
			for key, value in output.items():
				print key
				wr.writerow([key, value[type1], value[type2]])


