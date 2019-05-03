"""
extract topological features from networks
"""

import argparse
import os.path
import pifsUtils
import re

if __name__=="__main__":
	'''
	start the extraction
	'''
	parser = argparse.ArgumentParser(description="This script extract topological features from residue networks")
	parser.add_argument('--inDir', help = 'input directory with network sif files', type = str, required = True)
	parser.add_argument('--outDir', help = 'output directory of list of topological features', type = str, required = True)
	args = parser.parse_args()

	files = [f for f in os.listdir(args.inDir) if re.match('[\w]{4}_network.txt', f)]
	for file in files:
		pdbcode = file[:4]
		print pdbcode + ' ...'
		pifsUtils.extractOrganicFeatures(os.path.join(args.inDir, file), \
			os.path.join(args.outDir, pdbcode + '_TopolFeat.pifs'))
