import csv
from pifsDraw import *

#structDict = dict()
#with open('test_beta2', 'r') as instream:
#    read = csv.reader(instream, delimiter="\t")
#    for row in read:
#        print row[0]
#        protein = row[2].split(',')
#        other = row[3].split(',')
#        structDict[row[1]] = {"file": row[0], "protein": protein, "range": other}
#
for pdb in ['5sww', '5keg', '5td5', '6bux']:
    print pdb
    loadNetworkOnPyMOL('A3_xtal_6A/' + pdb + '_network.txt', pdb + '.pdb', \
        'A3_xtal_6A/' + pdb + '.pse', cutoff=6.0)

#for domain, entry in structDict.iteritems():
#    print domain
#    loadNetworkOnPyMOL('A3_beta2_6A/A3_beta2_6A_' + domain + '_network.txt', \
#        entry['file'], 'A3_beta2_6A/A3_beta2_6A_' + domain + '.pse', cutoff=10.0)
