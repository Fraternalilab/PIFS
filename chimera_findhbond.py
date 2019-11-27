#/usr/bin/python
import chimera
import glob
import subprocess
from chimera import runCommand

def findHbond_chimera(pdb):
    fn = pdb.split('/')[1][:-4]
    runCommand('open '+ pdb)
    runCommand('select :.b')
    runCommand('findhbond selRestrict cross intramodel true intermodel false saveFile ' + fn + '.hbonds')
    runCommand('close all')
    subprocess.call('grep -v HOH ' + fn + '.hbonds > tmp', shell=True)
    subprocess.call('mv tmp ' + fn + '.hbonds', shell=True)
    return None

if __name__ == "__main__":
    '''
    to be invoked through UCSF Chimera. Open the Chimera program --> Tools --> General Controls --> IDLE. Run this script.
    Remember to set "domain" to the name of the domain to be processed - de-silence the next line and modify accordingly.
    '''
    # domain = sys.argv[1]
    for pdb in ['5keg', '5sww', '5td5', '6bux']:
        tarfile = 'pdb/' + domain + '_graft' + pdb + '.tar.gz'
        comparisonFile = 'out/' + domain + '_graft' + pdb + '.summary'
        # first the non-restrained ones
        subprocess.call('mkdir ' + domain + '_graft' + pdb, shell=True )
        subprocess.call('tar -xzf ' + tarfile + ' -C ' + domain + '_graft' + pdb, shell=True)
        pdbfiles = glob.glob(domain + '_graft' + pdb + '/' + domain + '_graft' + pdb + '*_dna_*.pdb')
        orig = domain + '_graft' + pdb + '/' + domain + '_graft' + pdb + '.pdb'
        findHbond_chimera(orig)
        for pdbfile in pdbfiles:
            findHbond_chimera(pdbfile)
        subprocess.call('tar -czf ' + domain + '_graft' + pdb + '_hbonds.tar.gz ' +  domain + '_graft' + pdb + '*.hbonds', shell=True)
        subprocess.call('rm -rf '+ domain + '_graft' + pdb + '/', shell=True)
