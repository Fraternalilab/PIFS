#/usr/bin/python
import chimera
import glob
import subprocess
from chimera import runCommand

def findHbond_chimera(pdb):
    runCommand('open '+ pdb)
    runCommand('select :.b')
    runCommand('findhbond selRestrict cross intramodel true intermodel false saveFile ' + pdb[:-4] + '.hbonds')
    runCommand('close all')
    subprocess.call('grep -v HOH ' + pdb[:-4] + '.hbonds > tmp', shell=True)
    subprocess.call('mv tmp ' + pdb[:-4] + '.hbonds', shell=True)
    return None

if __name__ == "__main__":
    '''
    to be invoked through UCSF Chimera. Open the Chimera program --> Tools --> General Controls --> IDLE. Run this script.
    Remember to set "domain" to the name of the domain to be processed - de-silence the next line and modify accordingly.
    '''
    # domain = sys.argv[1]
    for pdb in ['5keg', '5sww', '5td5', '6bux']:
        tarfile = domain + '_graft' + pdb + '.tar.gz'
        tarfile_restrained = domain + '_graft' + pdb + '_restrained.tar.gz'
        comparisonFile = domain + '_graft' + pdb + '.summary'
        comparisonFile_restrained = domain + '_graft' + pdb + '_restrained.summary'
        # first the non-restrained ones
        subprocess.call('tar -xzf ' + tarfile, shell=True)
        pdbfiles = glob.glob(domain + '_graft' + pdb + '*_dna_*.pdb')
        orig = domain + '_graft' + pdb + '.pdb'
        findHbond_chimera(orig)
        for pdbfile in pdbfiles:
            findHbond_chimera(pdbfile)
        subprocess.call('tar -czf ' + domain + '_graft' + pdb + '_hbonds.tar.gz ' +  domain + '_graft' + pdb + '*.hbonds', shell=True)
        subprocess.call('rm '+ domain + '_graft' + pdb + '*.hbonds', shell=True)
        subprocess.call('rm ' + domain + '_graft' + pdb + '*.pdb', shell=True)
        # then the restrained ones
        print tarfile_restrained
        subprocess.call('tar -xzf ' + tarfile_restrained, shell=True)
        pdbfiles = glob.glob(domain + '_graft' + pdb + '*_dna_*.pdb')
        orig = domain + '_graft' + pdb + '.pdb'
        findHbond_chimera(orig)
        for pdbfile in pdbfiles:
            findHbond_chimera(pdbfile)
        subprocess.call('tar -czf ' + domain + '_graft' + pdb + '_restrained_hbonds.tar.gz ' +  domain + '_graft' + pdb + '*.hbonds', shell=True)
        subprocess.call('rm '+ domain + '_graft' + pdb + '*.hbonds', shell=True)
        subprocess.call('rm ' + domain + '_graft' + pdb + '*.pdb', shell=True)

 
