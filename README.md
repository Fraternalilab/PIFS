# PIFS
Protein Interface Fingerprinting with Subgraphs

*Joseph Ng, November 2019*

## Introduction

PIFS is currently a collection of Python scripts developed to fingerprint protein interfaces of interest, using "subgraphs" (i.e. small residue networks, with C-alphas connected based on a user-defined distance threshold). It interfaces with atomic coordinate files (in Protein Data Bank [PDB] format; i.e. `*.pdb`) directly, offers functionalities to encode and visualise subgraphs, evaluation of interfaces, as well as detection of problematic cases such as entanglement of chains. 

## Content

In this repository the codes for the functions, sample scripts to process files, as well as sample input are provided:

| Task | Functions | Sample scripts of heuristics | Sample input | Sample output |
| ---- | --------- | ---------- | ------------ | ------------- |
| Encode residue subgraphs for protein-X interface (where X can be protein/DNA/RNA) | `pifsUtils.py` | `pifs.py` | |
| Encode residue subgraphs for a given region of the protein | `pifsUtils.py` | `pifs_SeqRegion.py` | |
| Visualise residue subgraphs on PyMOL | `pifsDraw.py` | `pifsDrawNetworks.py` | |
| Detect pi-pi stacking contacts between chains | `pifsCalcInteractions.py` | `A3_interchain.py` | |
| Detect hydrogen bonds between chains (interfacing with UCSF Chimera) | `chimera_findhbond.py` | *NA* | |
| Detect problematic interfaces (entanglement of chains) | `pifsEntanglement.py` | `A3_entangle.py` | |
| Visualise problematic interfaces (entanglement of chains) | `pifsDrawEntanglement.py` | *NA (see below)* | |

## Usage

### 1. Encode residue subgraphs for protein-X interface (where X can be protein/DNA/RNA)

```
# example usage; please change --popsDir and --popscompDir as appropriate (see below description)
python pifs.py --structList pifs_testInput --inDir pdb \
  --networkOutDir out --popscompOutDir out \
  --popsDir ~/pops --popscompDir ~/POPSCOMP --errorFile out/out.err
```

This will generate, for each input structure, one residue network centred around the interface between two specified chains. [POPSCOMP](https://github.com/Fraternalilab/POPSCOMPlegacy) and [POPS](https://github.com/Fraternalilab/POPSlegacy) is called to identify such interface. Their locations could be specified by options `--popscompDir` and `--popsDir`. 

This code requires as an input the PDB structures (the directory could be specified by `--inDir`), and a list indicating which chains are to be considered. An example is given in the file `pifs_testInput`, of the following format:

```
# tab-delimited
# name_of_structures  chain1  chain2
FavStructure1	A	B
FavStructure2	A	C
```

The commented lines are NOT to be included in the actual file. Each line corresponds to one input file. For example, the first line would entail analysing the interface between chains A and B in the pdb file `FavStructure1.pdb`.

A detailed description of options could be given by `python pifs.py --help`.


### 2. Encode residue subgraphs for a given region of the protein 

```
# example usage; please change --popsDir and --popscompDir as appropriate (see below description)
python pifs_SeqRegion.py --structList pifsRegion_testInput --inDir pdb \
  --networkOutDir out --errorFile out/out.err
```

Same as (1), but here the region of interest is given in the text-based input file (see `pifsRegion_testInput` for an example), of the following format:
```
# tab-delimited
# name_of_structures	chain	region
FavStructure1	A	20,33
FavStructure2	B	43,66
```

The commented lines are NOT to be included in teh actual file. Each line corresponds to one input file. For example, the first line would entail giving the residue network for the region between residues 20 and 33 of chain A of `FavStructure1.pdb`.

A detailed description of options could be given by `python pifs_SeqRegion.py --help`.

### 3. Visualise residue subgraphs on PyMOL

This could be achieved by a function `loadNetworkOnPyMOL` in `pifsDraw.py` which interfaces with PyMOL and create and PyMOL session file with an input residue graph visualised in 3D. An example usage could be found in `pifsDrawNetworks.py`.

### 4. Detect pi-pi stacking contacts between chains

```
python A3_interchain.py A3A
```
would perform the calculation for the domain A3A.

This is intended to be performed over large *ensembles* of protein-DNA complexes - here compressed in a tar file. The code would loop through every structure in the ensemble and detect pi-pi stacking between the protein and DNA chains. As an input the range of DNA bases to be considered has to be specified. `A3_DNA_Definition` gives a self-explanatory example of how this range could be specified.

### 5. Detect hydrogen bonds between chains (interfacing with UCSF Chimera)

`chimera_findhbond.py` was written to detect hydrogen bonds, again over large ensembles of protein-DNA complexes. This is intended to be run using the Python shell provided on UCSF Chimera, so that users would indicate the name of the domain (or protein; as indicated in the filenames of the pdb files) on the code, and run the module via the Chimera GUI & Python Shell. The output would be one text-based file listing hydrogen bonds between the protein and DNA chains, for each pdb file in the ensemble. Please consult the manual of UCSF Chimera on running custom Python scripts.

### 6. Detect problematic interfaces (entanglement of chains) 

```
python A3_interchain.py A3A
```
would perform the calculation for the domain A3A.

These codes were written to deal with specifically the "entanglement" of nucleic acid chain inside a protein loop, which was observed in certain protein-DNA structural poses obtained via *in silico* modelling. As input, a definition file for the range of DNA bases to consider (as above; see point no. 4) has to be specified, as well as another text-based file indicating the residue range for each loop on the protein chain to be examined. These codes used geometrical methdos to detect automatically such cases of nucleic acid entanglement. See Ch. 5 of the thesis (Ng, 2019) for further information.

### 7. Visualise problematic interfaces (entanglement of chains)

`pifsDrawEntanglement.py` was written as a plugin to PyMOL to visualise a case of entanglement as defined by the user. Ch. 5 of the thesis (Ng, 2019) - especially listing 5.5 and the discussion surrounding it offers a more explicit explanation. The code provides an additional PyMOL function `entangle` which could be invoked on the PyMOL command line - see below for an example:

```
# !/ usr / bin / pymol
run pifsDrawEntanglement.py
# Open a PDB file for visualisation
# Say the structural object is called ’test1’,
# an example of combination to depict :
entangle test1 , p1=22 , p2=28 , p3=35 , n1=1 , n2=2
``` 

This would depict a triangle defined on the protein loop (here with residues 22, 28 and 35) and an arrow defined on the nucleic acid chain (in this case from base position 1 to 2) - to aid visual inspection of individual entanglement cases.
