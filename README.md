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

1. Encode residue subgraphs for protein-X interface (where X can be protein/DNA/RNA)

```
python pifs.py
```

This will generate, for each input structure, one residue network centred around the interface between two specified chains. [POPSCOMP](https://github.com/Fraternalilab/POPSCOMPlegacy) and [POPS](https://github.com/Fraternalilab/POPSlegacy) is called to identify such interface. Their locations could be specified by options `--popscompDir` and `--popsDir`. 

This code requires as an input the PDB structures (the directory could be specified by `--inDir`), and a list indicating which chains are to be considered. An example is given in the file `pifs_testInput`, of the following format:

```
# tab-delimited
# name_of_structures  chain1  chain2
myFavStructure1  A B
myFavStructure2  A C
```

The commented lines are NOT to be included in the actual file. Each line corresponds to one input file. For example, the first line would entail analysing the interface between chains A and B in the pdb file `myFavStructure1.pdb`.

A detailed description of options could be given by `python pifs.py --help`.


2. Encode residue subgraphs for a given region of the protein 




