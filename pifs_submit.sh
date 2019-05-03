#!/bin/bash
#£type=rna
#sbatch --nodelist=beo-06 --export=file=`echo $type`,networkDir=`echo $type`_network/,pifsDir=`echo $type`_pifs/,pifsTag=`echo $type` --mem-per-cpu=4G pifs.sh

type=dna
for id in {1..2}; do
	sbatch --nodelist=beo-09 --export=file=`echo $type`_`echo $id`,networkDir=`echo $type`_network/,pifsDir=`echo $type`_pifs/,pifsTag=`echo $type` --mem-per-cpu=4G pifs.sh
done
