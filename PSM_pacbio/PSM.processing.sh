#!/bin/bash

mkdir psm

names="Brec1 D7L D7R Cre"

parallel "exonerate -m affine:bestfit -E yes -t fasta/PSM_{1}.fasta -q psm_references/{1}.fa --verbose 0 --showsugar no --showcigar no --showvulgar no --showalignment no --percent 80 --ryo '>%ti\tmissmatch:%em\tlength:%tl\tcigar:%C\n%tas\n' | seqtk seq > psm/{1}.fa" ::: $names
