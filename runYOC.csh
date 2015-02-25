#!/bin/tcsh
## This file is part of YOC.
## 
## YOC is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## YOC is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with YOC.  If not, see <http://www.gnu.org/licenses/>

# author:Raluca Uricaru

# description: runs YOC with the two scripts: runYass.csh and runOC.csh

# input: <fasta1> <fasta2> <evalue> <overlap_ratio>

# output:

# 1. runYass.csh
# The result files can be found in <yass_output> directory:
#	- 2 Yass output files (.yass, .yass1)
#	- .mat file -- the format needed by the chaining algorithm

# 2. runOC.csh
# The result files can be found in <results> directory
#	- output_YOC_ov<overlap%>.txt containing the coverages for the genome pair
#	- a .chn file containing the chained fragments


if ( $#argv < 2 || $#argv > 4 || $#argv == 3) then
    echo "Usage: <1st_genome.fasta> <2nd_genome.fasta> <evalue> <overlap_ratio>"
    echo "Output: (1) YASS (.mat, .yass), (2) OverlapChainer (.chn), and (3) coverage information (output_YOC_ov<overlap>.txt)"
    exit 0
endif

set fasta1 = $argv[1] #first genome.fasta
set fasta2 = $argv[2] #second genome.fasta
	
if ( $#argv == 2 ) then
	set evalue = 10
	set overlap_ratio = 10
else
	set evalue = $argv[3] #evalue for YASS, 10 by default
	set overlap_ratio = $argv[4] #overlap ratio for OverlapChainer, 10 by default
endif

echo "Running YOC.v1 on " $fasta1 "and " $fasta2 " with Evalue=" $evalue " and overlap=" $overlap_ratio " ..."

# run YASS with the Evalue that was provided
echo "1. run YASS on the 2 genomes with Evalue =" $evalue
./runYass.csh $fasta1 $fasta2 $evalue
echo "Finished running YASS"

# run OC with overlap ratio 10 on the .mat file

set dirseq_yass = "yass_output"

set fasta1_p = ${fasta1:t}         # filename without path
set fasta1_pe = ${fasta1_p:r}      # filename without path nor extension

set fasta2_p = ${fasta2:t}         # filename without path
set fasta2_pe = ${fasta2_p:r}      # filename without path nor extension

echo "2. run OverlapChainer on YASS result and obtaining chaining of fragments in Chainer format"
./runOC.csh $fasta1_pe $fasta2_pe $dirseq_yass/$fasta1_pe"_"$fasta2_pe"_"$evalue".mat" $overlap_ratio
echo "Finished running OverlapChainer"
echo "Results were generated in yass_output (for YASS) and results folder (YOC chaining results and coverage statistics)"

exit 0;
