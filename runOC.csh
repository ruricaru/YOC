#!/usr/bin/env tcsh

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

# description: for a .mat file we run OverlapChainer with an overlap ratio given as a parameter

# output: a .chn file and a file containing the coverage results

if ( $#argv < 3 ) then
    echo "Usage: " ${1} " strain1 name"
    echo "	" ${2} "strain2 name"
    echo "	" ${3} "mat file"
    echo "	" ${4} "overlap ratio (optional, default 10)"
    echo "         OC output .chn and a file with the coverage information"
    exit 0
endif

set bact1 = $1
set bact2 = $2
set mat_file = $3
set mat_file_p = ${mat_file:t}         # filename without path
set mat_file_pe = ${mat_file_p:r}      # filename without path nor extension

set overlap_ratio = $4
if ( $4 == "" ) then
	set overlap_ratio = 10
endif

set comp_cov = "computeCoverage.awk"
#set comp_cov_sameas_filtered = "comp_Chainer+OC_coverage_sameas_filtered.awk"

set results_dir = "results/"
if ( !(-d $results_dir)) then
	mkdir $results_dir
endif

set results_YOC = $results_dir"YOC_ov"$overlap_ratio"/"
if ( !(-d $results_YOC) ) then
	mkdir $results_YOC
endif

# The OC program may decide that the .mat file has a problem of circularity and thus the output may be shifted. You know if this is the case by looking in the .chn file for >H 2 Circular followed by the left shift and the right shift. 
# The coordinates of the fragments in this case are shifted, based on the coordinates in the second fasta file. 
# In order to obtain the true coordinates on the second fasta file, we substract the right shift from all x2 and y2 from fragments with x2 >= right shift (where x2 = beginning of the fragment on the second fasta file, y2 end of the fragment), and add the left shift to x2 and y2 of all the other fragments. 
# The coordinates on the first fasta file do not change. 

set OC = "java -cp OverlapChainer/bin tests.OverlapChainer "$overlap_ratio
echo $OC

#Keep this part only to be consistent with previous versions BUT DO NOT USE COVERAGES, GAPS LENGTHS, MPs FROM HERE
set output_YOC = $results_dir"output_YOC_ov"$overlap_ratio".txt"
set output_YOC_p = ${output_YOC:t}         # filename without path
set output_YOC_pe = ${output_YOC_p:r}      # filename without path nor extension

set output_chain_YOC = $results_YOC$mat_file_pe".chn"
rm -f $output_chain_YOC

${OC} $mat_file >> $output_chain_YOC
echo $output_chain_YOC".couv" >> $output_YOC
echo $bact1" "$bact2" "$overlap_ratio >> $output_YOC
awk -v type=2 -f $comp_cov $output_chain_YOC >> $output_YOC

