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

# description: runs yass on the 2 fasta files with the Evalue that was provided

# output: 2 Yass output files (.yass and .yass1)

# Note: the script reverses the first fasta file if necessary and keeps the reversed file:
# 	- if the coverage of the reversed fragments with respect to the yass output is larger than the coverage of the direct fragments, then reverse the first fasta file. 
# 	- the reversed cases are output in the file Inversed_genomes.txt, and the reversed fasta files are marked with a "_".
# 	- if you already know which cases have to be reversed just put your Inversed_genomes.txt file.


if ( $#argv < 3 ) then
    echo "Too few parameters..."
    echo "Usage: " ${1} " fasta file 1"
    echo "	" ${2} "fasta file 2"
    echo "	" ${3} "E-value"
    echo "	" ${4} "reversed targets file (optional)"
    echo "         Yass output: strain1_strain2_evalue.yass (with yass -d 2),  strain1_strain2_evalue.yass1 (with yass -d 1)"
    exit 0
endif

set yass = "./yass_1.14"

set dirseq_yass = "yass_output"
if ( !(-d $dirseq_yass) ) then
	mkdir $dirseq_yass
endif

set yass_blast2chainer_mp = "./yass_blast2chainer 0 2"

set bacteria1 = $1
set bacteria1_p = ${bacteria1:t}         # filename without path
set bacteria1_pe = ${bacteria1_p:r}      # filename without path nor extension

set bacteria2 = $2
set bacteria2_p = ${bacteria2:t}         # filename without path
set bacteria2_pe = ${bacteria2_p:r}      # filename without path nor extension

set target = $bacteria1
set query = $bacteria2


foreach evalue ($3)

	set output_yass = ${bacteria1_pe}"_"${bacteria2_pe}"_"$evalue

	if ( -f ${dirseq_yass}"/"${output_yass}".yass" ) then 
		rm ${dirseq_yass}"/"${output_yass}".yass" 
		rm ${dirseq_yass}"/"${output_yass}".yass1" 
		rm ${dirseq_yass}"/"${output_yass}".mat" 
	endif	
	
	if ( $4 != "" ) then
		if ( `grep ${bacteria1_pe}"_"${bacteria2_pe} $4` != ""   && -f ${target}"_" )  then
			${yass} ${target}"_" ${query} -C 5,-4,-2,-4,-4,5,-4,-2,-2,-4,5,-4,-4,-2,-4,5 -p "#@_##_##_#__@_###,#_##@___##___#___#@#_#" -E $evalue -O 1000000 -d 2 -o ${dirseq_yass}"/"${output_yass}".yass"	
			${yass} ${target}"_" ${query} -C 5,-4,-2,-4,-4,5,-4,-2,-2,-4,5,-4,-4,-2,-4,5 -p "#@_##_##_#__@_###,#_##@___##___#___#@#_#" -E $evalue -O 1000000 -d 1 -o ${dirseq_yass}"/"${output_yass}".yass1"
		else
			if ( `grep ${bacteria1_pe}"_"${bacteria2_pe} $4` != "" ) then
				awk -v fastaFile=${target} -f inverseGenome.awk ${target}
				${yass} ${target}"_" ${query} -C 5,-4,-2,-4,-4,5,-4,-2,-2,-4,5,-4,-4,-2,-4,5 -p "#@_##_##_#__@_###,#_##@___##___#___#@#_#" -E $evalue -O 1000000 -d 2 -o ${dirseq_yass}"/"${output_yass}".yass"	
				${yass} ${target}"_" ${query} -C 5,-4,-2,-4,-4,5,-4,-2,-2,-4,5,-4,-4,-2,-4,5 -p "#@_##_##_#__@_###,#_##@___##___#___#@#_#" -E $evalue -O 1000000 -d 1 -o ${dirseq_yass}"/"${output_yass}".yass1"
			else
				${yass} ${target} ${query} -C 5,-4,-2,-4,-4,5,-4,-2,-2,-4,5,-4,-4,-2,-4,5 -p "#@_##_##_#__@_###,#_##@___##___#___#@#_#" -E $evalue -O 1000000 -d 2 -o ${dirseq_yass}"/"${output_yass}".yass"
				${yass} ${target} ${query} -C 5,-4,-2,-4,-4,5,-4,-2,-2,-4,5,-4,-4,-2,-4,5 -p "#@_##_##_#__@_###,#_##@___##___#___#@#_#" -E $evalue -O 1000000 -d 1 -o ${dirseq_yass}"/"${output_yass}".yass1"
			endif
		endif
	else
		#test whether most of yass fragments are inversed 
		#if this is the case, inverse one of the genome sequences
		#re-run yass with the new genome sequences
		${yass} ${target} ${query} -C 5,-4,-2,-4,-4,5,-4,-2,-2,-4,5,-4,-4,-2,-4,5 -p "#@_##_##_#__@_###,#_##@___##___#___#@#_#" -E $evalue -O 1000000 -d 2 -o ${dirseq_yass}"/"${output_yass}".yass"

		awk -v fastaFile=${target} -v inversedGenomes="Inversed_genomes.txt" -f testInversedGenomes.awk ${dirseq_yass}"/"${output_yass}".yass"

		if ( -f ${target}"_" && `grep ${bacteria1_pe}"_"${bacteria2_pe} Inversed_genomes.txt` != "" ) then
			rm ${dirseq_yass}"/"${output_yass}".yass"
			${yass} ${target}"_" ${query} -C 5,-4,-2,-4,-4,5,-4,-2,-2,-4,5,-4,-4,-2,-4,5 -p "#@_##_##_#__@_###,#_##@___##___#___#@#_#" -E $evalue -O 1000000 -d 2 -o ${dirseq_yass}"/"${output_yass}".yass"
			${yass} ${target}"_" ${query} -C 5,-4,-2,-4,-4,5,-4,-2,-2,-4,5,-4,-4,-2,-4,5 -p "#@_##_##_#__@_###,#_##@___##___#___#@#_#" -E $evalue -O 1000000 -d 1 -o ${dirseq_yass}"/"${output_yass}".yass1"
		else
			${yass} ${target} ${query} -C 5,-4,-2,-4,-4,5,-4,-2,-2,-4,5,-4,-4,-2,-4,5 -p "#@_##_##_#__@_###,#_##@___##___#___#@#_#" -E $evalue -O 1000000 -d 1 -o ${dirseq_yass}"/"${output_yass}".yass1"
		endif
	endif

	${yass_blast2chainer_mp} ${dirseq_yass}"/"${output_yass}".yass"	  ${dirseq_yass}"/"${output_yass}".mat" 
	
end
exit 0

