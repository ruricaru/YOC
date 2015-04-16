#!/usr/bin/env awk -f
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
# creation date: 24.09.2007

# description: genome_length.awk

#input: 
#    fasta genome sequence file

# output:
#    genome length


BEGIN{
  genome_length = 0
}

(NR > 1){
  genome_length += length($1);
}

END{
  print genome_length;
}
