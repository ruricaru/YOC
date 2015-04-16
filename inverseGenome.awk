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

# description: innverseGenome.awk
# input: fasta file to reverse and complement

END{

system("head -n 1 " fastaFile " > " fastaFile"_");

system("  awk '{ a[NR]=$0; } END { for(i=NR; i>1; --i) { for (j=length(a[i]); j; --j) { ch = substr(a[i], j, 1); if (ch==\"t\" || ch==\"T\") {ch=\"A\";} else if (ch==\"a\" || ch==\"A\") {ch=\"T\";}  else if (ch==\"g\" || ch==\"G\") {ch=\"C\";}  else if (ch==\"c\" || ch==\"C\") {ch=\"G\";} else ch = toupper(ch); printf(\"%c\",ch) }  } } ' " fastaFile " >> " fastaFile"_");
}

