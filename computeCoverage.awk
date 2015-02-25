#!/usr/bin/awk -f
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

# description: computeCoverage.awk
# input: 
#    chainer output .chn

# output:
# 1. coverage
#    (total length of fragments taken in the chain, without overlaps)
# 2. overlaps
# 3. fragment lengths distribution
# 4. total fragment length per length group



BEGIN{
#  FS = "[ ,\[ \]]*";
  Couverture[1] = Couverture[2] = 0;
  Gaps[1] = Gaps[2] = 0;
  score = 0;

  
  if ( type == 0 || type == 3)
  {
	#vmatch chainer type = 0
	#vmatch OC type = 3
	no_fragments_300 = 0;
	no_fragments_500 = 0;
	no_fragments_sup500 = 0;

	length1_fragments_300 = 0;
	length1_fragments_500 = 0;
	length1_fragments_sup500 = 0;

	length2_fragments_300 = 0;
	length2_fragments_500 = 0;
	length2_fragments_sup500 = 0;
  }
  else
  if ( type == 1  || type == 2 )
  {
	#yass/blast chainer type = 1
	#yass/blast OC type = 2
	no_fragments_1kb = 0;
	no_fragments_5kb = 0;
	no_fragments_10kb = 0;
	no_fragments_sup10kb = 0;

	length1_fragments_1kb = 0;
	length1_fragments_5kb = 0;
	length1_fragments_10kb = 0;
	length1_fragments_sup10kb = 0;

	length2_fragments_1kb = 0;
	length2_fragments_5kb = 0;
	length2_fragments_10kb = 0;
	length2_fragments_sup10kb = 0;
  }

  if ( type == 2 || type == 3)
  {
	weight_sum = 0;
	Overlap[1] = 0;
	Overlap[2] = 0;
  }
  
  #first type of insertion: 
  #gap = 0
  #gap < 20 
  #gap < 0
  #on one of the sequences
  Insert1[1] = 0;
  Insert1[2] = 0;
  no1 = 0;
  
  #second type of insertion: 
  #min_gap < 10% of the 2 fragments which flank the gap
  Insert2[1] = 0;
  Insert2[2] = 0;
  no2 = 0;

  #third type of insertion:
  #it takes in accoumpt the context of the gaps
  Insert3[1] = 0;
  Insert3[2] = 0;
  no3 = 0;

  prec1 = 0;
  prec2 = 0;
  prec3 = 0;
  prec4 = 0;

  #minimum length of an insertion
  MIN_INSERT = 50;
}

# score
(NR == 2){
  split($1, a, "#");
  score = a[2];
  #print a[2];
}


# lignes d'entete: imprimer
(NR > 2){


  if ($1 == "total"){
	time = $3;
  }
  else
  {
	  weight_sum += $3;

	  gsub(/[,\[\]]/, " ");

	  if (debug) {
	    print NF;
	    print ;
	  }

	  if ( $2 == $1 )
	  {
		LgCour[1] = 0;
		LgCour[2] = 0;
	  }
	  else{
	  	LgCour[1] = $2 - $1 +1;
	  	LgCour[2] = $4 - $3 +1;
	  }
	  
	  gap1 = 0;
	  gap2 = 0;

	  if (NR == 3)
	  {
	    if ( type != 2 && type != 3)
	    {   	   
		Gaps[1] = $1;	   
		Gaps[2] = $3;	    
	    }
	    else
	    {
	     	Gaps[1] = 0;
		Gaps[2] = 0; 
	    }

	    prec1 = $1;
	    prec2 = $2;	
	    prec3 = $3;
	    prec4 = $4;	
	  }
	  else
	  {
	    if ( type != 2 && type != 3)
	    {	     
		if ( $1 != $2 ){
			Gaps[1] += $1 -$2 -1;
			Gaps[2] += $3 - $4 -1;	
		}
		gap1 = prec1 - $2 -1;
		gap2 = prec3 - $4 -1;    
	    }
	    else
	    {
		if ( $2 >= prec1 )
			Overlap[1] += $2 - prec1 +1;
		else
			Gaps[1] += prec1 - $2 -1;

		gap1 = prec1 - $2;

		if ( $4 >= prec3 )
			Overlap[2] += $4 - prec3 +1;
		else
			Gaps[2] += prec3 - $4 -1;    

		gap2 = prec3 - $4;
	    }
	  }
	  
	  if ( LgCour[1] < 0 ) LgCour[1] = (-1) * LgCour[1];
	  if ( LgCour[2] < 0 ) LgCour[2] = (-1) * LgCour[2];

	  
	  Couverture[1] += LgCour[1];
	  Couverture[2] += LgCour[2];
	  
	  #analyse the gaps on the 2 sequences in order to decide whether we deal with an insertion or not
	  if ( gap1>=0 && gap2>=0 )
	  {
		if ( gap1 == 0 && gap2 >= MIN_INSERT )
		{
			#insert = gap2 on seq2
			Insert1[2] += gap2;
			no1++;
			#print "insert2 ", gap1, "(", lg1_prec, LgCour[1], ")", gap2, "(", lg2_prec, LgCour[2], ")" ;
			#print  $2, prec1, $4, prec3 ;
		
		}
		else
		if ( gap2 == 0 && gap1 >= MIN_INSERT )
		{
			#insert = gap1 on seq1
			Insert1[1] += gap1;
			no1++;
			#print "insert1 ", gap1, "(", lg1_prec, LgCour[1], ")", gap2, "(", lg2_prec, LgCour[2], ")" ;
			#print  $2, prec1, $4, prec3 ;
		}
		else
		{
			#the 2 gaps > 0
			#we choose the smallest gap
			#we verify whether it is less than 10% from the smallest fragment 
			#(between the 2 fragments that flank our gap)
			lg1_prec = prec2 - prec1;
			lg2_prec = prec4 - prec3;
		
		
			(gap1 > gap2) ? min_gap = gap2 : min_gap = gap1;
			if ( min_gap == gap1 )
				max_gap = gap2;
			else
				max_gap = gap1;
		
			#if our smallest gap is less than half of the biggest
			if ( max_gap - min_gap >= MIN_INSERT)
			{
			#if our smallest gap is less than 20
			if ( min_gap <= 20 )
			{
				if ( min_gap == gap1 )
				{
					if ( gap1 <= lg1_prec/10 && gap1 <= LgCour[1]/10 )
					{
						Insert1[2] += gap2 - gap1;
						#print "insert2 ", gap1, "(", lg1_prec, LgCour[1], ")", gap2, "(", lg2_prec, LgCour[2], ")" ;
						#print  $2, prec1, $4, prec3;
						no1++;
					}
				}
				else
				{
					if ( gap2 <= lg2_prec/10 && gap2 <= LgCour[2]/10 )
					{
						Insert1[1] += gap1 - gap2;
						#print "insert1 ", gap1, "(", lg1_prec, LgCour[1], ")", gap2, "(", lg2_prec, LgCour[2], ")" ;
						#print  $2, prec1, $4, prec3;
						no1++;
					}
				}
			}
			else
			{
				if ( min_gap == gap1 )
				{
					if ( gap1 <= lg1_prec/10 && gap1 <= LgCour[1]/10 )
					{
						Insert2[2] += gap2 - gap1;
						#print "2: insert2 ", gap1, "(", lg1_prec, LgCour[1], ")", gap2, "(", lg2_prec, LgCour[2], ")" ;
						#print  $2, prec1, $4, prec3;
						no2++;
					}
				}
				else
				{
					if ( gap2 <= lg2_prec/10 && gap2 <= LgCour[2]/10 )
					{
						Insert2[1] += gap1 - gap2;
						#print "2: insert1 ", gap1, "(", lg1_prec, LgCour[1], ")", gap2, "(", lg2_prec, LgCour[2], ")" ;
						#print  $2, prec1, $4, prec3;
						no2++;
					}
				}
			}}
		}
	  }
	  else
	  if ( gap1 < 0 && gap2 > 0 && gap2 >= MIN_INSERT )
	  {
	  	Insert1[2] += gap2;
		#print "insert2 ", gap1, "(", lg1_prec, LgCour[1], ")", gap2, "(", lg2_prec, LgCour[2], ")" ;
		#print  $2, prec1, $4, prec3 ;
		no1++;
	  }
	  else
	  if ( gap2 < 0 && gap1 > 0 && gap1 >= MIN_INSERT)
	  {
		Insert1[1] += gap1;
		#print "insert1 ", gap1, "(", lg1_prec, LgCour[1], ")", gap2, "(", lg2_prec, LgCour[2], ")" ;
		#print  $2, prec1, $4, prec3 ;
		no1++;
	  }
		  

	  if ( LgCour[1] != 0 && LgCour[1] != 0 )
	  {
	 
	   # we keep every fragment that we find
	   #first fragment = pos 1
	  fragment_x1[NR-3] = $1;
	  fragment_x2[NR-3] = $3;
	  fragment_y1[NR-3] = $2;
	  fragment_y2[NR-3] = $4;

	  if ( type == 0 || type == 3 )
		{
		#vmatch
			if ( LgCour[1] <= 300 && LgCour[2] <= 300 )
			{
				no_fragments_300++;
				length1_fragments_300 += LgCour[1];
				length2_fragments_300 += LgCour[2];
			}
			else
			if ( (LgCour[1] > 300 && LgCour[1] <= 500 && LgCour[2] <= 500) || 
	 		     (LgCour[2] > 300 && LgCour[2] <= 500 && LgCour[1] <= 500) )
			{
				no_fragments_500++;
				length1_fragments_500 += LgCour[1];
				length2_fragments_500 += LgCour[2];
			}
			else
			{
				no_fragments_sup500++;
				length1_fragments_sup500 += LgCour[1];
				length2_fragments_sup500 += LgCour[2];
			}
		}
		else
		{
		#yass
			if ( LgCour[1] <= 1000 && LgCour[2] <= 1000 )
			{
				no_fragments_1kb++;
				length1_fragments_1kb += LgCour[1];
				length2_fragments_1kb += LgCour[2];
			}
			else
			if ( (LgCour[1] > 1000 && LgCour[1] <= 5000 && LgCour[2] <= 5000) || 
	 		     (LgCour[2] > 1000 && LgCour[2] <= 5000 && LgCour[1] <= 5000) )
			{
				no_fragments_5kb++;
				length1_fragments_5kb += LgCour[1];
				length2_fragments_5kb += LgCour[2];
			}
			else
			if ( (LgCour[1] > 5000 && LgCour[1] <= 10000 && LgCour[2] <= 10000) || 
	 		     (LgCour[2] > 5000 && LgCour[2] <= 10000 && LgCour[1] <= 10000) )
			{
				no_fragments_10kb++;
				length1_fragments_10kb += LgCour[1];
				length2_fragments_10kb += LgCour[2];
			}
			else
			{
				no_fragments_sup10kb++;
				length1_fragments_sup10kb += LgCour[1];
				length2_fragments_sup10kb += LgCour[2];
			}
		}
	  }

	  prec1 = $1;
	  prec2 = $2;
	  prec3 = $3;
	  prec4 = $4;
  }
  
}
END{
  if ( type == 2 || type == 3 )
  {
	Couverture[1] -= Overlap[1];
	Couverture[2] -= Overlap[2];
  }

  #we take context into account
  #NO_F fragments at left and NO_F at right
  NO_F = 5
  Insert3[1] = 0;
  Insert3[2] = 0;

  for (i=1; i<NR-4; i++)
  {
	if ( i<=NO_F )
		min = 1;
	else
		min = i - NO_F + 1;
	if ( i + NO_F <= NR-4 )
		max = i+NO_F;
	else
		max = NR-4;

	coverage_before_1 = 0;
	coverage_after_1 = 0;
	coverage_before_2 = 0;
	coverage_after_2 = 0;

	len_before_1 = fragment_y1[min] - fragment_x1[i];
	len_before_2 = fragment_y2[min] - fragment_x2[i];
	
	coverage_before_1 += fragment_y1[i] - fragment_x1[i1];
	coverage_before_2 += fragment_y2[i] - fragment_x2[i];

	for ( j=min; j<i; j++ )
	{
		x1 = fragment_x1[j];
		y1 = fragment_y1[j+1];
		
		coverage_before_1 += fragment_y1[j] - fragment_x1[j];

		if ( x1 < y1 )
			coverage_before_1 += x1 -y1;
			
		x2 = fragment_x2[j];
		y2 = fragment_y2[j+1];

		coverage_before_2 += fragment_y2[j] - fragment_x2[j];
		
		if ( x2 < y2 )
			coverage_before_2 += x2 - y2;
	}
 
	len_after_1 = fragment_y1[i+1] - fragment_x1[max];
	len_after_2 = fragment_y2[i+1] - fragment_x2[max];

	coverage_after_1 += fragment_y1[max] - fragment_x1[max];
	coverage_after_2 += fragment_y2[max] - fragment_x2[max];

	for ( j=i+1; j<max; j++ )
	{
		x1 = fragment_x1[j];
		y1 = fragment_y1[j+1];
		
		coverage_after_1 += fragment_y1[j] - fragment_x1[j];
		
		if ( x1 < y1 )
			coverage_after_1 += x1 -y1;
			
		x2 = fragment_x2[j];
		y2 = fragment_y2[j+1];

		coverage_after_2 += fragment_y2[j] - fragment_x2[j];
		
		if ( x2 < y2 )
			coverage_after_2 += x2-y2;	
	}

	percentage = 7/10;
	cov_percentage = 20;
	if (coverage_before_1 >= percentage * len_before_1 &&
	    coverage_before_2 >= percentage * len_before_2 &&
	    coverage_after_1 >= percentage * len_after_1 &&
	    coverage_after_2 >= percentage * len_after_2)
	{
		gap1 = fragment_x1[i] - fragment_y1[i+1];
		gap2 = fragment_x2[i] - fragment_y2[i+1];
		
		(gap1 > gap2) ? min_gap = gap2 : min_gap = gap1;
		if ( min_gap == gap1 )
			max_gap = gap2;
		else
			max_gap = gap1;
		
		#if our smallest gap is less than half of the biggest
		if ( max_gap-min_gap >= MIN_INSERT)
		{
		if ( min_gap > 20 )
		{
			if ( min_gap == gap1 )
			{
				if ( gap1 <= coverage_before_1/cov_percentage && gap1 <= coverage_after_1/cov_percentage )
				{
					Insert3[2] += gap2 - gap1;
					#print "3: insert2 ", gap1, "(", coverage_before_1, coverage_after_1, ")", gap2, "(", coverage_before_2, coverage_after_2, ")" ;
					#print  fragment_y1[i+1] , fragment_x1[i], fragment_y2[i+1], fragment_x2[i];
					no3++;
				}
			}
			else
			{
				if ( gap2 <= coverage_before_2/cov_percentage && gap2 <= coverage_after_2/cov_percentage )
				{
					Insert3[1] += gap1 - gap2;
					#print "3: insert1 ", gap1, "(", coverage_before_1, coverage_after_1, ")", gap2, "(", coverage_before_2, coverage_after_2, ")" ;
					#print  fragment_y1[i+1] , fragment_x1[i], fragment_y2[i+1], fragment_x2[i];
					no3++;
				}
			}
		}}
	}
  }

  print "Coverages ", Couverture[1], Couverture[2];
  #print "Insert1 ", Insert1[1], Insert1[2], no1;
  #print "Insert2 ", Insert2[1], Insert2[2], no2;
  #print "Insert3 ", Insert3[1], Insert3[2], no3;
  #print "Gaps ", Gaps[1], Gaps[2], "Total_gaps ", Gaps[1]+Gaps[2];
  if ( type == 2 || type == 3 ) 
	print "Overlaps ", Overlap[1], Overlap[2];
  if ( type != 2 && type != 3 ) 
  {
	#print "Score ", score;
  	# types de chainage: g, gg
  	#print "GG Weight_sum ", score + Gaps[1] + Gaps[2] - 2;
  	# print "G Weight_sum ", score - 2;
  	#print "CHAINER output";
  }
  else
  { 
	#print "Weight_sum ", weight_sum;
	#print ”Chaining output";
  }
  print "Fragments ", NR-4;
  if ( type == 0 || type == 3 )
  {
	print "300 500";
	print no_fragments_300, " ", no_fragments_500, " ", no_fragments_sup500;
	print length1_fragments_300, " ", length1_fragments_500, " ", length1_fragments_sup500;
	print length2_fragments_300, " ", length2_fragments_500, " ", length2_fragments_sup500;
  }
  else
  {
	print "1000 5000 10000";
	print no_fragments_1kb, " ", no_fragments_5kb, " ", no_fragments_10kb, " ", no_fragments_sup10kb;
	print length1_fragments_1kb, " ", length1_fragments_5kb, " ", length1_fragments_10kb, " ", length1_fragments_sup10kb;
	print length2_fragments_1kb, " ", length2_fragments_5kb, " ", length2_fragments_10kb, " ", length2_fragments_sup10kb;
  }
  print "Time ", time;
  
}
