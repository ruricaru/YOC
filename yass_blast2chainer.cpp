/*
This file is part of YOC.

YOC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

YOC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with YOC.  If not, see <http://www.gnu.org/licenses/>
*/

/* 
 * author:Raluca Uricaru
 *
 * description: extracts information from YASS and blast files and builds chaining input files
 * ./yass_blast2chainer.cpp <type_of_weight> <type_of_fragments> <input_file>
 * input: 
 *  <type_of_weight> = 0 (nb of identities) or 1 (alignment length)
 *  <type_of_fragment> = 1 (for old versions of YASS), 2 (for YASS or blast) 
 *  <input_file> = yass file or blast file
 * output: <input_file_name>.mat file (chainer input file)
 *
 * In a chainer input file fragments are given as follows
 *  #<fragment_weight>
 *  [start_seq1,end_seq1] [start_seq2,end_seq2]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char **argv) 
{
	//the yass file
	FILE *yass_blast;
	int i,j;
	char buffer[512], *s;
	float identity, evalue, bitscore;
	int alignment_length, mismatches, gap, q_start, q_end, s_start, s_end;
	int type_of_weight;

	/*
		1 for wrong yass
		2 for blast or for yass with corrected %id, alignment_length, mismatches and gap openings
	*/
	int type_of_software;
	int alignment_length1, alignment_length2, gaps_long_seq, gaps_commun;

	//the chainer file
	char *chainer_name = (char*)malloc(300*sizeof(char));
	FILE *chainer;

	s = (char*)malloc(300*sizeof(char));

	//reading the input parameters
	if ( argc < 4 || argc > 5 )
	{
		printf("The number of parameters doesn't correspond %d", argc);
		exit(1);
	}

	if ( (yass_blast = fopen(argv[3], "r")) == NULL )	
	{
		printf("Failed to open yass file %s\n", argv[3]);
		exit(1);
	}

	type_of_weight = atoi(argv[1]);

	type_of_software = atoi(argv[2]);

	if ( argc == 4 )
	{
		strncpy(chainer_name, argv[3], strlen(argv[3]) - strlen(strchr(argv[3], '.')));
		chainer_name[strlen(argv[3]) - strlen(strchr(argv[3], '.'))] = '\0';
		strcat(chainer_name, ".mat");
	}
	else
	{
		strcpy(chainer_name, argv[4]);
	}

	//open chainer file
	if ( (chainer = fopen(chainer_name, "w+")) == NULL )	
	{
		printf("Failed to open chainer file %s\n", chainer_name);
		exit(1);
	}

	//sscanf("93.356", "%f", &identity);
	/*identity = atof("93.25");
	printf("%f", identity);*/

	//build the chainer file
	fprintf(chainer, ">CHA 2\n");
	while ( fgets(buffer, 512, yass_blast) != NULL )
	{
		if ( buffer[0] == '#' )
			continue;

		/* # Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score */

		// query id
		s = strtok(buffer, "\t");
		if ( !s )
		{
			//printf("Wrong yass file\n");
			exit(1);
		}

		// subject id
		s = strtok(NULL, "\t");
		if ( !s )
		{
			//printf("Wrong yass file\n");
			exit(1);
		}

		// % identity
		s = strtok(NULL, "\t");
		if ( !s )
		{
			//printf("Wrong yass file\n");
			exit(1);
		}
		sscanf(s++, "%f", &identity);

		// alignment length
		s = strtok(NULL, "\t");
		if ( !s )
		{
			//printf("Wrong yass file\n");
			exit(1);
		}
		sscanf(s++, "%d", &alignment_length);

		// mismatches
		s = strtok(NULL, "\t");
		if ( !s )
		{
			//printf("Wrong yass file\n");
			exit(1);
		}
		sscanf(s++, "%d", &mismatches);

		// gap openings
		s = strtok(NULL, "\t");
		if ( !s )
		{
			//printf("Wrong yass file\n");
			exit(1);
		}
		sscanf(s++, "%d", &gap);

		// q.start
		s = strtok(NULL, "\t");
		if ( !s )
		{
			//printf("Wrong yass file\n");
			exit(1);
		}
		sscanf(s++, "%d", &q_start);

		// q.end
		s = strtok(NULL, "\t");
		if ( !s )
		{
			//printf("Wrong yass file\n");
			exit(1);
		}
		sscanf(s++, "%d", &q_end);

		// s.start
		s = strtok(NULL, "\t");
		if ( !s )
		{
			//printf("Wrong yass file\n");
			exit(1);
		}
		sscanf(s++, "%d", &s_start);

		// s.end
		s = strtok(NULL, "\t");
		if ( !s )
		{
			//printf("Wrong yass file\n");
			exit(1);
		}
		sscanf(s++, "%d", &s_end);

		// e-value
		s = strtok(NULL, "\t");
		if ( !s )
		{
			//printf("Wrong yass file\n");
			exit(1);
		}
		sscanf(s++, "%g", &evalue);

		// bit score
		s = strtok(NULL, "\t");
		if ( !s )
		{
			//printf("Wrong yass file\n");
			exit(1);
		}
		sscanf(s++, "%g", &bitscore);
		
		//fprintf(chainer, "#%f\n", identity);
		//weight = alignment_length-mismatches;
		//fprintf(chainer, "#%d\n", weight);
		//fprintf(chainer, "#%d\n", alignment_length);
		
		/*
			wrong yass
		*/
		if ( type_of_software == 1 )
		{
			alignment_length1 = abs(q_start - q_end);
			alignment_length2 = abs(s_start - s_end);
			if ( alignment_length1 + 1 == alignment_length )
				gaps_long_seq = alignment_length1 - alignment_length2;
			else
			if ( alignment_length2 + 1 == alignment_length )
				gaps_long_seq = alignment_length2 - alignment_length1;
			else
			{
				/*printf("Error");
				exit(1);*/
			}
			gaps_commun = (gap - gaps_long_seq)/2;

			if ( type_of_weight == 0 )
			{
				//nb of identities
				fprintf(chainer, "#%d\n", alignment_length-mismatches-gaps_commun-gaps_long_seq);
			}
			else
			if ( type_of_weight == 1 )
			{
				//true alignment length
				fprintf(chainer, "#%d\n", alignment_length-gaps_commun-gaps_long_seq);
			}

			if ( q_start <= q_end )
			{
				fprintf(chainer, "[%d,%d] ", q_start, q_end);
			}
			else
			{
				fprintf(chainer, "[%d,%d] ", q_end, q_start);
			}
			if ( s_start <= s_end )
			{
				fprintf(chainer, "[%d,%d]\n", s_start, s_end);
			}
			else
			{
				fprintf(chainer, "[%d,%d]\n", s_end, s_start);
			}
		}
		else
		{
			if ( type_of_weight == 0 )
			{
				//nb of identities
				int matches = identity/100*alignment_length;
				fprintf(chainer, "#%d\n", matches);
			}
			else
			if ( type_of_weight == 1 )
			{
				//true alignment length
				fprintf(chainer, "#%d\n", alignment_length);
			}

			if ( q_start <= q_end )
			{
				fprintf(chainer, "[%d,%d] ", q_start, q_end);
			}
			else
			{
				fprintf(chainer, "[%d,%d] ", q_end, q_start);
			}
			if ( s_start <= s_end )
			{
				fprintf(chainer, "[%d,%d]\n", s_start, s_end);
			}
			else
			{
				fprintf(chainer, "[%d,%d]\n", s_end, s_start);
			}
		}
		
	}

	fclose(yass_blast);
	fclose(chainer);
	return 1;
}
