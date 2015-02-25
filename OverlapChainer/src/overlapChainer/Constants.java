package overlapChainer;
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

import java.io.File;

/**
 * @author Raluca Uricaru
 *
 */

public class Constants {
	public static String fasta1, fasta2;
	public static File circular_genomes;
	
	public static final boolean DEBUG = false;
	public static final boolean OUTPUT_COMPLETE = false;
	public static final int NO_MAX_FRAG = 100000;
	
	public static int fragWeightMin = 10000000;
	public static int fragWeightMax = 0;	
}
