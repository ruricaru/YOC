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

import java.util.Comparator;

/**
 * @author Raluca Uricaru
 *
 */

public class PointXComparator implements Comparator<Point>
{
	public final int compare ( Point p1, Point p2 )
	{
		int compare = p1.getCoord().getA()-p2.getCoord().getA();

		if (compare == 0)
		{
			if (p1.okLow && p2.okLow || !p1.okLow && !p2.okLow)
			{
				return p1.getCoord().getB()- p2.getCoord().getB();
			}
			else
				if ( p1.okLow )
					return -1;
				else 
					return 1;
		}
		else
			return compare;

	} // end compare

}
