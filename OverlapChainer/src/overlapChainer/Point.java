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

/**
 * @author Raluca Uricaru
 *
 */

public class Point implements Comparable<Point>{
	
	/*
	 * the x(A), y(B) coordinates of our point
	 */
	IntIntPair coord;
	
	/*
	 * index of the fragment in FragmentList to which the point corresponds  
	 */
	int indexFragment;
	
	/*
	 * true if it is the lowest corner
	 */
	boolean okLow;
	
	/*
	 * index of point in the list of points ordered by y
	 */
	int indexPoint_y;
	
	public Point (Point p)
	{
		coord = new IntIntPair(p.getCoord().getA(), p.getCoord().getB());
		this.indexFragment = p.indexFragment;
		okLow = p.okLow;
		indexPoint_y = p.indexPoint_y;
	}
	
	public Point (IntIntPair p, int indexFragment, boolean okLow)
	{
		coord = new IntIntPair(p.getA(), p.getB());
		
		this.indexFragment = indexFragment;
		
		this.okLow = okLow;	
		
		indexPoint_y = -1;
	}

	public IntIntPair getCoord() {
		return coord;
	}

	public void setCoord(IntIntPair coord) {
		this.coord = coord;
	}

	public int getIndexFragment() {
		return indexFragment;
	}

	public void setIndexFragment(int indexFragment) {
		this.indexFragment = indexFragment;
	}

	public boolean isOkLow() {
		return okLow;
	}

	public void setOkLow(boolean okLow) {
		this.okLow = okLow;
	}

	public int compareTo(Point arg0) {
		// TODO Auto-generated method stub
		
		return new PointYComparator().compare(this, arg0);
	}

	public int getIndexPoint_y() {
		return indexPoint_y;
	}

	public void setIndexPoint_y(int indexPoint_y) {
		this.indexPoint_y = indexPoint_y;
	}

}
