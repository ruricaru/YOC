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
 * @author raluca 
 *
 */

@SuppressWarnings("rawtypes")
public class Fragment implements Comparable{
	/*
	 * fragment limits in the first genome
	 */
	protected int x1;
	protected int y1;
	/*
	 * fragment limits in the second genome
	 */
	protected int x2;
	protected int y2;
	
	protected int weight;
	
	/*
	 * id% de YASS not used in this version
	 */
	protected double idp;
	
	protected int length1;
	protected int length2;
	
	/*
	 * index of original Fragments before sorting by x1
	 */
	//private int initialIndex;
	
	/*
	 * movingDirection in case of rearranging fragments
	 *  -1 = if it has been moved to the left
	 *  	 (use +leftShift on x2 and y2 in order to obtain the true values)
	 *  1  = if it has been moved to the right
	 *  	 (use -rightShift on x2 and y2 in order to obtain the true values)
	 *  0  = if it has not been moved
	 */
	protected int movingDirection;
	
	/*
	 * recursion level to which it corresponds
	 */
	protected int level;

	
	/*
	 * only for the sweep line
	 * weight of the best chain arriving in this fragment 
	 */
	private int weightChain;
	private int link;
	
	public Fragment() {
		super();
		this.x1 = 0;
		this.y1 = 0;
		this.x2 = 0;
		this.y2 = 0;
		this.weight = 0;
		//this.initialIndex = -1;
		this.movingDirection = 0;
		weightChain = weight;
		link = -1;
	}
	
	/**
	 * @param x1
	 * @param y1
	 * @param x2
	 * @param y2
	 */
	public Fragment(int x1, int y1, int x2, int y2) {
		super();
		this.x1 = x1;
		this.y1 = y1;
		this.x2 = x2;
		this.y2 = y2;
		this.weight = 0;
		//this.initialIndex = -1;
		this.movingDirection = 0;
		weightChain = weight;
		link = -1;
	}
	
	public Fragment(Fragment f) {
		super();
		this.x1 = f.x1;
		this.y1 = f.y1;
		this.x2 = f.x2;
		this.y2 = f.y2;
		this.length1 = f.length1;
		this.length2 = f.length2;
		this.weight = f.weight;
		this.level = f.level;
		//this.initialIndex = f.initialIndex;
		this.movingDirection = f.movingDirection;
		weightChain = f.weightChain;
		link = f.link;
	}
	
	public int getX1() {
		return x1;
	}
	public void setX1(int x1) {
		this.x1 = x1;
	}
	public int getX2() {
		return x2;
	}
	public void setX2(int x2) {
		this.x2 = x2;
	}
	public int getY1() {
		return y1;
	}
	public void setY1(int y1) {
		this.y1 = y1;
	}
	public int getY2() {
		return y2;
	}
	public void setY2(int y2) {
		this.y2 = y2;
	}
	
	public int getMovingDirection() {
		return this.movingDirection;
	}

	public void setMovingDirection(int movingDirection) {
		this.movingDirection = movingDirection;
	}

	public int getWeight() {
		return weight;
	}

	public void setWeight(int weight) {
		this.weight = weight;
		
		if ( weight > Constants.fragWeightMax )
			Constants.fragWeightMax = weight;
		
		if ( weight < Constants.fragWeightMin )
			Constants.fragWeightMin = weight;
			
	}
	
	public void computeLengths() {
		//TODO avant je mettais pas le +1
		this.length1 = this.y1 - this.x1 + 1;
		this.length2 = this.y2 - this.x2 + 1;
	}

	public int getLength1() {
		return length1;
	}

	public int getLength2() {
		return length2;
	}
	
	public int getLevel() {
		return level;
	}

	public void setLevel(int level) {
		this.level = level;
	}
	
	/*public int getInitialIndex() {
		return initialIndex;
	}

	public void setInitialIndex(int initialIndex) {
		this.initialIndex = initialIndex;
	}*/
	

	public int getWeightChain() {
		return weightChain;
	}

	public void setWeightChain(int weightChain) {
		this.weightChain = weightChain;
	}

	
/*	public int compareTo(Object o){
		if (! (o instanceof Fragment) )
			throw new ClassCastException();
		Fragment f = (Fragment)o;
		return (new Integer(this.x1)).compareTo(new Integer(f.getX1()));
	}
*/	
	
	public int compareTo(Object o){
		if (! (o instanceof Fragment) )
			throw new ClassCastException();
		Fragment f = (Fragment) o;
		return this.x1 - f.getX1();
		//return (new Integer(this.x1)).compareTo(new Integer(f.getX1()));
	}
	
	public boolean equals(Object o){
		if (! (o instanceof Fragment) )
			throw new ClassCastException();
		Fragment f = (Fragment)o;
		if (this.x1==f.x1 && this.x2==f.x2 && this.y1==f.y1 && this.y2==f.y2 && this.weight==f.weight)
			return true;
		else
			return false;
	}
	
	public String toString(){
		String fr_str = new String();
		fr_str = "("+this.x1 +","+this.y1+") ("+ +this.x2 +","+this.y2+")";
		return fr_str;
	}

	public double getIdp() {
		return idp;
	}

	public void setIdp(double idp) {
		this.idp = idp;
	}

	public int getLink() {
		return link;
	}

	public void setLink(int link) {
		this.link = link;
	}
	
	public boolean almostInferiorTo(Fragment f)
	{
		boolean ok = true;
		
		int dif1 = f.x1 - this.y1;
		int dif2 = f.x2 - this.y2;
		
		int lev;
		if ( f.level < this.level )
			lev = f.level;
		else
			lev = this.level;
		
		if ( dif1 <= 0 )
		{	
			/*
			 * verify the overlap with the pre-fragment in the chain
			 */
			
			if ( !(Math.abs(dif1) <= FragmentListBacktrack.percentage(f.length1, lev) && 
				 Math.abs(dif1) <= FragmentListBacktrack.percentage(this.length1, lev)))
			{
				ok = false;
			}
		}
		
		if ( dif2 <= 0 ) 
		{
			/*
			 * verify the overlap with the pre-fragment in the chain
			 */
			
			if ( !(Math.abs(dif2) <= FragmentListBacktrack.percentage(f.length2, lev) && 
					Math.abs(dif2) <= FragmentListBacktrack.percentage(this.length2, lev)))
			{
				ok = false;
			}
		}
		
		return ok;
	}
	
	
	/*
	 * used for dealing with both proportional and constant overlaps
	 * 
	 * if maximum overlap allowed >= 100 then it is a maximum constant overlap and not a ratio
	 */
	public boolean almostInferiorToUpgraded(Fragment f)
	{
		boolean ok = true;
		
		int dif1 = f.x1 - this.y1;
		int dif2 = f.x2 - this.y2;
		
		int lev;
		if ( f.level < this.level )
			lev = f.level;
		else
			lev = this.level;
		
		if ( f.y1 <= this.y1 || f.y2 <= this.y2 || this.x1 >= f.x1 || this.x2 >= f.x2 )
			return false;
		
		if ( dif1 <= 0 )
		{	
			/*
			 * verify the overlap with the pre-fragment in the chain
			 */
			
			if ( !(Math.abs(dif1) <= FragmentListBacktrack.percentage(f.length1, lev) && 
				 Math.abs(dif1) <= FragmentListBacktrack.percentage(this.length1, lev) &&
				 dif1 < f.length1 && dif1 < this.length1) )
			{
				ok = false;
			}
		}
		
		if ( dif2 <= 0 ) 
		{
			/*
			 * verify the overlap with the pre-fragment in the chain
			 */
			
			if ( !(Math.abs(dif2) <= FragmentListBacktrack.percentage(f.length2, lev) && 
					Math.abs(dif2) <= FragmentListBacktrack.percentage(this.length2, lev)  &&
					 dif2<f.length2 && dif2<this.length2))
			{
				ok = false;
			}
		}
		
		return ok;
	}
	
	/*
	 * Not used anymore
	 */
	public boolean inferiorButOverlapped(Fragment f)
	{
		boolean ok = true;
		
		int dif1 = f.x1 - this.y1;
		int dif2 = f.x2 - this.y2;
		
		int lev;
		if ( f.level < this.level )
			lev = f.level;
		else
			lev = this.level;
		
		if ( dif1 <= 0 )
		{	
			/*
			 * verify the overlap with the pre-fragment in the chain
			 */
			if ( !(Math.abs(dif1) <= FragmentListBacktrack.percentage(f.length1, lev) && 
				 Math.abs(dif1) <= FragmentListBacktrack.percentage(this.length1, lev)))
				ok = false;
		}
		
		if ( dif2 <= 0 ) 
		{
			/*
			 * verify the overlap with the pre-fragment in the chain
			 */
			if ( !(Math.abs(dif2) <= FragmentListBacktrack.percentage(f.length2, lev) && 
					Math.abs(dif2) <= FragmentListBacktrack.percentage(this.length2, lev)))
				ok = false;
		}
		
		if ( ok == false )
			return false;
		else
			if (dif1 <= 0|| dif2 <= 0)
				return true;
			else
				return false;
	}
	
	/*
	 * Not used anymore
	 */
	public int getBiggestAcceptedOverlap()
	{
		int overlap = 0;
		overlap = FragmentListBacktrack.percentage(this.length1, level) + FragmentListBacktrack.percentage(this.length2, level);
		return overlap;
	}
}
