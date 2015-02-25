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

public class Pair<X,Y> {
	private X a;
	private Y b;

	/**
	 * @param a
	 * @param b
	 */
	public Pair(X a, Y b) {
		super();
		this.a = a;
		this.b = b;
	}
	
	public X getA() {
		return a;
	}
	public void setA(X a) {
		this.a = a;
	}
	public Y getB() {
		return b;
	}
	public void setB(Y b) {
		this.b = b;
	}
	
	
}
