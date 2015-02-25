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

public class Marker {
	/*
	 * the index of the Fragment to which the marker is attached to in the FragmentList
	 */
	int originIndex;
	
	/*
	 * the index of the Point to which the marker is positioned on in the list of Points ordered on x
	 */
	int pointIndex;
	
	/*
	 * marked or not
	 */
	boolean marked;
	
	/*
	 * chosen to be deleted or not
	 */
	boolean semi_deleted;
	
	/*
	 * pointed to or not by another fragment
	 */
	boolean pointed;
	
	/*
	 * first marker position before in the list of markers
	 */
	int markerBefore;
	
	/*
	 * first marker position after in the list of markers ordered on y
	 */
	int markerAfter;
	
	
	/*
	 * weight
	 */
	int weight;
	
	public Marker()
	{
		pointIndex = -1;
		originIndex = -1;
		marked = false;
		markerBefore = -1;
		markerAfter = -1;
		semi_deleted = false;
		pointed = false;
	}
	
	
	public Marker(int point, int origin)
	{
		pointIndex = point;
		originIndex = origin;
		marked = false;
		markerBefore = -1;
		markerAfter = -1;
		semi_deleted = false;
		pointed = false;
	}
	
	
	public int getOriginIndex() {
		return originIndex;
	}

	public void setOriginIndex(int originIndex) {
		this.originIndex = originIndex;
	}

	public int getPointIndex() {
		return pointIndex;
	}

	public void setPointIndex(int pointIndex) {
		this.pointIndex = pointIndex;
	}

	public boolean isMarked() {
		return marked;
	}

	public void setMarked(boolean marked) {
		this.marked = marked;
	}

	public int getMarkerBefore() {
		return markerBefore;
	}

	public void setMarkerBefore(int markerBefore) {
		this.markerBefore = markerBefore;
	}

	public int getMarkerAfter() {
		return markerAfter;
	}

	public void setMarkerAfter(int markerAfter) {
		this.markerAfter = markerAfter;
	}


	public int getWeight() {
		return weight;
	}


	public void setWeight(int weight) {
		this.weight = weight;
	}


	public boolean isSemi_deleted() {
		return semi_deleted;
	}


	public void setSemi_deleted(boolean semi_deleted) {
		this.semi_deleted = semi_deleted;
	}


	public boolean isPointed() {
		return pointed;
	}


	public void setPointed(boolean pointed) {
		this.pointed = pointed;
	}
	
	
}
