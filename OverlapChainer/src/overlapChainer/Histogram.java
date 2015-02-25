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

//********************************************************************
//The setData(float f) finds in which bin the value falls for nBins between the given minimum and maximum values. 
//An integer array keeps track of the number of times the input
//value fell into a particular bin.
//********************************************************************	


//********************************************************************
public class Histogram {

	private int [] bins = null;
	private int nBins;
	private float xLow, xHigh;
	private float delBin;
	
	
	//----------------------------------------------------------------
	Histogram (int nBins, float xLow, float xHigh){
		
		this.nBins = nBins;
		this.xLow  = xLow;
		this.xHigh = xHigh;
		
		bins = new int[nBins];
		delBin = (xHigh-xLow)/(float)nBins;
	}
	
	//----------------------------------------------------------------
	// Extra constructor to allow for double values
	Histogram (int nBins, double xLow, double xHigh,int _width,int _height){
		this(nBins, (float) xLow, (float) xHigh);
	}
	
	//----------------------------------------------------------------
	void setData(double data){
		setData((float)data);
	}
	
	//----------------------------------------------------------------
	void setData(float data){
	
		if( data >= xLow && data < xHigh){
			int bin = (int)((data-xLow)/delBin);
			if(bin >=0 && bin < nBins) bins[bin]++;
		}
	}

	public int[] getBins() {
		return bins;
	}

	public float getDelBin() {
		return delBin;
	}    
	
}
