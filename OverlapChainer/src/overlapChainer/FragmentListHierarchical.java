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
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

/**
 * @author Raluca Uricaru
 *
 */

public class FragmentListHierarchical extends AbstractFragmentListNotRearranged {
	
	/*
	 * choose the level till which you test the circularity
	 */
	int level=-1;
	
	public FragmentListHierarchical(File fragmentFile) {
		super(fragmentFile);
		
		init();
	}
	
	/*
	 * in case of java heap space error
	 */
	public FragmentListHierarchical(File fragmentFile, boolean  ok) {
		super(fragmentFile, ok);
		
		init();
	}

	
	public void init(){
		/*
		 * compute the number of levels and their corresponding limits
		 */
		setLevels();
		
		/* create following:
		 * 	list of fragments for each level
		 * 	(they are kept empty if there are too many fragments > NO_MAX_FRAG)
		 */
		this.Fragments_k = new Vector<List<Integer>>(noRecLevels);
		this.Fragments_k.setSize(noRecLevels);
		
		this.noFragOnLevel = new Vector<Integer>(noRecLevels);
		this.noFragOnLevel.setSize(noRecLevels);
		
		for (int i=0; i<noRecLevels; i++)
		{
			this.Fragments_k.setElementAt(new ArrayList<Integer>(), i);
			this.noFragOnLevel.set(i, new Integer(0));
		}
		
		/*
		 * size of shifts on each step
		 * shift = the position at which a fragment is found in a sequence
		 * 		   minus the position in the other sequence
		 * sequence = the window that we want to align on a certain level
		 */
		limitShiftOnStep = new Integer[noStepsOnLevel];
		
		/*
		 * build every level fragment list (if not too many fragments)
		 * attribute to each fragment a recursion level
		 */
		selectFragments_k();
	}
	

	/*
	 * once we've done all tests in order to be sure of the circularity
	 * we can search for the cutting point
	 */
	@Override
	public boolean findCuttingPoint() {
		int cuttingPoint = -1;
		/*
		 * search the first significant fragment at the beginning of the first genome
		 * this could be the cutting point
		 */
		Fragment fFirst = Fragments.get(0);
		int idFirst=0;
		for (Fragment fragment : Fragments )
		{
			if ( fragment.getLevel() <= level )
			{
				fFirst = fragment;
				break;
			}	
			idFirst++;
		}
		
		/*
		 * search the first significant fragment at the end of the first genome
		 * this could be the cutting point
		 */
		Fragment fLast = Fragments.get(Fragments.size()-1);
		Fragment fragment;
		int idLast = Fragments.size()-1;
		boolean ok = false;
		for (int i=Fragments.size()-1; i>=0; i--)
		{
			fragment = Fragments.get(i);
			/*
			 * we need the distance on the second genome (between the 2 cutting points)
			 * smaller than 5% of the "length" of the second genome
			 * 
			 * we need the first cutting point after the second one
			 * on the second genome
			 * 
			 * we search a fragment with these properties
			 * with level <= level
			 */
			if ( fragment.getLevel() <= level && fFirst.getX2() >= fragment.getY2() &&
				 fFirst.getX2()-fragment.getY2()< (maxY2*5/100) )
			{
				idLast = i;
				fLast = fragment;
				ok = true;
				break;
			}	
		}
		
		//TODO
		if ( ok )
			cuttingPoint = idFirst;
		
		if (Constants.DEBUG)
		{
			System.out.println("cut point " + cuttingPoint + " (" + fFirst.getX1() + "," + fFirst.getY1() + ") (" + fFirst.getX2() + "," + fFirst.getY2()+")");
			System.out.println("cut point2 " + idLast + " (" + fLast.getX1() + "," + fLast.getY1() + ") (" + fLast.getX2() + "," + fLast.getY2()+")");
		}
			
		if ( cuttingPoint != -1 )
			return rearrangeFragments(cuttingPoint, idLast);
		return false;
	}

	
	
	/*
	 * compute the number of levels and their corresponding limits
	 */
	@Override
	public void setLevels() {
		
		int lev = 0;
		noRecLevels = 0;
		limitRecLevels_list = new ArrayList<Integer>();
		limitPercentage_list = new ArrayList<Integer>();
		
		int fWMax = Constants.fragWeightMax;
		/*if ( fWMax >= 100000 )
		{
			lev = fWMax/10000 * 10000;
			while ( lev >= 100000)
			{
				limitRecLevels_list.add(new Integer(lev));
				limitPercentage_list.add(new Integer(20));
				lev -= 10000;
				noRecLevels++;
			}
			fWMax = 99999;
		}*/
		
		if ( fWMax >= 10000 )
		{
			lev = fWMax/50000 * 50000;
			while ( lev >= 0)
			{
				if ( lev == 0 ) lev = 10000;
				limitRecLevels_list.add(new Integer(lev));
				limitPercentage_list.add(new Integer(20));
				lev -= 50000;
				noRecLevels++;
			}
			fWMax = 9999;
		}
		
		/*if (fWMax >= 10000)
		{
			lev = fWMax/2500 * 2500;
			while ( lev >= 10000)
			{
				limitRecLevels_list.add(new Integer(lev));
				limitPercentage_list.add(new Integer(20));
				lev -= 2500;
				noRecLevels++;
			}
			fWMax = 9999;
		}*/
		
		
		
		
		if (fWMax >= 1000)
		{
			lev = fWMax/1000 * 1000;
			while ( lev >= 1000)
			{
				limitRecLevels_list.add(new Integer(lev));
				limitPercentage_list.add(new Integer(15));
				lev -= 1000;
				noRecLevels++;
			}
			fWMax = 999;
		}	
		
		
		if (fWMax >= 100)
		{
			lev = fWMax/100 * 100;
			while ( lev >= 100)
			{
				limitRecLevels_list.add(new Integer(lev));
				limitPercentage_list.add(new Integer(10));
				lev -= 100;
				noRecLevels++;
			}
			fWMax = 99;
		}
		
		if (fWMax >= 10)
		{
			lev = fWMax/10 * 10;
			while ( lev >= 10)
			{
				limitRecLevels_list.add(new Integer(lev));
				limitPercentage_list.add(new Integer(5));
				lev -= 10;
				noRecLevels++;
			}
		}
		
		noRecLevels++;
		limitRecLevels_list.add(new Integer(0));
		limitPercentage_list.add(new Integer(1));
		
		/*for ( int i=0; i<noRecLevels; i++ )
			System.out.println(i+ "  " + limitRecLevels_list.get(i));*/
	}
	
	/*
	 * set the maximum size of shifts for every step of level <<level>>
	 */
	@Override
	public void setShifts(int level) {
		
		int limitRecLev = limitRecLevels_list.get(level).intValue();
		
		for (int step=0; step<AbstractFragmentListNotRearranged.noStepsOnLevel-1; step++)
			limitShiftOnStep[step] = new Integer(limitRecLev*(int)Math.pow(2, step));
		
		//alow any size of shifting
		limitShiftOnStep[AbstractFragmentListNotRearranged.noStepsOnLevel-1] = new Integer(5000000);

	}
	
	/*
	 * preliminary test of circularity
	 * histogram of shift lengths
	 * if bimodal and peak points are complementary on the size of the genome, then 
	 * we say that the genome "COULD BE" circular and cut in different positions
	 * 
	 * so, we continue testing
	 */
	protected boolean histogramTestCircular()
	{
		boolean okCircular = false;
		int max_level = noRecLevels/2+1;
		
		/*
		 * build list with shifting lengths for every fragment
		 */
		int pas=0, dec;
		int xLow=0, xHigh=0;
		int prev_xLow=0, prev_xHigh=0;
		ArrayList<Integer> shiftLengths = new ArrayList<Integer>();
		System.out.println();
		
		for ( Fragment fragment: Fragments)
		{
			if ( fragment.getLevel() <= max_level )
			{
				dec = fragment.getX1()-fragment.getX2();
				if ( pas == 0 )
				{
					xLow = dec;
					xHigh = dec;
					pas = 1;
				}
				else
				{
					if ( dec < xLow )
						xLow = dec;
					if ( dec > xHigh )
						xHigh = dec;
				}
				System.out.print(dec + " " );
				shiftLengths.add(new Integer(dec));
			}
		}
		System.out.println();
		System.out.println("xlow = " + xLow + " xhigh = " + xHigh);
		
		prev_xHigh = xHigh;
		prev_xLow = xLow;
		if ( xLow < 0 )
		{
			/*if ( Math.abs(xLow) < Math.abs(xHigh) )
				xLow = 0;
			else
			{
				xHigh = Math.abs(xLow);
				xLow = 0;
			}*/
	
			xHigh = Math.abs(xLow) + Math.abs(xHigh);
			xLow = 0;
		}
		System.out.println("xlow = " + xLow + " xhigh = " + xHigh);
		
		/*
		 * build histogram of shiftings
		 */
		int nBins = 20;
		
		Histogram hist = new Histogram(nBins, xLow, xHigh);
		
		for ( Integer shift : shiftLengths )
		{
			if ( prev_xLow < 0 )
				hist.setData(shift.floatValue() + Math.abs(prev_xLow));
			else
				hist.setData(shift.floatValue());
		}
		
		/*
		 * test whether we have a bimodal histogram
		 */
		
		/*
		 * find the 2 modes
		 */
		int max1=-1, max2=-1;
		int poz1=-1, poz2=-1;
		
		for ( int i=0; i<nBins; i++ )
		{
			System.out.println("(" + i + ", " + hist.getBins()[i] + ")");
			if ( hist.getBins()[i] > max1 )
			{
				max2 = max1;
				poz2 = poz1;
				
				max1 = hist.getBins()[i];
				poz1 = i;
			}
			else
			if ( hist.getBins()[i] > max2 )
			{
				max2 = hist.getBins()[i];
				poz2 = i;
			}	
		}
		
		System.out.println("bin1 = " + poz1 + " bin2 = " + poz2 + " " + (int)hist.getDelBin());
		
		/*
		 * test whether they are complementary on genome length
		 */
		int dec1_bin =  Math.abs(prev_xLow + poz1 * (int)hist.getDelBin());
		int dec2_bin =  Math.abs(prev_xLow + poz2 * (int)hist.getDelBin());
		
		//= 20% genome length
		int eps;
		int Y_max;
		
		if ( Fragments.get(this.noFragments-1).getY1() > Fragments.get(this.noFragments-1).getY2() )
			Y_max = Fragments.get(this.noFragments-1).getY1();
		else
			Y_max = Fragments.get(this.noFragments-1).getY2();
		
		eps = Y_max/5;
		
		System.out.println("dec1_bin = " + dec1_bin + " dec2_bin = " + dec2_bin + " Y_max "+Y_max + " eps " + eps);

		if ( dec1_bin+dec2_bin <= Y_max+eps
			 && dec1_bin+dec2_bin >= Y_max-eps )	
			okCircular = true;
		
		return okCircular;
	}
	
	
	/*
	 * 2nd test
	 * test whether we deal with circular genomes cut in different positions; 
	 * find the cutting point
	 */
	@SuppressWarnings("unchecked")
	@Override
	public int testCircular() {
		
		/*
		 * we need at least 50 and at most 100 fragments
		 * in order to test circularity
		 */
		int no_max_frag = 100;
		int no_min_frag = 50;
		
		if ( Fragments.size()<50 )
			no_min_frag = Fragments.size();
		
		int okCircular_fragment = -1;
		Fragment fragment, prev_fragment;
		int begin_g2;
		LinkedHashMap<Integer,Pair<Integer,Double>> NoRearrangements = new LinkedHashMap<Integer,Pair<Integer,Double>>(Fragments.size());
	
		int nof = 0;
		while ( nof<no_max_frag && level < noRecLevels-1 )
		{
			level++;
			nof += this.noFragOnLevel.get(level);	
		}
		
		if ( nof > no_max_frag )
		{
			level--;
			nof -= this.noFragOnLevel.get(level+1);
		}	
		if ( nof < no_min_frag && level < noRecLevels)
		{
			level++;
			nof += this.noFragOnLevel.get(level);			
		}
		
		if ( Constants.DEBUG )
		{
			System.out.println("\nTest circularity till level " + level);
			System.out.println(nof + " fragments tested of " + Fragments.size());
		}
		
		int countF, countT;
		for (int i=0; i<Fragments.size(); i++){
			
			fragment = this.Fragments.get(i);
			if ( fragment.getLevel() > level ) continue;
			
			begin_g2 = fragment.getX2();
			
			/*
			 * count fragments that we find before "fragment" in the first genome
			 * but after it, in the second one
			 */
			countF = countT = 0;
			for (int j=0; j<i; j++)
			{
				prev_fragment = this.Fragments.get(j);
				if ( prev_fragment.getLevel() >= level ) continue;
				
				countT++;
				if ( prev_fragment.getX2() > begin_g2 )
					countF++;
			}
				
			//System.out.print("(" + countF+ ", " + countT + ", " + fragment.getX1()+")");
	
			if ( countT>0 )
				NoRearrangements.put(new Integer(i), new Pair(new Integer(countT), new Double((double)countF/countT)));
			else
				NoRearrangements.put(new Integer(i), new Pair(new Integer(0), new Double(0)));
		}
		
		/*
		 * decide whether we deal with the problem of circularity or not
		 */
		okCircular_fragment = -1;
		double epsilon = 0.1;
		double val_max = 0;
		int fmax = 0;
		int NoBigRearrangements = 0;
		
		for (Map.Entry<Integer,Pair<Integer,Double>> e : NoRearrangements.entrySet() )
		{
			if ( e.getValue().getB() >= 0.9 )
				NoBigRearrangements++;
			
			if ( okCircular_fragment == -1 )
			{
				val_max = e.getValue().getB();
				okCircular_fragment = e.getKey();
				fmax = e.getValue().getA();
				continue;
			}
			
			if ( e.getValue().getB() > val_max )
			{
				if ( e.getValue().getB() == 1)
				{
					val_max = e.getValue().getB();
					okCircular_fragment = e.getKey();
					fmax = e.getValue().getA();
				}
				else
				if ( e.getValue().getB() - val_max <= epsilon && 
						 Fragments.get(okCircular_fragment).getX2() <
						 Fragments.get(e.getKey()).getX2() )
						continue;
				else
				{
					val_max = e.getValue().getB();
					okCircular_fragment = e.getKey();
					fmax = e.getValue().getA();
				}
			}
			else
			{
				if ( e.getValue().getB() == 1 && 
						Fragments.get(okCircular_fragment).getX2() >
						Fragments.get(e.getKey()).getX2()  )
				{
					val_max = e.getValue().getB();
					okCircular_fragment = e.getKey();
					fmax = e.getValue().getA();
				}
				else
				if ( val_max - e.getValue().getB() <= epsilon && 
						 Fragments.get(okCircular_fragment).getX2() >
						 Fragments.get(e.getKey()).getX2()   )	
				{
					val_max = e.getValue().getB();
					okCircular_fragment = e.getKey();
					fmax = e.getValue().getA();
				}
			}
			//System.out.print(e.getValue().getB() + "=(" + Fragments.get(e.getKey()).getX1() + ")");
		}
		
		if ( okCircular_fragment == 0 )
			okCircular_fragment = -1;
		
		
		if ( (okCircular_fragment != -1 && val_max==1 && NoBigRearrangements >=2) || (val_max>=0.9 && NoBigRearrangements >=2) )
		{
			/*System.out.println("\n\nProblems due to circularity. Wrong cut: fragment " + 
								Fragments.get(okCircular_fragment));*/
		}
		else
		{
			//System.out.println("\n\nNO problems due to circularity");
			okCircular_fragment = -1;
		}
		
		NoRearrangements.clear();
		return okCircular_fragment;
	}

}
