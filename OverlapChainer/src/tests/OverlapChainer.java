package tests;
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

import overlapChainer.*;

/**
 * @author raluca
 *
 */

public class OverlapChainer {
	
	public static void main(String[] args) throws Throwable{

		if ( args.length > 2 )
			throw new Exception("Wrong number of parameters!");

		
		int overlap_ratio = 10; //overlap ratio default used if only one parameter given instead of 2
		if ( args.length > 1 )
		{
			// if overlap_ratio >= 100 then it is a maximum constant overlap and not a ratio
			overlap_ratio = new Integer(args[0]).intValue();
		}
		
		FragmentListSweepline FragmentL;
		try{
			FragmentL = new FragmentListSweepline(new File(args[1]), overlap_ratio);
		}
		catch (OutOfMemoryError e)
		{
			/*
			 * Java heap space overlapped
			 * do not read all the fragments from the .mat file
			 * we take MAX_NB_FRAG 
			 */
			System.out.println("Java Heap Space Error");
			FragmentL = null;
			System.gc();
			FragmentL = new FragmentListSweepline(new File(args[1]), false, overlap_ratio);
		}
				
				
		/*
		 * test whether we deal with circular genomes cut in different positions
		 * find the cut point (for first levels)
		 */
			/*
			 * assure circularity
			 */
			if ( FragmentL.testCircular() != -1)
			{	
				/*
				 * we deal with wrong cutting point of circular genomes
				 * -> rearrange fragments
				 */	 
				if ( FragmentL.findCuttingPoint() )
				{
					//System.out.println("Circular");
				}
			}
		
		long t1 = System.currentTimeMillis();
			
		FragmentL.init2();		
		FragmentL.setLevels();
		FragmentL.selectFragments_k();
		
		/*
		 * Passing the sweep line
		 */
		FragmentL.chainWithOverlaps();
		
		long t2 = System.currentTimeMillis();
		
		if (FragmentL.getFragments().get(FragmentL.getChain().get(0)).getMovingDirection() != 0 )
			System.out.println(">H 2 Circular " + FragmentL.getLeftShift() + " " + FragmentL.getRightShift() + "\n#0");
		else
			System.out.println(">H 2\n#0");
		
		
		//the normal output
		FragmentL.Chain_toString();	
		
		float t = (float)(t2-t1) / 1000F;
		System.out.println("total time: " + t + "\n\n");
	}
}



