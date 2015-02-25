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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.Vector;

/**
 * @author Raluca Uricaru
 *
 */

public class FragmentListSweepline extends FragmentListHierarchical {
	
	/*
	 *  structures for sweep line algorithm
	 *  keep points ordered on x
	 */
	ArrayList<Point> Points;
	
	/*
	 * we need a list of markers (activated or not)
	 * 
	 * we shall keep the markers in a vector of size = 2*nb_of_fragments = nb_of_points
	 * for reallocation problems
	 * 
	 * keep markers ordered on y
	 */
	Vector<Marker> Markers;
	
	/*
	 * overlap ratio
	 */
	int overlap_ratio;
	
	/*
	 * choose the level till which you test the circularity
	 */
	int level=-1;
	
	public FragmentListSweepline(File fragmentFile, int overlap_ratio_) {
		super(fragmentFile);	
		overlap_ratio = overlap_ratio_;
	}
	
	/*
	 * in case of java heap space error
	 */
	public FragmentListSweepline(File fragmentFile, boolean  ok, int overlap_ratio_) {
		super(fragmentFile, ok);	
		overlap_ratio = overlap_ratio_;
	}

	
	public void init(){
		
		/*
		 * compute the number of levels and their corresponding limits but only for 
		 * determining whether genomes are circular or not
		 * for the rest, levels are given as a parameter of the program
		 */
		setLevelsForCircularity();
		
		this.noFragOnLevel = new Vector<Integer>(noRecLevels);
		this.noFragOnLevel.setSize(noRecLevels);
		
		for (int i=0; i<noRecLevels; i++)
		{
			this.noFragOnLevel.set(i, new Integer(0));
		}
	
		
		/*
		 * build every level fragment list (if not too many fragments)
		 * attribute to each fragment a recursion level
		 */
		selectFragments_k();
	}
	
	
	/*
	 * compute the number of levels and their corresponding limits
	 */
	public void setLevelsForCircularity() {
		
		int lev = 0;
		noRecLevels = 0;
		limitRecLevels_list = new ArrayList<Integer>();
		limitPercentage_list = new ArrayList<Integer>();
		
		int fWMax = Constants.fragWeightMax;
		if ( fWMax >= 100000 )
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
		}
		
		if (fWMax >= 10000)
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
		}
	
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
	 * for the hierarchical chaining based on the sweep line
	 * we do not fix several levels: one level with one overlap ratio
	 */
	@Override
	public void setLevels() {
		
		noRecLevels = 1;
		limitRecLevels_list = new ArrayList<Integer>();
		limitPercentage_list = new ArrayList<Integer>();
		
		limitRecLevels_list.add(new Integer(0));
		limitPercentage_list.add(new Integer(overlap_ratio));
	}
	
	
	public void init2()
	{
		/*
		 * builds points and markers list
		 */
		Points = new ArrayList<Point>(Fragments.size()*2);
		Markers = new Vector<Marker>(Fragments.size()*2);
		Markers.setSize(Fragments.size()*2);
		
		
		/*
		 * create points list directly ordered on y
		 */
		int index = 0;
		for (Fragment f : Fragments)
		{
			Point pLow = new Point(new IntIntPair(f.getX1(), f.getX2()), index, true);
			Point pHigh = new Point(new IntIntPair(f.getY1(), f.getY2()), index, false);
			index++;
			
			int i = Collections.binarySearch(Points, pLow);
			
			if ( i < 0 )
				Points.add(-i-1, pLow);
			else
				Points.add(i, pLow);
			
			i = Collections.binarySearch(Points, pHigh);
			
			if ( i < 0 )
				Points.add(-i-1, pHigh);
			else
				Points.add(i, pHigh);
		}
		
		/*
		 * perturbe corner points  on y coordinates
		 * keep index_y for each point
		 */
		int perturbation_y = 0;
		Points.get(0).setIndexPoint_y(0);
		
		for (int i=1; i<Points.size(); i++)
		{
			Point prev = Points.get(i-1);
			Point current = Points.get(i);
			
			current.setIndexPoint_y(i);
			
			if ( prev.getCoord().getB() - perturbation_y == current.getCoord().getB())
			{
				//System.out.println(perturbation_y + " " + Fragments.get(prev.getIndexFragment()) + " " + Fragments.get(current.getIndexFragment()));
				perturbation_y++;
				current.setCoord(new IntIntPair(current.getCoord().getA(), current.getCoord().getB() + perturbation_y));
			}			
			else
			{
				current.setCoord(new IntIntPair(current.getCoord().getA(), current.getCoord().getB()+perturbation_y));
			}
			
		}
		
		Collections.sort(Points, new PointXComparator());

		
		/*
		 * perturbe corner points  on x coordinates
		 *  
		 * +
		 * 
		 * create marker list (with non active markers)
		 */
		int perturbation_x = 0;
		Point prev = Points.get(0);
		
		/*
		 * create non active marker for prev
		 */
		Marker m = new Marker(0, prev.indexFragment);
		m.setWeight(Fragments.get(prev.indexFragment).getWeight());
		Markers.set(prev.getIndexPoint_y(), m);
		
		for (int i=1; i<Points.size(); i++)
		{
			Point current = Points.get(i);
			
			/*
			 * create non active marker for current
			 */
			m = new Marker(i, current.indexFragment);
			m.setWeight(Fragments.get(current.indexFragment).getWeight());
			Markers.set(current.getIndexPoint_y(), m);
	
			
			if (prev.getCoord().getA() - perturbation_x == current.getCoord().getA())
			{
				perturbation_x++;
				current.setCoord(new IntIntPair(current.getCoord().getA()+perturbation_x, current.getCoord().getB()));
				//System.out.println(Fragments.get(prev.getIndexFragment()) + " " + Fragments.get(current.getIndexFragment()));
			}		
			else
			{
				current.setCoord(new IntIntPair(current.getCoord().getA()+perturbation_x, current.getCoord().getB()));
			}
			prev = current;
		}
		
		/*for (int i=0; i<Points.size(); i++)
		{
			Point current = Points.get(i);
			System.out.print("["+current.getCoord().getA()+","+current.getCoord().getB()+"]=" + current.getIndexFragment() + ","+current.okLow+ ","+current.getIndexPoint_y()+"   ");
		}*/
		
		//System.out.println("\nFinished initialisations...");
	
	}
	
	
	/*
	 * recursive step k
	 * fragments with weight bigger than limit on step k
	 * (fragment list previously ordered by their initial position in the first genome)
	 */
	public void selectFragments_k ()
	{
		Fragment fragment = new Fragment();
		int k;
		int size1, size2;
		
		
		for ( int index=0; index<Fragments.size(); index++ ){
			fragment = this.Fragments.get(index);
			/*
			 * we try to find the proper recursion level for the fragment
			 */
			if ( noRecLevels == 1 )
			{
				this.Fragments.get(index).setLevel(0);
				this.noFragOnLevel.set(0, this.noFragOnLevel.get(0).intValue()+1);
			}
			else
			{
				for (k=0; k<noRecLevels; k++)
				{
					size1 = limitRecLevels_list.get(k).intValue();
					size2 = 10000000;
					if ( k > 0 )
					{
						size2 = limitRecLevels_list.get(k-1).intValue();
								
						if ( fragment.getWeight()>= size1 && fragment.getWeight()<size2 )
						{
							this.Fragments.get(index).setLevel(k);
							this.noFragOnLevel.set(k, this.noFragOnLevel.get(k).intValue()+1);
							break;
						}
					}
				}
			}
		}
	}
	

	
	/*
	 * using the method of the sweepline, we try to find the best chain (don't accept overlaps)
	 */
	public void chainWithoutOverlaps()
	{
		/*
		 * we have all points sorted on x in Points
		 * we have non-activated markers sorted on y in Markers (for which we have access to their position in Points
		 * creating a marker for us means activating a marker in Markers
		 */
		
		
		int firstActivatedMarkerIndex = -1;
		int lastActivatedMarkerIndex = -1;
		
		int nb_markers=0;
		int pace = -1;
		for (Point p : Points)
		{
			pace++;
			/*
			 * without overlaps
			 * 
			 * I "search first activated marker below p" phase
			 */
			
			/*
			 * search the non-activated marker corresponding to p
			 */
			Marker mp = Markers.get(p.getIndexPoint_y());
			
			/*
			 * search the first activated marker m below p
			 */
			Marker m;
			if ( firstActivatedMarkerIndex == -1 || firstActivatedMarkerIndex > p.getIndexPoint_y() )
			{
				m = null;
			}
			else
			{
				int last = lastActivatedMarkerIndex;
				while (last >= p.getIndexPoint_y())
					last = Markers.get(last).getMarkerBefore();
				
				if ( last == -1 )
					m = null;
				else
					m = Markers.get(last);
			}
			
			/*
			 * II "test whether p is a lower corner or an upper corner" phase
			 */
			Fragment frag_p = Fragments.get(p.getIndexFragment());
		
			if ( p.okLow )
			{
				/*
				 * if lower corner
				 */
				if ( m != null )
				{
					frag_p.setWeightChain(frag_p.getLength1()+frag_p.getLength2()+ m.getWeight());
					//frag_p.setWeightChain(frag_p.getWeight() + m.getWeight());
					frag_p.setLink(m.getOriginIndex());
				}
				else
				{
					frag_p.setWeightChain(frag_p.getLength1()+frag_p.getLength2());
					//frag_p.setWeightChain(frag_p.getWeight());
					frag_p.setLink(-1);
				}
			}
			else
			{
				/*
				 * if upper corner
				 */
				if ((m!=null && frag_p.getWeightChain() > m.getWeight()) || m==null)
				{
					/*
					 * activate marker mp
					 */
					nb_markers++;
					mp.setMarked(true);
					mp.setWeight(frag_p.getWeightChain());
					
					if ( m != null )
					{
						mp.setMarkerBefore(Points.get(m.getPointIndex()).getIndexPoint_y());
						mp.setMarkerAfter(m.getMarkerAfter());
						//ICIIIIIIIIII
						if ( m.getMarkerAfter() != -1 )
							Markers.get(m.getMarkerAfter()).setMarkerBefore(p.indexPoint_y);
						m.setMarkerAfter(p.getIndexPoint_y());
					}
					else
					{
						mp.setMarkerBefore(-1);
						mp.setMarkerAfter(firstActivatedMarkerIndex);
						//ICIIIIIIIIII
						if (firstActivatedMarkerIndex != -1)
							Markers.get(firstActivatedMarkerIndex).setMarkerBefore(p.indexPoint_y);
					}
								
					/*
					 * update firstActivatedMarkerIndex, lastActivatedMarkerIndex
					 */
					if (p.getIndexPoint_y() < firstActivatedMarkerIndex || firstActivatedMarkerIndex == -1 )
						firstActivatedMarkerIndex = p.getIndexPoint_y();
					if (p.getIndexPoint_y() > lastActivatedMarkerIndex || lastActivatedMarkerIndex == -1 )
						lastActivatedMarkerIndex = p.getIndexPoint_y();
					
					/*
					 * III "remove all mi higher and lighter than mp" phase
					 */
					int i = mp.getMarkerAfter();
					while (i != -1 && i<=lastActivatedMarkerIndex)
					{
						Marker mi = Markers.get(i);
						if ( mi.getWeight() < mp.getWeight() )
						{
							nb_markers--;
							mi.setMarked(false);
							if ( mi.getMarkerBefore() != -1 )
								Markers.get(mi.getMarkerBefore()).setMarkerAfter(mi.getMarkerAfter());
							
							if ( mi.getMarkerAfter() != -1 )
								Markers.get(mi.getMarkerAfter()).setMarkerBefore(mi.getMarkerBefore());
							
							if ( i== firstActivatedMarkerIndex )
								firstActivatedMarkerIndex = mi.getMarkerAfter();
							
							if ( i == lastActivatedMarkerIndex )
								lastActivatedMarkerIndex = mi.getMarkerBefore();
							
							
							i = mi.getMarkerAfter();
							
							mi.setMarkerAfter(-1);
							mi.setMarkerBefore(-1);
										
						}
						else
							i = mi.getMarkerAfter();
							
					}
				}
				
			}
		}
		
		try{
			BufferedWriter bf = new BufferedWriter (new FileWriter("markers", true));
			bf.append(nb_markers+" "+fragmentFile.getName()+"\n");
			bf.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		
		int i = firstActivatedMarkerIndex;
		int imax=-1, weightmax=0;
		
		/*while ( i!=-1 && i<= lastActivatedMarkerIndex )
		{
			Marker m = Markers.get(i);
			if ( weightmax < m.getWeight() )
			{
				weightmax = m.getWeight();
				imax = i;
			}
			i = m.getMarkerAfter();
			//System.out.println(Fragments.get(m.getOriginIndex()));
		}*/
		
		imax = lastActivatedMarkerIndex;
		
		/*
		 * backtrack the best chain	
		 */
		i = Markers.get(imax).getOriginIndex();
		while ( i > -1 )
		{
			Chain.add(i);
			i = Fragments.get(i).getLink();
		}
		
		Collections.sort(Chain);
	}
	
	
	/*
	 * using the method of the sweepline, we try to find the best chain (accept overlaps)
	 */
	public void chainWithOverlaps_wrong_deleted_markers()
	{
		/*
		 * we have all points sorted on x in Points
		 * we have non-activated markers sorted on y in Markers (for which we have access to their position in Points
		 * creating a marker for us means activating a marker in Markers
		 */
		
		
		int firstActivatedMarkerIndex = -1;
		int lastActivatedMarkerIndex = -1;
		
		
		int pace = -1;
		
		int nb_markers = 0;
		
		
		/*
		 * list of low points opened and not already closed
		 */
		LinkedList<Integer> opened_lows = new LinkedList<Integer>();

		for (Point p : Points)
		{
			pace++;
			int max_weight = 0;
			int max_index = -1;
			
			/*
			 * search the non-activated marker corresponding to p
			 */
			Marker mp = Markers.get(p.getIndexPoint_y());
			Fragment frag_p = Fragments.get(p.getIndexFragment());
			
		
			/*
			 * update the opened_lows list
			 */
			if (p.isOkLow())
			{
				/*
				 * add the newly opened lows to the list
				 */
				opened_lows.add(p.indexPoint_y);
			}
			
			/*
			 * I "search first activated marker almost below p" phase
			 */
				
			Marker m;
			int marker_index_before_mp = -1;
			
			if ( firstActivatedMarkerIndex == -1 )
			{
				m = null;
			}
			else
			{
				if (p.okLow)
				{
					/*
					 * search the first activated marker m "almost" below p
					 */
					int first = firstActivatedMarkerIndex;
					max_weight = 0;
					max_index = -1;
					while (first != -1 && first <= lastActivatedMarkerIndex)
					{
						Fragment first_frag = Fragments.get(Markers.get(first).getOriginIndex());
							
						if ( first_frag.almostInferiorTo(frag_p) )
						{
							int overlap1 = frag_p.getX1() - first_frag.getY1();
							int overlap2 = frag_p.getX2() - first_frag.getY2();
							
							if ( overlap1 > 0)
								overlap1 = 0;
							if ( overlap2 > 0)
								overlap2 = 0;
							
							int chain_len = Markers.get(first).getWeight() + frag_p.getLength1() + frag_p.getLength2() + overlap1 +overlap2;
							
							if ( chain_len > max_weight )
							{
								max_weight = chain_len;
								max_index = first;
							}
						}
					
						first = Markers.get(first).getMarkerAfter();
					}
					
				}
				else
				{
					//up corner
					
					/*
					 * compute the best marker for the upper corner of the current fragment
					 */
					int max_weight_up = 0;
					int max_index_up = -1;
					
					int first = firstActivatedMarkerIndex;
					while (first != -1 && first <= lastActivatedMarkerIndex && first <= p.getIndexPoint_y())
					{
						
						Fragment first_frag = Fragments.get(Markers.get(first).getOriginIndex());
						
						
						int overlap1 = frag_p.getX1() - first_frag.getY1();
						int overlap2 = frag_p.getX2() - first_frag.getY2();
						
						if ( overlap1 > 0)
							overlap1 = 0;
						if ( overlap2 > 0)
							overlap2 = 0;
						
						int chain_len = Markers.get(first).getWeight() + frag_p.getLength1() + frag_p.getLength2() + overlap1 +overlap2;
							
						if ( chain_len > max_weight_up )
						{
							max_weight_up = chain_len;
							max_index_up = first;
						}
						
						////////////////
					 	marker_index_before_mp = first;
					 	
						first = Markers.get(first).getMarkerAfter();
					}
					
					
					max_index = max_index_up;
					max_weight = max_weight_up;	
				}
					
				
				if ( max_index == -1 )
				{
					m = null;
				}
				else
					m = Markers.get(max_index);
			}
			
			/*
			 * II "test whether p is a lower corner or an upper corner" phase
			 */
			
		
			if ( p.okLow )
			{
				/*
				 * if lower corner
				 */
				if ( m != null )
				{
					frag_p.setWeightChain(max_weight);
					frag_p.setLink(m.getOriginIndex());
				}
				else
				{
					frag_p.setWeightChain(frag_p.getLength1()+frag_p.getLength2());
					frag_p.setLink(-1);
				}
			}
			else
			{
				/*
				 * if upper corner
				 */
				boolean ok = false;
				
				if (m == null)
					ok = true;
				else
				if (frag_p.getWeightChain() > m.getWeight())
					ok = true;
				
				if (!ok)
				{
					/*
					 * search the corresponding "low" for which we are closing the "up" and delete it from the list
					 */		
					
					int k = 0;
					for (int idx: opened_lows)
					{
						Marker marker = Markers.get(idx);
						int frag_index = marker.getOriginIndex();
						
						if ( frag_index == mp.getOriginIndex())
						{
							opened_lows.remove(k);
							break;
						}
						k++;
					}
				}
				else
				{			
					/*
					 * activate marker mp
					 */
					mp.setMarked(true);
					mp.setWeight(frag_p.getWeightChain());
					
					nb_markers++;
					
					///////////
					if (marker_index_before_mp == -1)
					{
						mp.setMarkerBefore(-1);
						mp.setMarkerAfter(firstActivatedMarkerIndex);
						if (firstActivatedMarkerIndex != -1)
							Markers.get(firstActivatedMarkerIndex).setMarkerBefore(p.indexPoint_y);
					}
					else
					{
						mp.setMarkerBefore(marker_index_before_mp);
						mp.setMarkerAfter(Markers.get(marker_index_before_mp).getMarkerAfter());
						if ( Markers.get(marker_index_before_mp).getMarkerAfter() != -1 )
							Markers.get(Markers.get(marker_index_before_mp).getMarkerAfter()).setMarkerBefore(p.indexPoint_y);
						Markers.get(marker_index_before_mp).setMarkerAfter(p.getIndexPoint_y());
			
					}
								
					/*
					 * update firstActivatedMarkerIndex, lastActivatedMarkerIndex
					 */
					if (p.getIndexPoint_y() < firstActivatedMarkerIndex || firstActivatedMarkerIndex == -1 )
						firstActivatedMarkerIndex = p.getIndexPoint_y();
					if (p.getIndexPoint_y() > lastActivatedMarkerIndex || lastActivatedMarkerIndex == -1 )
						lastActivatedMarkerIndex = p.getIndexPoint_y();
					
					
					/*
					 * search the corresponding "low" for which we are closing the "up" and delete it from the list
					 * 
					 * for "lows" before, starting with the last one test whether frag_p isAlmostInferior to frag_low 
					 * and if the weight of frag_low > weight of frag_p with frag_low
					 */		
					
					int k = 0;
					int remove_index = -1;
				
					for (int idx : opened_lows)
					{
						Marker marker = Markers.get(idx);
						int frag_index = marker.getOriginIndex();
						Fragment frag = Fragments.get(frag_index);
						
						
						if ( frag_index == mp.getOriginIndex())
						{
							remove_index = k;
						}
						else	
						if ( frag_p.almostInferiorTo(frag))
						{	
							int overlap1 = frag.getX1() - frag_p.getY1();
							int overlap2 = frag.getX2() - frag_p.getY2();
								
							
							if ( overlap1 > 0)
								overlap1 = 0;
							if ( overlap2 > 0)
								overlap2 = 0;
							
							int chain_len = frag_p.getWeightChain() + frag.getLength1() + frag.getLength2() + overlap1 +overlap2;
							
							
							
							if ( chain_len > frag.getWeightChain() )
							{
								frag.setWeightChain(chain_len);
								frag.setLink(mp.getOriginIndex());
							}
						}	
							
						k++;
					}
					if ( remove_index != -1 )
						opened_lows.remove(remove_index);
				
					/*
					 * III "remove all mi higher and lighter than W(mp) - biggest accepted overlap for origin(mp)" phase
					 * this doesn't work (see example on paper)
					 * 
					 * remove all mi higher than mp that have either
					 * 		dif_mi_mp_ups>0 and dif_mi_mp_ups + overlap_mp <= overlap_mi and mi.getWeight() < frag_p.Link.getWeight()
					 * or
					 * 		dif_mi_mp_ups>0 and mi.getWeight() < mp.getWeight()
					 */
					int i = mp.getMarkerAfter();
					
					int limit_weight;
					if ( frag_p.getLink() != -1 )
						limit_weight = Fragments.get(frag_p.getLink()).getWeightChain();
					else
						limit_weight = 0;
					
					while (i != -1 && i<=lastActivatedMarkerIndex)
					{
						Marker mi = Markers.get(i);
						
						boolean test = false;
						
						Fragment frag_i = Fragments.get(mi.getOriginIndex());
							
						int overlap_mi = FragmentListBacktrack.percentage(frag_i.length2, frag_i.level) ;
						int overlap_mp = FragmentListBacktrack.percentage(frag_p.length2, frag_p.level) ;
						
						int dif_mi_mp_ups = frag_i.y2 - frag_p.y2;
						
						if ( dif_mi_mp_ups >0 && dif_mi_mp_ups + overlap_mp <= overlap_mi )
						{
							if ( mi.getWeight() < limit_weight )
								test = true;
						}
						else
							if (dif_mi_mp_ups >0 && mi.getWeight() < mp.getWeight() )
								test = true; 
							
						if ( test )
						{			
							nb_markers--;
			
							mi.setMarked(false);			
							
							if ( mi.getMarkerBefore() != -1 )
								Markers.get(mi.getMarkerBefore()).setMarkerAfter(mi.getMarkerAfter());
							
							if ( mi.getMarkerAfter() != -1 )
								Markers.get(mi.getMarkerAfter()).setMarkerBefore(mi.getMarkerBefore());
							
							if ( i== firstActivatedMarkerIndex )
								firstActivatedMarkerIndex = mi.getMarkerAfter();
							
							if ( i == lastActivatedMarkerIndex )
								lastActivatedMarkerIndex = mi.getMarkerBefore();
							
							i = mi.getMarkerAfter();
							
							mi.setMarkerAfter(-1);
							mi.setMarkerBefore(-1);
										
						}
						else
							i = mi.getMarkerAfter();
							
					}
				}
				
			}
		}
		/*try{
			BufferedWriter bf = new BufferedWriter (new FileWriter("markers", true));
			bf.append(nb_markers+" "+fragmentFile.getName()+"\n");
			bf.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}*/
		
		int i = firstActivatedMarkerIndex;
		int imax=-1, weightmax=-1;
		
		while ( i!=-1 && i<= lastActivatedMarkerIndex )
		{
			Marker m = Markers.get(i);
			if ( weightmax < m.getWeight() )
			{
				weightmax = m.getWeight();
				imax = i;
			}
			i = m.getMarkerAfter();
		}
		
		//imax = lastActivatedMarkerIndex;
		
		/*
		 * backtrack the best chain	
		 */

		i = Markers.get(imax).getOriginIndex();
		while ( i > -1 )
		{
			Chain.add(i);
			i = Fragments.get(i).getLink();
		}
		
		Collections.sort(Chain);

	}
	
	
	/*
	 * using the method of the sweepline, we try to find the best chain (accept overlaps)
	 */
	public void chainWithOverlaps()
	{
		/*
		 * we have all points sorted on x in Points
		 * we have non-activated markers sorted on y in Markers (for which we have access to their position in Points
		 * creating a marker for us means activating a marker in Markers
		 */
		
		int firstActivatedMarkerIndex = -1;
		int lastActivatedMarkerIndex = -1;
		
		
		int pace = -1;
		
		int nb_markers = 0;
		
		
		/*
		 * list of low points opened and not already closed
		 */
		LinkedList<Integer> opened_lows = new LinkedList<Integer>();

		for (Point p : Points)
		{
			pace++;
			int max_weight = 0;
			int max_index = -1;
			
			/*
			 * search the non-activated marker corresponding to p
			 */
			Marker mp = Markers.get(p.getIndexPoint_y());
			Fragment frag_p = Fragments.get(p.getIndexFragment());
			
		
			/*
			 * update the opened_lows list
			 */
			if (p.isOkLow())
			{
				/*
				 * add the newly opened lows to the list
				 */
				opened_lows.add(p.indexPoint_y);
			}
			
			/*
			 * I "search first activated marker almost below p" phase
			 */
				
			Marker m;
			int marker_index_before_mp = -1;
			
			if ( firstActivatedMarkerIndex == -1 )
			{
				m = null;
			}
			else
			{
				if (p.okLow)
				{
					/*
					 * search the first activated marker m "almost" below p
					 */
					int first = firstActivatedMarkerIndex;
					max_weight = 0;
					max_index = -1;							
					while (first != -1 && first <= lastActivatedMarkerIndex)
					{
						Fragment first_frag = Fragments.get(Markers.get(first).getOriginIndex());
							
						if ( first_frag.almostInferiorToUpgraded(frag_p) )
						{
							int overlap1 = frag_p.getX1() - first_frag.getY1();
							int overlap2 = frag_p.getX2() - first_frag.getY2();
							
							if ( overlap1 > 0)
								overlap1 = 0;
							if ( overlap2 > 0)
								overlap2 = 0;
							
							int chain_len = Markers.get(first).getWeight() + frag_p.getLength1() + frag_p.getLength2() + overlap1 +overlap2;
							
							if ( chain_len > max_weight )
							{
								max_weight = chain_len;
								max_index = first;
							}
						}
					
						first = Markers.get(first).getMarkerAfter();
					}
					
				}
				else
				{
					//up corner
					
					/*
					 * compute the best marker for the upper corner of the current fragment
					 */
					int max_weight_up = 0;
					int max_index_up = -1;
					
					int first = firstActivatedMarkerIndex;
					while (first != -1 && first <= lastActivatedMarkerIndex && first <= p.getIndexPoint_y())
					{
						
						Fragment first_frag = Fragments.get(Markers.get(first).getOriginIndex());
						
						
						int overlap1 = frag_p.getX1() - first_frag.getY1();
						int overlap2 = frag_p.getX2() - first_frag.getY2();
						
						if ( overlap1 > 0)
							overlap1 = 0;
						if ( overlap2 > 0)
							overlap2 = 0;
						
						int chain_len = Markers.get(first).getWeight() + frag_p.getLength1() + frag_p.getLength2() + overlap1 +overlap2;
							
						if ( chain_len > max_weight_up )
						{
							max_weight_up = chain_len;
							max_index_up = first;
						}
						
						////////////////
					 	marker_index_before_mp = first;
					 	
						first = Markers.get(first).getMarkerAfter();
					}
					
					
					max_index = max_index_up;
					max_weight = max_weight_up;	
				}
					
				
				if ( max_index == -1 )
				{
					m = null;
				}
				else
					m = Markers.get(max_index);
			}
			
			/*
			 * II "test whether p is a lower corner or an upper corner" phase
			 */
			
		
			if ( p.okLow )
			{
				/*
				 * if lower corner
				 */
				if ( m != null )
				{
					frag_p.setWeightChain(max_weight);
					frag_p.setLink(m.getOriginIndex());
				}
				else
				{
					frag_p.setWeightChain(frag_p.getLength1()+frag_p.getLength2());
					frag_p.setLink(-1);
				}
			}
			else
			{
				/*
				 * if upper corner
				 */
				boolean ok = false;
				
				if (m == null)
					ok = true;
				else
				if (frag_p.getWeightChain() > m.getWeight())
					ok = true;
				
				/*
				 * search the corresponding "low" for which we are closing the "up" and delete it from the list
				 * 
				 * for "lows" before, starting with the last one test whether frag_p isAlmostInferior to frag_low 
				 * and if the weight of frag_low > weight of frag_p with frag_low
				 */		
				
				int k = 0;
				int remove_index = -1;
			
				for (int idx : opened_lows)
				{
					Marker marker = Markers.get(idx);
					int frag_index = marker.getOriginIndex();
					Fragment frag = Fragments.get(frag_index);
					
					
					if ( frag_index == mp.getOriginIndex())
					{
						remove_index = k;
					}
					else	
					if ( frag_p.almostInferiorToUpgraded(frag))
					{	
						int overlap1 = frag.getX1() - frag_p.getY1();
						int overlap2 = frag.getX2() - frag_p.getY2();
							
						
						if ( overlap1 > 0)
							overlap1 = 0;
						if ( overlap2 > 0)
							overlap2 = 0;
						
						int chain_len = frag_p.getWeightChain() + frag.getLength1() + frag.getLength2() + overlap1 +overlap2;
						
						
						
						if ( chain_len > frag.getWeightChain() )
						{
							frag.setWeightChain(chain_len);
							frag.setLink(mp.getOriginIndex());
						}
					}	
						
					k++;
				}
				if ( remove_index != -1 )
					opened_lows.remove(remove_index);
				
				
				if ( ok )
				{
					/*
					 * activate marker mp
					 */
					mp.setMarked(true);
					mp.setWeight(frag_p.getWeightChain());
					
					nb_markers++;
					
					///////////
					if (marker_index_before_mp == -1)
					{
						mp.setMarkerBefore(-1);
						mp.setMarkerAfter(firstActivatedMarkerIndex);
						if (firstActivatedMarkerIndex != -1)
							Markers.get(firstActivatedMarkerIndex).setMarkerBefore(p.indexPoint_y);
					}
					else
					{
						mp.setMarkerBefore(marker_index_before_mp);
						mp.setMarkerAfter(Markers.get(marker_index_before_mp).getMarkerAfter());
						if ( Markers.get(marker_index_before_mp).getMarkerAfter() != -1 )
							Markers.get(Markers.get(marker_index_before_mp).getMarkerAfter()).setMarkerBefore(p.indexPoint_y);
						Markers.get(marker_index_before_mp).setMarkerAfter(p.getIndexPoint_y());
			
					}
								
					/*
					 * update firstActivatedMarkerIndex, lastActivatedMarkerIndex
					 */
					if (p.getIndexPoint_y() < firstActivatedMarkerIndex || firstActivatedMarkerIndex == -1 )
						firstActivatedMarkerIndex = p.getIndexPoint_y();
					if (p.getIndexPoint_y() > lastActivatedMarkerIndex || lastActivatedMarkerIndex == -1 )
						lastActivatedMarkerIndex = p.getIndexPoint_y();
					
					
				
					/*
					 * III "remove all mi higher 
					 * with W(mi) lighter than W(mp)
					 * and height(frag_i) < height(frag_p)
					 */
					int i = mp.getMarkerAfter();
					
					int limit_weight = mp.getWeight();
					
					while (i != -1 && i<=lastActivatedMarkerIndex)
					{
						Marker mi = Markers.get(i);
						
						Fragment frag_i = Fragments.get(mi.getOriginIndex());
								
						int dif_mi_mp_ups = frag_i.y2 - frag_p.y2;

							
						if ( dif_mi_mp_ups >= 0 && frag_i.length2 < frag_p.length2 && mi.getWeight() < limit_weight )
						{			
							nb_markers--;
			
							mi.setMarked(false);			
							
							if ( mi.getMarkerBefore() != -1 )
								Markers.get(mi.getMarkerBefore()).setMarkerAfter(mi.getMarkerAfter());
							
							if ( mi.getMarkerAfter() != -1 )
								Markers.get(mi.getMarkerAfter()).setMarkerBefore(mi.getMarkerBefore());
							
							if ( i== firstActivatedMarkerIndex )
								firstActivatedMarkerIndex = mi.getMarkerAfter();
							
							if ( i == lastActivatedMarkerIndex )
								lastActivatedMarkerIndex = mi.getMarkerBefore();
							
							i = mi.getMarkerAfter();
							
							mi.setMarkerAfter(-1);
							mi.setMarkerBefore(-1);
										
						}
						else
							i = mi.getMarkerAfter();
					}		
				}
				
			}
		}
		/*try{
			BufferedWriter bf = new BufferedWriter (new FileWriter("markers", true));
			bf.append(nb_markers+" "+fragmentFile.getName()+"\n");
			bf.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}*/
		
		int i = firstActivatedMarkerIndex;
		int imax=-1, weightmax=-1;
		
		while ( i!=-1 && i<= lastActivatedMarkerIndex )
		{
			Marker m = Markers.get(i);
			if ( weightmax < m.getWeight() )
			{
				weightmax = m.getWeight();
				imax = i;
			}
			i = m.getMarkerAfter();
		}
		
		//imax = lastActivatedMarkerIndex;
		
		/*
		 * backtrack the best chain	
		 */

		i = Markers.get(imax).getOriginIndex();
		while ( i > -1 )
		{
			Chain.add(i);
			i = Fragments.get(i).getLink();
		}
		
		Collections.sort(Chain);

	}
}

