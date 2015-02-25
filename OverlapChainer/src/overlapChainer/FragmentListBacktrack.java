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
import java.util.Collections;
import java.util.List;
import java.util.Vector;

/**
 * @author Raluca Uricaru
 *
 */

public class FragmentListBacktrack extends FragmentListHierarchical {

	public static final int NO_MAX_FRAG_BKTR = 30;
	public static final boolean fix_Nb_Frag_Bktr = true;
	
	ArrayList<Integer> Frag_for_backtr;
	ArrayList<Integer> Sol_backtr;
	int couv1, couv2;
	
	public FragmentListBacktrack(File fragmentFile) {
		super(fragmentFile);
	}
	
	/*
	 * in case of java heap space error
	 */
	public FragmentListBacktrack(File fragmentFile, boolean ok) {
		super(fragmentFile, ok);
	}
	
	public void init(){
		
		if ( !fix_Nb_Frag_Bktr )
		{
			/*
			 * compute the number of levels and their corresponding limits
			 */
			setLevels();
			
		}
		else
		{
			/*
			 * build the list of Fragments that will be considered
			 * for the backtracking
			 */
			Frag_for_backtr = new ArrayList<Integer>();
			
			for ( int i=0; i<Fragments.size(); i++)
			{
				Fragment f = Fragments.get(i);
				if ( Frag_for_backtr == null || Frag_for_backtr.size() == 0 )
					Frag_for_backtr.add(new Integer(i));
				else
				{
					if ( f.getWeight() > Fragments.get(Frag_for_backtr.get(0)).getWeight()
							||
						 f.getWeight() == Fragments.get(Frag_for_backtr.get(0)).getWeight()	&&
						 f.getLength1()+f.getLength2() > Fragments.get(Frag_for_backtr.get(0)).getLength1() + Fragments.get(Frag_for_backtr.get(0)).getLength2())
					{
						int j=1;
						while ( j<Frag_for_backtr.size() && (f.getWeight() > Fragments.get(Frag_for_backtr.get(j)).getWeight()
							||
						 f.getWeight() == Fragments.get(Frag_for_backtr.get(j)).getWeight()	&&
						 f.getLength1()+f.getLength2() > Fragments.get(Frag_for_backtr.get(j)).getLength1() + Fragments.get(Frag_for_backtr.get(j)).getLength2()))
							j++;
					
						Frag_for_backtr.add(j, new Integer(i));
						if (Frag_for_backtr.size() > NO_MAX_FRAG_BKTR )
							Frag_for_backtr.remove(0);
					}
				}
			}
			
			setLevels_for_fix_Nb_Frag_Bktr(Fragments.get(Frag_for_backtr.get(0)).getWeight());
		}
		
		
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
	 * compute the number of levels and their corresponding limits
	 */
	public void setLevels_for_fix_Nb_Frag_Bktr(int limit_level0) {
		
		int lev = 0;
		noRecLevels = 0;
		limitRecLevels_list = new ArrayList<Integer>();
		limitPercentage_list = new ArrayList<Integer>();
		
		limitRecLevels_list.add(new Integer(limit_level0));
		if (limit_level0 > 10000)
			limitPercentage_list.add(new Integer(20));
		else
			limitPercentage_list.add(new Integer(10));
		noRecLevels++;
		
		int fWMax = limit_level0-1;
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
	 * compute the number of levels and their corresponding limits
	 */
	@Override
	public void setLevels() {
		
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

	
	private boolean verifySolution (ArrayList<Integer> result, int [] couv){
		
		int couv1_aux = couv[0];
		int couv2_aux = couv[1];
		Fragment fragment_before, fragment_last;
		int dif1,dif2;
		
		boolean ok = true;
		
		if ( result.size() > 1 )
		{
			fragment_last = Fragments.get(Frag_for_backtr.get(result.get(result.size()-1)));
			fragment_before = Fragments.get(Frag_for_backtr.get(result.get(result.size()-2)));
			
			if ( fragment_before.getX2() > fragment_last.getX2() )
				ok = false;
			else
			{
				dif1 = fragment_last.getX1() - fragment_before.getY1();
				dif2 = fragment_last.getX2() - fragment_before.getY2();
				
				int lev;
				if ( fragment_last.getLevel() < fragment_before.getLevel() )
					lev = fragment_last.getLevel();
				else
					lev = fragment_before.getLevel();
			
				if ( dif1 < 0 )
				{
					couv1_aux += dif1;
					
					/*
					 * verify the overlap with the pre-fragment in the chain
					 */
					if ( !(Math.abs(dif1) <= percentage(fragment_before.getLength1(), lev) && 
						 Math.abs(dif1) <= percentage(fragment_last.getLength1(), lev)))
						ok = false;
				}
				
				if ( dif2 < 0 ) 
				{
					couv2_aux += dif2;
					
					/*
					 * verify the overlap with the pre-fragment in the chain
					 */
					if ( !(Math.abs(dif2) <= percentage(fragment_before.getLength2(), lev) && 
							Math.abs(dif2) <= percentage(fragment_last.getLength2(), lev)))
						ok = false;
				}
				
				couv1_aux += fragment_last.getLength1();
				couv2_aux += fragment_last.getLength2();
			}
		}
		else
		{
			fragment_last = Fragments.get(Frag_for_backtr.get(result.get(result.size()-1)));
			couv1_aux = fragment_last.getLength1();
			couv2_aux = fragment_last.getLength2();
		}
		
		if ( ok )
		{
			couv[0] = couv1_aux;
			couv[1] = couv2_aux;
		}
		

		if ( ok && couv1_aux + couv2_aux > couv1+couv2 )
		{
			couv1 = couv1_aux;
			couv2 = couv2_aux;
			
			Sol_backtr.clear();
			Sol_backtr.addAll(result);
			return true;
		}
		
		return ok;
	}

	
	/*
	 *  p←1;
		st[p] ← 0;
		cat timp p>0 executa
		  inceput
		     daca  <mai exista valori neincercate pe nivelul p> atunci
		           inceput
		             st[p] ← <o noua valoare din multimea solutiilor posibile>
		             daca valid (p) returneaza true atunci
		               daca <solutia este finala> atunci
		                  apel tipar (p)
		               altfel
		                 inceput
		                   p ← p+1;
		                   st[p] ← 0;
		                 sfarsit
		             sfarsit
		             altfel
		              p ←p-1;
		  sfarsit

	 */
	private void try_backtr_solution ()
	{
		ArrayList<Integer> result = new ArrayList<Integer>();
		int[] couv = {0,0};
		
		int index = 0;
		result.add(-1);
		int val;
		
		while (index > -1)
		{
			
			val = result.get(index);
			if ( val < Frag_for_backtr.size()-1 )
			{
				val++;
				
				result.set(index, val);
				
				if ( verifySolution(result, couv) )
				{
					/*TODO Delete this in case it doesn't work
					 * if all fragments correspond, get out directly
					 */
					if ( result.size() == Frag_for_backtr.size() )
						return;
					
					/*System.out.print("oui ");
					for ( int i=0; i<result.size(); i++)
						System.out.print(result.get(i)+" ");
					System.out.println();*/
					
					index++;
					result.add(val);	
				}
				/*else
				{
					System.out.print("non ");
					for ( int i=0; i<result.size(); i++)
						System.out.print(result.get(i)+" ");
					System.out.println();
				}*/
			}
			else
			{
				result.remove(index);
				index--;
				
				if ( index > 0 )
				{
					//coverage up to date
					Fragment f2 = Fragments.get(Frag_for_backtr.get(result.get(index)));
					Fragment f1 = Fragments.get(Frag_for_backtr.get(result.get(index-1)));
					int dif1 = f2.getX1() - f1.getY1();
					int dif2 = f2.getX2() - f1.getY2();
					if ( dif1 < 0 )
						couv[0] -= dif1;
					if ( dif2 < 0 )
						couv[1] -= dif2;
					couv[0] -= f2.getLength1();
					couv[1] -= f2.getLength2();
				}
				else
				if ( index == 0)
				{
					couv[0] = 0;
					couv[1] = 0;
				}		
			}
		}
	}

	public int buildChainWithBigWeightFragments_backtr(){
		
		int lev = 0;
		int noFrag = 0;
		
		if ( !fix_Nb_Frag_Bktr )
		{		
			while (noFrag < NO_MAX_FRAG_BKTR && lev < noRecLevels ) 
			{
				noFrag += this.noFragOnLevel.get(lev);
				lev++;
			}
			
			lev--;
			if ( noFrag > NO_MAX_FRAG_BKTR )
			{
				noFrag -= this.noFragOnLevel.get(lev);
				lev--;
			}
			
			if ( Constants.DEBUG )
			{
				System.out.println("Limit level " + limitRecLevels_list.get(lev));
				System.out.println("Level " + lev + " from " + noRecLevels);
				System.out.println("Fragments for backtracking " + noFrag);
			}
			
			if ( lev >= 0 )
			{
				/*
				 * build the list of Fragments that will be considered
				 * for the backtracking
				 */
				Frag_for_backtr = new ArrayList<Integer>();
				if ( noFragments <= Constants.NO_MAX_FRAG )
				{
					for (int i=0; i<=lev; i++)
						Frag_for_backtr.addAll(this.Fragments_k.get(i));
				}
				else
				{
					for (int index = 0; index<Fragments.size(); index++ )
					{
						Fragment f = Fragments.get(index);
						if ( f.getLevel() <= lev )
							Frag_for_backtr.add(new Integer(index));				
					}
				}
			}
		}
			
		
		if ( lev >= 0 )
		{	
			/*
			System.out.println(Frag_for_backtr.size());
			for ( int i=0; i<NO_MAX_FRAG_BKTR; i++)
				System.out.println(Fragments.get(Frag_for_backtr.get(i)) + " " + Fragments.get(Frag_for_backtr.get(i)).getWeight());
			*/
			
			Collections.sort(Frag_for_backtr);
			
			Sol_backtr = new ArrayList<Integer>();
			couv1 = couv2 = 0;
				
			try_backtr_solution();
			
			if ( Constants.DEBUG )
				System.out.println(couv1 + " " + couv2);
			
			for ( int i=0; i< Sol_backtr.size(); i++)
			{
				Fragment f = Fragments.get(Frag_for_backtr.get(Sol_backtr.get(i)));
				Chain.add(Frag_for_backtr.get(Sol_backtr.get(i)));
				if ( Constants.DEBUG ) System.out.println(f);
			}
			
			if ( Constants.DEBUG ) System.out.println();
			
		}
		
		return lev+1;
	}
}
