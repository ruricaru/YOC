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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Vector;
import java.util.regex.Pattern;

/**
 * @author Raluca Uricaru
 *
 */

public abstract class AbstractFragmentListNotRearranged extends AbstractFragmentList{
	protected static final int NO_MAX_FRAG_READ = 800000;

	protected Vector<List<Integer>> Fragments_k; 
	protected Vector<Integer> noFragOnLevel;

	/*
	 * if we deal with problems due to circularity,
	 * we rearrange the initial list of fragments
	 * we need to keep the leftShift and the rightShift values
	 * on the second genome
	 */
	protected int leftShift;
	protected int rightShift;

	public static int noRecLevels;
	public static ArrayList<Integer> limitRecLevels_list;
	public static ArrayList<Integer> limitPercentage_list;

	public static int noStepsOnLevel = 5;
	public static Integer[] limitShiftOnStep;

	public Vector<List<Integer>> getFragments_k() {
		return Fragments_k;
	}
	public List<Integer> getFragments_k( int k ) {	
		return Fragments_k.get(k);
	}
	public Vector<Integer> getNoFragOnLevel() {
		return noFragOnLevel;
	}
	public void setNoFragOnLevel(Vector<Integer> noFragOnLevel) {
		this.noFragOnLevel = noFragOnLevel;
	}
	public int getLeftShift() {
		return leftShift;
	}
	public void setLeftShift(int leftShift) {
		this.leftShift = leftShift;
	}
	public int getRightShift() {
		return rightShift;
	}
	public void setRightShift(int rightShift) {
		this.rightShift = rightShift;
	}	

	protected int No_Max_Frag_Read()
	{
		return NO_MAX_FRAG_READ;
	}

	public AbstractFragmentListNotRearranged (){

	}


	public AbstractFragmentListNotRearranged (File fragmentFile) {
		super(fragmentFile);

		this.leftShift = 0;
		this.rightShift = 0;

		/*
		 * sort fragments on x position in the first genome
		 */
		sortFragments_x1();
	}	

	/*
	 * in case of java heap space error
	 */
	public AbstractFragmentListNotRearranged (File fragmentFile, boolean ok) {
		super(fragmentFile, ok);

		this.leftShift = 0;
		this.rightShift = 0;

		/*
		 * sort fragments on x position in the first genome
		 */
		sortFragments_x1();
	}	

	/*
	 * sort fragments on x in the first genome
	 */

	@SuppressWarnings("unchecked")
	public void sortFragments_x1 ()
	{
		Collections.sort(this.Fragments);
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

		for (k=0; k<noRecLevels; k++)
			this.getFragments_k(k).clear();

		for ( int index=0; index<Fragments.size(); index++ ){
			fragment = this.Fragments.get(index);
			/*
			 * we try to find the proper recursion level for the fragment
			 */
			for (k=0; k<noRecLevels; k++)
			{
				size1 = limitRecLevels_list.get(k).intValue();
				size2 = 10000000;
				if ( k > 0 )
					size2 = limitRecLevels_list.get(k-1).intValue();

				if ( fragment.getWeight()>= size1 && fragment.getWeight()<size2 )
				{
					if ( noFragments <= Constants.NO_MAX_FRAG )
						this.getFragments_k(k).add(new Integer(index));
					this.Fragments.get(index).setLevel(k);
					this.noFragOnLevel.set(k, this.noFragOnLevel.get(k).intValue()+1);
					break;
				}
			}
		}
	}

	public static int percentage (int x, int k){
		int overlap = limitPercentage_list.get(k).intValue();

		// constant and not ratio
		if ( overlap >= 100 )
			return overlap/100;
		else
			return (int) ((double) x/100 * overlap );
	}

	/*
	 * recursive level k
	 * add to the chain fragments that are already in the Fragments_k list
	 * that don't overlap more than the accepted percentage
	 * between fragment in position "begin" and fragment in position "end" (global list of fragments)
	 * index_chain = the position in the chain after which we insert 
	 * step = on every level we execute several steps of fragment filling 
	 * 		  on every step we try to add some more fragments from the current level
	 * 		  from one step to another, we accept some more "shifting"
	 */
	public void addToChain_Fragments_k (int begin, int end, int k, int index_chain, int step)
	{
		int index;

		Fragment fragment = new Fragment();
		Fragment before = new Fragment();
		Fragment after = new Fragment();

		List<Integer> fragmentList = new ArrayList<Integer>();
		int fragmentListSize = 0;

		int aux_k = k;
		if ( k == noRecLevels )
		{
			//on the last level we allow any of the fragments
			aux_k = k-1;
			fragmentListSize = Fragments.size();
		}
		else
		{
			fragmentList = this.getFragments_k(k);
			if ( this.getFragments_k(k).isEmpty() )
				fragmentListSize = 0;
			else
				fragmentListSize = fragmentList.size();
		}

		int dif1, dif2;
		int no_f=0, no_f_c=0;

		/*
		 * k-lists were not created
		 * due to the big number of fragments
		 */
		if ( this.noFragments > Constants.NO_MAX_FRAG )
		{
			if ( begin == -1 )
				begin = 0;
			if ( end == -1 )
				end = this.noFragments-1;

			for (int i=begin; i<=end; i++)
			{
				fragment = this.Fragments.get(i);

				/*
				 * test the level of the fragment
				 */
				if ( k != noRecLevels )
					if ( fragment.getLevel() != k )
						continue;

				no_f_c++;

				/*
				 * test the shifting condition on this step
				 */
				if ( k != noRecLevels )
				{
					int shift = limitShiftOnStep[step].intValue();

					/*
					 * x positions relative to the current window
					 * window = formed by 2 fragments of lower level
					 */
					int x1_relative=fragment.getX1(), x2_relative=fragment.getX2();

					if ( this.Chain.size() == 0 || index_chain == -1 || k==0 )
					{
						x1_relative = fragment.getX1();
						x2_relative = fragment.getX2();
					}
					else
					{
						int previous_level_fragment = index_chain;
						while ( this.Fragments.get(Chain.get(previous_level_fragment)).getLevel() >= k && previous_level_fragment>0)
							previous_level_fragment--;			
						if ( this.Fragments.get(Chain.get(previous_level_fragment)).getLevel() < k )
						{
							before = this.Fragments.get(Chain.get(previous_level_fragment));
							x1_relative = fragment.getX1() - before.getX1();
							x2_relative = fragment.getX2() - before.getX2();
						}
						else
						{
							x1_relative = fragment.getX1();
							x2_relative = fragment.getX2();
						}
					}

					if ( shift <= Math.abs(x1_relative-x2_relative) )
						continue;
				}

				if ( this.Chain.size() == 0 )
				{
					this.Chain.add(new Integer(i));
					no_f++;
					index_chain++;
				}
				else
					if ( index_chain > -1 )
					{
						before = this.Fragments.get(Chain.get(index_chain));
						dif1 = fragment.getX1() - before.getY1();
						dif2 = fragment.getX2() - before.getY2();

						/*
						 * verify the overlap with the pre-fragment in the chain
						 */
						if ( ((dif1>=0) || (dif1<0 && 
								Math.abs(dif1) <= percentage(before.getLength1(), aux_k) && 
								Math.abs(dif1) <= percentage(fragment.getLength1(), aux_k))) &&
								((dif2>=0) || (dif2<0 && 
										Math.abs(dif2) <= percentage(before.getLength2(), aux_k) && 
										Math.abs(dif2) <= percentage(fragment.getLength2(), aux_k))) )
						{
							/*
							 * verify the overlap with the post-fragment in the chain
							 */
							if (this.Chain.size() > 1 && index_chain<this.Chain.size()-1)
							{
								after = this.Fragments.get(Chain.get(index_chain+1));
								dif1 = after.getX1() - fragment.getY1();
								dif2 = after.getX2() - fragment.getY2();

								if ( ((dif1>=0) || (dif1<0 && 
										Math.abs(dif1) <= percentage(after.getLength1(), aux_k) && 
										Math.abs(dif1) <= percentage(fragment.getLength1(), aux_k))) &&
										((dif2>=0) || (dif2<0 && 
												Math.abs(dif2) <= percentage(after.getLength2(), aux_k) && 
												Math.abs(dif2) <= percentage(fragment.getLength2(), aux_k))) )
								{
									/*
									 * add the fragment to the chain in the index_chain+1 position
									 * translate all the other fragments in the chain
									 */
									this.Chain.add(index_chain+1, new Integer(i));
									index_chain++;
									no_f++;
								}
								else
								{
									if (Constants.DEBUG)
									{
										System.out.print(fragment.getX1() + "\t" + fragment.getY1() + "\t\t");
										System.out.print(fragment.getX2() + "\t" + fragment.getY2() + "\n");
									}
								}
							}
							else
							{
								this.Chain.add(new Integer(i));
								index_chain++;
								no_f++;
							}
						}
						else
						{
							if (Constants.DEBUG)
							{
								System.out.print(fragment.getX1() + "\t" + fragment.getY1() + "\t\t");
								System.out.print(fragment.getX2() + "\t" + fragment.getY2() + "\n");
								System.out.println(dif1 + " " + dif2);
								System.out.println(fragment.getLength1() + " " + fragment.getLength2());
								System.out.println(before.getLength1() + " " + before.getLength2());
							}
						}
					}
					else
					{
						//index_chain == -1
						after = this.Fragments.get(Chain.get(0));
						dif1 = after.getX1() - fragment.getY1();
						dif2 = after.getX2() - fragment.getY2();

						if ( ((dif1>=0) || (dif1<0 && 
								Math.abs(dif1) <= percentage(after.getLength1(), aux_k) && 
								Math.abs(dif1) <= percentage(fragment.getLength1(), aux_k))) &&
								((dif2>=0) || (dif2<0 && 
										Math.abs(dif2) <= percentage(after.getLength2(), aux_k) && 
										Math.abs(dif2) <= percentage(fragment.getLength2(), aux_k))) )
						{
							/*
							 * add the fragment to the chain in the first position
							 * translate all the other fragments in the chain
							 */
							this.Chain.add(0, new Integer(i));
							index_chain++;
							no_f++;
						}
						else
						{
							if (Constants.DEBUG)
							{
								System.out.print(fragment.getX1() + "\t" + fragment.getY1() + "\t\t");
								System.out.print(fragment.getX2() + "\t" + fragment.getY2() + "\n");
							}
						}
					}
			}
		}
		else
		{
			//start adding fragments
			for (int i=0; i<fragmentListSize; i++){

				if ( k == noRecLevels )
					index = i;
				else
					index = fragmentList.get(i).intValue();
				/*
				 * get the current fragment
				 */
				fragment = this.Fragments.get(index);

				/*
				 * we are out of the interval
				 */
				if (begin!=-1 && end!=-1 && ( index < begin || index > end ))
					continue;
				else
					if ( begin == -1 && index > end )
						continue;
					else
						if ( end == -1 && index < begin )
							continue;

				no_f_c++;

				/*
				 * test the shifting condition on this step
				 */
				if ( k < noRecLevels )
				{
					int shift = limitShiftOnStep[step].intValue();

					/*
					 * x positions relative to the current window
					 * window = formed by 2 fragments of lower level
					 */
					int x1_relative=fragment.getX1(), x2_relative=fragment.getX2();

					if ( this.Chain.size() == 0 || index_chain == -1 || k==0 )
					{
						x1_relative = fragment.getX1();
						x2_relative = fragment.getX2();
					}
					else
					{
						int previous_level_fragment = index_chain;
						while ( this.Fragments.get(Chain.get(previous_level_fragment)).getLevel() >= k && previous_level_fragment>0)
							previous_level_fragment--;			
						if ( this.Fragments.get(Chain.get(previous_level_fragment)).getLevel() < k )
						{
							before = this.Fragments.get(Chain.get(previous_level_fragment));
							x1_relative = fragment.getX1() - before.getX1();
							x2_relative = fragment.getX2() - before.getX2();
						}
						else
						{
							x1_relative = fragment.getX1();
							x2_relative = fragment.getX2();
						}
					}

					if ( shift <= Math.abs(x1_relative-x2_relative) )
						continue;
				}

				if ( this.Chain.size() == 0 )
				{
					this.Chain.add(new Integer(index));
					no_f++;
					index_chain++;
				}
				else
					if ( index_chain > -1 )
					{
						before = this.Fragments.get(Chain.get(index_chain));
						dif1 = fragment.getX1() - before.getY1();
						dif2 = fragment.getX2() - before.getY2();

						/*
						 * verify the overlap with the pre-fragment in the chain
						 */
						if ( ((dif1>=0) || (dif1<0 && 
								Math.abs(dif1) <= percentage(before.getLength1(), aux_k) && 
								Math.abs(dif1) <= percentage(fragment.getLength1(), aux_k))) &&
								((dif2>=0) || (dif2<0 && 
										Math.abs(dif2) <= percentage(before.getLength2(), aux_k) && 
										Math.abs(dif2) <= percentage(fragment.getLength2(), aux_k))) )
						{
							/*
							 * verify the overlap with the post-fragment in the chain
							 */
							if (this.Chain.size() > 1 && index_chain<this.Chain.size()-1)
							{
								after = this.Fragments.get(Chain.get(index_chain+1));
								dif1 = after.getX1() - fragment.getY1();
								dif2 = after.getX2() - fragment.getY2();

								if ( ((dif1>=0) || (dif1<0 && 
										Math.abs(dif1) <= percentage(after.getLength1(), aux_k) && 
										Math.abs(dif1) <= percentage(fragment.getLength1(), aux_k))) &&
										((dif2>=0) || (dif2<0 && 
												Math.abs(dif2) <= percentage(after.getLength2(), aux_k) && 
												Math.abs(dif2) <= percentage(fragment.getLength2(), aux_k))) )
								{
									/*
									 * add the fragment to the chain in the index_chain+1 position
									 * translate all the other fragments in the chain
									 */
									this.Chain.add(index_chain+1, new Integer(index));
									index_chain++;
									no_f++;
								}
								else
								{
									if (Constants.DEBUG)
									{
										System.out.print(fragment.getX1() + "\t" + fragment.getY1() + "\t\t");
										System.out.print(fragment.getX2() + "\t" + fragment.getY2() + "\n");
									}
								}
							}
							else
							{
								this.Chain.add(new Integer(index));
								index_chain++;
								no_f++;
							}
						}
						else
						{
							if (Constants.DEBUG)
							{
								System.out.print(fragment.getX1() + "\t" + fragment.getY1() + "\t\t");
								System.out.print(fragment.getX2() + "\t" + fragment.getY2() + "\n");
								System.out.println(dif1 + " " + dif2);
								System.out.println(fragment.getLength1() + " " + fragment.getLength2());
								System.out.println(before.getLength1() + " " + before.getLength2());
							}
						}
					}
					else
					{
						//index_chain == -1
						after = this.Fragments.get(Chain.get(0));
						dif1 = after.getX1() - fragment.getY1();
						dif2 = after.getX2() - fragment.getY2();

						if ( ((dif1>=0) || (dif1<0 && 
								Math.abs(dif1) <= percentage(after.getLength1(), aux_k) && 
								Math.abs(dif1) <= percentage(fragment.getLength1(), aux_k))) &&
								((dif2>=0) || (dif2<0 && 
										Math.abs(dif2) <= percentage(after.getLength2(), aux_k) && 
										Math.abs(dif2) <= percentage(fragment.getLength2(), aux_k))) )
						{
							/*
							 * add the fragment to the chain in the first position
							 * translate all the other fragments in the chain
							 */
							this.Chain.add(0, new Integer(index));
							index_chain++;
							no_f++;
						}
						else
						{
							if (Constants.DEBUG)
							{
								System.out.print(fragment.getX1() + "\t" + fragment.getY1() + "\t\t");
								System.out.print(fragment.getX2() + "\t" + fragment.getY2() + "\n");
							}
						}
					}
			}
		}
	}


	public void circulariseFastaFile (String f1, String f2, int leftShift, int rightShift )
	{	
		//starting from the f2 file, obtain a circularised file called 

		File fasta1 = new File(f1);
		File fasta2 = new File(f2);
		String new_fasta_file_name = fasta1.getName().substring(0, fasta1.getName().lastIndexOf('.')).concat("_");
		new_fasta_file_name = new_fasta_file_name.concat(fasta2.getName().substring(0, fasta2.getName().lastIndexOf('.')));
		String new_fasta_file_path = fasta2.getAbsolutePath().substring(0, fasta2.getAbsolutePath().lastIndexOf('/')+1);

		//TODO
		/*
		 * we need to produce the circularised fasta file for every evalue
		 */
		String fname = this.fragmentFile.getName();
		String eval_s = "";
		if ( fname.lastIndexOf("mp") != -1)
		{
			fname = fname.substring(0, fname.lastIndexOf("mp")-1);
			eval_s = fname.substring(fname.lastIndexOf('_')+1, fname.length());
		}
		else
			eval_s = fname.substring(fname.lastIndexOf('_')+1, fname.lastIndexOf('.'));

		double eval = new Double (eval_s).doubleValue();
		new_fasta_file_name = new_fasta_file_name.concat("_"+eval_s+".fa");

		File new_fasta_file = new File (new_fasta_file_path.concat(new_fasta_file_name));

		try{
			BufferedWriter bf = new BufferedWriter(new FileWriter(Constants.circular_genomes, true));
			bf.append(fragmentFile.getName()+"\n");
			bf.close();

			bf = new BufferedWriter(new FileWriter(new_fasta_file));
			BufferedReader br = new BufferedReader(new FileReader(f2));

			String line = br.readLine();
			bf.append(line+"\n");
			int noChars = 0;

			BufferedWriter end = new BufferedWriter(new FileWriter("end.txt"));


			while ((line = br.readLine())!=null)
			{
				noChars += line.length();
				if ( noChars <= leftShift )
				{
					end.append(line+"\n");
				}

				if ( noChars > leftShift )
				{
					if ( noChars - line.length() > leftShift )
						bf.append(line+"\n");
					else
					{
						end.append(line.substring(0, leftShift%line.length()));
						bf.append(line.substring(leftShift%line.length())+"\n");
					}
				}
			}

			end.close();
			BufferedReader end_br = new BufferedReader(new FileReader("end.txt"));
			while ((line = end_br.readLine())!=null)
			{
				bf.append(line+"\n");
			}
			end_br.close(); bf.close(); br.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
			System.exit(0);
		}
	}



	public String Chain_toString(){
		String output = new String("");
		Fragment fragment, after;
		int coverage1=0, coverage2=0;
		int dif1, dif2;
		int original_x2, original_y2;

		Integer pozI;
		fragment = this.Fragments.get(Chain.get(Chain.size()-1));

		// if circular
		if ( fragment.getMovingDirection() != 0 && Constants.fasta1 != null && 
				Constants.fasta2 != null && Constants.circular_genomes != null)
		{	
			// obtain the second fasta file in a "circularised" manner
			circulariseFastaFile(Constants.fasta1, Constants.fasta2, leftShift, rightShift);
		}

		System.out.print("[" + (fragment.getY1()+1) + "," + (fragment.getY1()+1) + "] ");
		System.out.print("[" + (fragment.getY2()+1) + "," + (fragment.getY2()+1) + "] 0\n");


		for ( int index=Chain.size()-1; index>=0; index--)
		{
			pozI = Chain.get(index);
			fragment = this.Fragments.get(pozI.intValue());

			System.out.print("[" + fragment.getX1() + "," + fragment.getY1() + "] ");
			System.out.print("[" + fragment.getX2() + "," + fragment.getY2() + "] " + fragment.getWeight() +"\n");

			if ( index > 0 )
			{
				after = this.Fragments.get(Chain.get(index-1).intValue());
				dif1 = fragment.getX1() - after.getY1()-1;
				dif2 = fragment.getX2() - after.getY2()-1;
				if ( dif1 < 0 ) coverage1 += dif1;
				if ( dif2 < 0 ) coverage2 += dif2;
			}
			coverage1 += fragment.getLength1();
			coverage2 += fragment.getLength2();
		}

		System.out.print("[0,0] [0,0] 0\n");

		if (Constants.OUTPUT_COMPLETE)
		{
			System.out.print("Coverage 1st genome: " + coverage1 + "\n");
			System.out.print("Coverage 2nd genome: " + coverage2 + "\n");
		}

		return output;

	}

	/*
	 * fragment list for recursive step k
	 */ 
	public String Fragments_k_toString (int k)
	{
		List<Integer> fragmentList = this.getFragments_k(k);
		String output = new String("");

		if ( !fragmentList.isEmpty() )
		{
			Fragment fragment;	
			//System.out.println("\nFragments on level " + k + ": " + fragmentList.size());
			for ( Integer pozI: fragmentList )
			{
				fragment = this.Fragments.get(pozI.intValue());
				System.out.print("[" + fragment.getX1() + "," + fragment.getY1() + "] ");
				System.out.print("[" + fragment.getX2() + "," + fragment.getY2() + "] " + fragment.getWeight() +"\n");
			}
		}
		return output;
	}

	/*
	 * once we know the right cutting point, we rearrange fragments
	 */
	public boolean rearrangeFragments(int cutPos1, int cutPos2)
	{	
		Fragment cutFragment1 = Fragments.get(cutPos1);
		Fragment cutFragment2 = Fragments.get(cutPos2);


		/*System.out.println(cutFragment1.getX1() + " " + cutFragment1.getY1() +
				" " + cutFragment1.getX2() + " " + cutFragment1.getY2());
		System.out.println(cutFragment2.getX1() + " " + cutFragment2.getY1() +
				" " + cutFragment2.getX2() + " " + cutFragment2.getY2());
		 */

		//we verify whether the 2 cutting points are separated by several fragments or not
		if ( cutPos1 == 0 || (cutFragment1.getX2() <= cutFragment2.getY2()) )
			leftShift = cutFragment1.getX2()-1;
		else
			leftShift = cutFragment2.getY2()-1;

		/*
		 * compute the shifts on the second genome
		 */
		rightShift = maxY2-leftShift + 1;
		Fragment f;
		for (int i=0; i<Fragments.size(); i++)
		{
			f = Fragments.get(i);
			if ( Fragments.get(i).getX2()>=leftShift)
			{			
				f.setX2(f.getX2()-leftShift);
				f.setY2(f.getY2()-leftShift);
				f.setMovingDirection(-1);
			}
			else
			{
				f.setX2(f.getX2()+rightShift);
				f.setY2(f.getY2()+rightShift);
				f.setMovingDirection(1);
			}
			Fragments.set(i, f);
		}
		/*if ( this.noFragments <= NO_MAX_FRAG )
			selectFragments_k();*/

		return true;
	}

	public abstract void setLevels();
	public abstract void setShifts(int level);
	public abstract int testCircular();
	public abstract boolean findCuttingPoint();

}
