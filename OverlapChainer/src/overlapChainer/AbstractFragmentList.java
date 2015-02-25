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

package overlapChainer;

import java.io.BufferedReader;

import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

/**
 * @author Raluca Uricaru
 *
 */

public abstract class AbstractFragmentList {
	
	protected File fragmentFile;
	protected int noFragments;
	protected List<Fragment> Fragments; 
	
	// biggest values for Y1 and Y2 between all fragments
	protected int maxY1, maxY2;
	
	protected List<Integer> Chain;
	
	public List<Integer> getChain( ) {
		return Chain;
	}
	public File getFragmentFile() {
		return fragmentFile;
	}
	public void setFragmentFile(File fragmentFile) {
		this.fragmentFile = fragmentFile;
	}
	public List<Fragment> getFragments() {
		return Fragments;
	}
	public void setFragments(List<Fragment> fragments) {
		Fragments = fragments;
	}
	public int getNoFragments() {
		return noFragments;
	}
	public void setNoFragments(int noFragments) {
		this.noFragments = noFragments;
	}
	
	public AbstractFragmentList () {
		
	}
	
	public AbstractFragmentList (File fragmentFile) {
		
		Fragments = new ArrayList<Fragment>();
		this.noFragments = 0;
		this.fragmentFile = fragmentFile;
		
		Chain = new ArrayList<Integer>(); 
		
		maxY1 = maxY2 = 0;
		
		/*
		 * read fragments from .mat file
		 */
		readFragmentFile();
	}
	
	public AbstractFragmentList (File fragmentFile, boolean ok) {
		
		Fragments = new ArrayList<Fragment>();
		this.noFragments = 0;
		this.fragmentFile = fragmentFile;
		
		Chain = new ArrayList<Integer>(); 
		
		maxY1 = maxY2 = 0;
		
		/*
		 * read fragments from .mat file
		 */
		this.readFragmentFile(ok);
	}
	
	/*
	 * read fragments from .mat file
	 * input file for Chainer and Hierarchical
	 */
	public void readFragmentFile() {
		int noLine = 0;
		int noTok = 0;
		String line = new String ();
		
		try{
			BufferedReader br = new BufferedReader(new FileReader(this.fragmentFile));	
			String tokens[];
			
			/*
			 * read the first line >CHA 2
			 */
			line = br.readLine();
		
			/*
			 * start reading the fragments
			 */
			Fragment frag = new Fragment();
			Pattern p = Pattern.compile(",|\\s|\\[|\\]");
			while ( (line = br.readLine()) != null ){
					noLine++;
					if ( noLine%2 != 0 )
					{
						frag = new Fragment();
						frag.setWeight((new Double (Double.parseDouble(line.substring(1)))).intValue());
					}
					else
					{
						tokens = line.split(p.pattern());
						for (noTok=1; noTok<tokens.length; noTok++){
								if ( noTok == 1 )
									frag.setX1(new Integer(tokens[noTok]));
								else
								if ( noTok == 2 )
									frag.setY1(new Integer(tokens[noTok]));
								else
								if ( noTok == 5 )
									frag.setX2(new Integer(tokens[noTok]));
								else
								if ( noTok == 6 )
									frag.setY2(new Integer(tokens[noTok]));
						}
						if ( frag.getY1() > maxY1 )
							maxY1 = frag.getY1();
						if ( frag.getY2() > maxY2 )
							maxY2 = frag.getY2();
						frag.setIdp(0);
						frag.computeLengths();
						this.noFragments++;
						this.Fragments.add(frag);
					}
			}
			br.close();
		}
		catch (Exception e){
			System.out.println("Reading error: "+e.toString());
			System.out.println(line);
		}					
	}
	
	
	protected abstract int No_Max_Frag_Read();
	
	
	/*
	 * read fragments from .mat file
	 * input file for Chainer and Hierarchical
	 * we call this function if we encountered java heap error
	 * we read max NO_MAX_FRAG_READ frag from .mat file
	 */
	public void readFragmentFile( boolean ok ) {
		int noLine = 0;
		int noTok = 0;
		
		int noMaxFragRead = this.No_Max_Frag_Read();
			
		//System.out.println(noMaxFragRead);
		
		try{
			BufferedReader br = new BufferedReader(new FileReader(this.fragmentFile));	
			String line = new String ();
			String tokens[];
			
			/*
			 * read the first line >CHA 2
			 */
			line = br.readLine();
			
			/*
			 * start reading the fragments
			 */
			Fragment frag = new Fragment();
			Pattern p = Pattern.compile(",|\\s|\\[|\\]");
			while ( (line = br.readLine()) != null && Fragments.size()<noMaxFragRead ){
					noLine++;
					if ( noLine%2 != 0 )
					{
						frag = new Fragment();
						frag.setWeight((new Double (Double.parseDouble(line.substring(1)))).intValue());
					}
					else
					{
						tokens = line.split(p.pattern());
						for (noTok=1; noTok<tokens.length; noTok++){
								if ( noTok == 1 )
									frag.setX1(new Integer(tokens[noTok]));
								else
								if ( noTok == 2 )
									frag.setY1(new Integer(tokens[noTok]));
								else
								if ( noTok == 5 )
									frag.setX2(new Integer(tokens[noTok]));
								else
								if ( noTok == 6 )
									frag.setY2(new Integer(tokens[noTok]));
						}
						if ( frag.getY1() > maxY1 )
							maxY1 = frag.getY1();
						if ( frag.getY2() > maxY2 )
							maxY2 = frag.getY2();
						frag.setIdp(0);
						frag.computeLengths();
						this.noFragments++;
						this.Fragments.add(frag);
					}
			}
			br.close();
		}
		catch (Exception e){
			System.out.println("Reading error: "+e.toString());
		}					
	}
	
	
	public String Fragments_toString()
	{
		String output = new String("");
		//System.out.println("\nTotal number of fragments: " + Fragments.size());
		
		for ( Fragment fragment: this.Fragments )
		{
			//System.out.print(fragment.getLevel() + "\t");
			System.out.print("[" + fragment.getX1() + "," + fragment.getY1() + "] ");
			System.out.print("[" + fragment.getX2() + "," + fragment.getY2() + "] " + fragment.getWeight() +"\n");
			
		}
		return output;
	}
	
	public abstract String Chain_toString();

	
}
