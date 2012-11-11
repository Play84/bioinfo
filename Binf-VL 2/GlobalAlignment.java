import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
/**
 * For a course ('Grundlagen der Bioinformatik')
 * Needleman-Wunsch algorithm implementation (iirc), for simple global sequence alignment.
 * Returns alignment score and, optionally all optimal alignments.
 * 
 * Score computation is O(|S|+|T|) memory and O(|S|*|T|) time. Alignment computation needs at least quadratic memory in this implementation.  
 * 
 * The score of an alignment is:
 * -m for a single position insertion or deletion
 * +s(a,b) for a single base mutation of base A to base B, or a single amino acid replacement
 * 
 * Harcoded values
 * m=1
 * s(a,b) = 1 if a=b, else -1
 * Sloppy architecture
 * 
 * Note: In the following "probe" refers to the lower sequence in the alignment and the left side of the matrix.
 * The "template" refers to the upper sequence in the alignment and the top side of the matrix.
 * This is all convention; both sequences are probably about equally long (i.e. probe is not shorter)
 * 
 * @author kehwan
 */

enum EditscriptEvent {
	DEL/* - */,		//an arrow pointing at the cell to the left. A gap in the probe. A deletion in the evolution of Template into Probe. (or an insertion in the probe)
	MUT /* \ */,	//an arrow pointing at the upper-left. A match. A mutation.
	INS /* | */		//an arrow pointing up. A gap in the template. An insertion in the template. (or a deletion in the probe)
};

public class GlobalAlignment {
	
	static final int indel_cost = -1;
	static ScoringMatrix S;
	/**
	 * Scoring function (matrix) of a mutation from nucleotide a to nucleotide b.
	 * 
	 * @param a: Nucleic acid base (ascii char, case sensitive)
	 * @param b: Nucleic acid base (as above)
	 * @return
	 */
	static int s(char a, char b){
			//Looking up the score with the scoring Matrix object, S.
			return S.s(a,b);
		}
	
	static String T;	//t1...tm. The first of two sequences to align. Goes on top of alignment. The 'Template'.
	static String P;	//p1...pn. The second of two sequences to align. Goes in bottom of alignment. The 'Probe'
	static int[][] D;	//Score matrix. D(i,j) = Score of highest scoring alignment of t1..tj with p1...pj
	static EnumSet<EditscriptEvent>[][] A;	//A(i,j) is a subset of {MUT, INS, DEL}. 
	
	//Used by the function that prints the optimal alignments using the data from A
	static int i;
	static int j;
	static Alignment alignment;
	
	public static void main(String[] args) throws IOException{
		boolean makeAlignments;
		
		if (args.length < 2){
			System.out.println("Parameter 1 required: Name of a FASTA file with two sequences.");
			System.out.println("Parameter 2 required: Name of a file with a scoring matrix, BLOSUM format.");
			System.out.println("Optional third parameter -onlyScore disables alignment computation (uses less memory)");
			System.exit(0);
		}
		
		if (args.length	== 3 && args[2].equals("-onlyScore"))
			makeAlignments = false;
		else
			makeAlignments = true;
		
		//Read in scoring Matrix from file in argument 2.
		S = new ScoringMatrix(new File(args[1]));
		
		//This block will read in the first sequence in the fasta file (first command line argument) into T and the second into P
		//(Sequences are just ascii strings here)
		//TODO: A correct fasta parser.
		//TODO: Is it really true that the garbage collector does not collect Strings? Then I shouldn't use a String line = like I do. Oh well.
		{
			String line;
			BufferedReader br = new BufferedReader(new FileReader(args[0]));
			StringBuilder t = new StringBuilder();
	
			line = br.readLine();

			//Skip to beginning of the first sequence
			while (line != null && (line.length() == 0 || line.charAt(0) != '>'))
				line = br.readLine();
			
			//Read until it is over (blank line or >) 
			while((line = br.readLine()) != null && line.length() != 0 && line.charAt(0) != '>')
				t.append(line);
			
			//and write sequence into T (template)
			T = t.toString();
			
			//now read next sequence
			t = new StringBuilder();
			
			//Skip to the beginning of sequence
			while (line != null && (line.length() == 0 || line.charAt(0) != '>'))
				line = br.readLine();
			
			//Read until it is over (blank line or >) 
			while((line = br.readLine()) != null && line.length() != 0 && line.charAt(0) != '>')
				t.append(line);
			
			P = t.toString();
			
			br.close();
		}
		//done reading sequences into T and P
		
		int m = T.length();
		int n = P.length();
		
		if (makeAlignments && (/* (n*m) is large */ n * m > Math.pow(10, 9)))	//We would require a really big matrix
			System.out.println("Warning: Memory use will be > 1GB due to long inputs. Use 2nd parameter -onlyScore to not calculate alignments and reduce memory use to O(|Seq1|+|Seq2|)");
		
		//The following section will calculate the scores of the optimal alignments, and optionally the backtracking information.
		D = new int[n+1][]; //distance matrix
		
		if (makeAlignments)
			A = new EnumSet[n+1][m+1]; //If A[i][j] contains, say, {MUT, DEL} then that means that an optimal alignment of T1...Tj with P1...Pi can end with either a mutation or a deletion-in-template, i.e. the corresponding types of alignment column (base-base, or base-gap, iirc)
										//It can be visualized as sets of arrows, or equivalently, as a set of last columns.
		
		//Row 0: alignments of prefixes of the template with ""
		D[0] = new int[m+1];
		for (int j=0;j<=m; j++){
			D[0][j] = j * indel_cost;
			
			if (makeAlignments)
				A[0][j] = EnumSet.of(EditscriptEvent.DEL); //note this corresponds to { <- }, set of left-arrow only
		}
	
		//Fill out the matrix row by row. Only save the row currently being filled and the one above. Discard the rest.
		//Fill out the A matrix compeltely though, if alignments are enabled.
		
		for (int i=1;i<=n;i++){	//row i
			D[i] = new int[m+1];
			//Column 0 is for alignments of prefixes of the probe with ""
			D[i][0] = i * indel_cost;
			if (makeAlignments)
				A[i][0] = EnumSet.of(EditscriptEvent.INS); //set of up-arrow only. 
			
			//Fill out row i
			for (int j=1;j<=m;j++){	//cell (i,j)
				
				//Possibility 1 for the optimal alignment of t1...tj with p1...pi is:
				// ~~~ optimal alignment of t1...tj-1~~~   A			(note: A is T[j] in the general case)
				// ~~~          with        p1...pi-1~~~   C			(note: C is P[i])
				// This corresponds to T[j] having been mutated into P[i] during the evolution of T into P.
				int score_if_match = D[i-1][j-1] + s(P.charAt(i/*string index starts at 0*/-1), T.charAt(j/*string index starts at 0*/-1));
				
				//2: Possibility 2 is:
				// ~~~ optimal alignment of t1...tj-1~~~   A
				// ~~~          with        p1...pi  ~~~   -
				// This corresponds to T[j] having been deleted.
				int score_if_del = D[i][j-1] + indel_cost;
				
				//3: Possibility 3 is:
				// ~~~ optimal alignment of t1...tj    ~~~   -
				// ~~~          with        p1...pi-1  ~~~   C
				// This corresponds to P[j] having been inserted.
				int score_if_ins = D[i-1][j] + indel_cost;
				
				//Calculate optimal score
				D[i][j] = Math.max(Math.max(score_if_match, score_if_del), score_if_ins);
				
				if (makeAlignments){
					//Record what the optimal alignment of t1...tj with p1...pi was.
					A[i][j] = EnumSet.noneOf(EditscriptEvent.class);
					if (score_if_match == D[i][j])
						A[i][j].add(EditscriptEvent.MUT);	//Could it have been possibility 1?
					if (score_if_ins == D[i][j])
						A[i][j].add(EditscriptEvent.INS);	//or possibility 2?
					if (score_if_del == D[i][j])
						A[i][j].add(EditscriptEvent.DEL);	//or possbility 3?
				}
			}
			//done computing scores in row i. free up memory used by previous row
			if (i > 0)
				D[i-1] = null;
		}
		//Write out alignment score.
		System.out.println("The best possible alignment score is: " + D[n][m]);
		//Output optimal alignments
		if (makeAlignments){
			System.out.println("The alignments with this score are:\n");
			i = n;
			j = m;
			alignment = new Alignment();			
		    recursivelyPrintAlignments();
		}
		else
			System.out.println("You have disabled alignment output.");
	}

	/**
	 * Print all optimal alignments, corresponding to all paths from (i,j) to (0,0) 
	 * in the graph determined by the 'arrows' in matrix A 
	 * 
	 * This sloppy function uses many global variables:
	 * A, T, P, i, j, alignment		(at the time of this writing)
	 * 
	 * This allows for a recursion where the stack frames save nothing. 
	 * (i.e. no memory is allocated in stackframes for locals such as alignment, or even i, j)
	 */
	static void recursivelyPrintAlignments(){
		if (i==0 && j == 0){
			alignment.print();
			System.out.println();
			return;
		}
		else {	
			if (A[i][j].contains(EditscriptEvent.DEL /* <- */)){
				alignment.append(T.charAt(j-1), '_');
				j--;
				recursivelyPrintAlignments();
				//now backtracking...
				j++;	
				alignment.shorten();
			}
			
			if (A[i][j].contains(EditscriptEvent.MUT /* \ */)){
				alignment.append(T.charAt(j-1), P.charAt(i-1));
				i--; j--;
				recursivelyPrintAlignments();
				//now backtracking...
				i++; j++; 
				alignment.shorten();
			}
			
			if (A[i][j].contains(EditscriptEvent.INS /* | */)){
				alignment.append('_', P.charAt(i-1));
				i--;
				recursivelyPrintAlignments();
				//now backtracking...
				i++;	
				alignment.shorten();
			}
			
			return;
		}
	}
}	
	
	/*
	 * The overhead of managing all those function calls is considerable. Edit: No, it's absurd!
	 * Will this run out of memory when both sequences are large, but have a small number of optimal alignments?
	 * Yes, currently for each column in the alignment there is a function call put on the stack.
	 * That is not nice, even if we avoid holding a deep copy of the current alignment as a local variable. (that would case O(oh no(|T|, |S|)) memory use)
	 */
	 //TODO: Let's make it that for every *branch* there is a function call on the stack.
	
	//TODO: Find a library for finding paths in implicit directed acyclic graphs
	//I want to write the code below, which is equivalent to the code in my program.
	
	/*This matrix A specifies a directed acyclic graph. Now:
	 * For each path in the graph from node (n,m) to (0,0):
	 * 		 alignment = new Alignment();
	 * 		 for each node (i,j) on the path (in path order), where e is the type of the edge through which we entered (i,j):
	 *  		 if edgetype == diagonal
	 *  			alignment.add(T[j+1], P[i+1])
	 *  	     else if edgetype == left
	 *              alignment.add(T[j+1], '_')
	 *           else if edgetype == up
	 *              alignment.add('_', P[j+1])
	 *       
	 *       alignment.print()
	 */