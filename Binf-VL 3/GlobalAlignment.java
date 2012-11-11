import java.io.IOException;
import java.util.*;
/**
 * For a course ('Grundlagen der Bioinformatik')
 * Needleman-Wunsch algorithm implementation (iirc), for simple global sequence alignment.
 * Returns alignment score and, optionally all optimal alignments.
 * 
 * Score computation is O(|S|+|T|) memory and O(|S|*|T|) time. Alignment computation needs at least quadratic memory.  
 * 
 * The score of an alignment is:
 * -1 for a single position insertion or deletion
 * +S(a,b) for a single mutation of A to B - Reads a blosum substitution matrix into S
 * 
 * Sloppy architecture
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
	 * Scoring function (scoring matrix) of a mutation from nucleotide a to nucleotide b.
	 * 
	 * @param a: Nucleic acid base (ascii char, case sensitive)
	 * @param b: Nucleic acid base (as above)
	 * @return
	 */
	static int s(char a, char b){
			//Looking up the score has been outsourced to the substitution Matrix object.
			return S.d(a,b);
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
			System.out.println("Parameter 2 required: Name of a file with a substiution matrix, BLOSUM format.");
			System.out.println("Optional third parameter -onlyScore disables alignment computation (uses less memory)");
			System.exit(0);
		}
		
		if (args.length	== 3 && args[2].equals("-onlyScore"))
			makeAlignments = false;
		else
			makeAlignments = true;
		
		//Read in Substiution Matrix from file in argument 2.
		S = new ScoringMatrix(args[1]);
		
		//Read in the first sequence in the fasta file (first command line argument) into T and the second into P
		//(Sequences are just ascii strings here)
	
		ArrayList<String> sequences = FASTA.readAllSequences(args[0]);
		T = sequences.get(0);
		P = sequences.get(1);
		int m = T.length();
		int n = P.length();
		sequences = null;
		
		if (makeAlignments && (n * m > Math.pow(10, 9)))	//We would require more than 10^9 cells in the matrix
			System.out.println("Warning: Memory use will be > 1GB due to long inputs. Use 2nd parameter -onlyScore to not calculate alignments and reduce memory use to O(|Seq1|+|Seq2|)");
		
		//Calculate scores of the optimal alignments, and the backtracking information.
		D = new int[n+1][];
		
		if (makeAlignments)
			A = new EnumSet[n+1][m+1];
		
		//Row 0: alignments of prefixes of the template with ""
		D[0] = new int[m+1];
		for (int j=0;j<=m; j++){
			D[0][j] = j * indel_cost;
			
			if (makeAlignments)
				A[0][j] = EnumSet.of(EditscriptEvent.DEL);
		}
	
		//Fill out the matrix row by row. Only save the row currently being filled and the one above. Discard the rest.
		//Fill out the A matrix compeltely though, if alignments are enabled.
		
		for (int i=1;i<=n;i++){	//row i
			D[i] = new int[m+1];
			//Column 0 is for alignments of prefixes of the probe with ""
			D[i][0] = i * indel_cost;
			if (makeAlignments)
				A[i][0] = EnumSet.of(EditscriptEvent.INS);
			
			//Fill out row i
			for (int j=1;j<=m;j++){	//cell (i,j)
				
				//Possibility 1 for the optimal alignment of t1...tj with p1...pi is:
				// ~~~ optimal alignment of t1...tj-1~~~   A
				// ~~~          with        p1...pi-1~~~   T
				// This corresponds to T[j] having been mutated into P[i] during the evolution of T into P.
				int score_if_match = D[i-1][j-1] + s(P.charAt(i/*string index starts at 0*/-1), T.charAt(j/*string index starts at 0*/-1));
				
				//2: Possibility 2 is:
				// ~~~ optimal alignment of t1...tj-1~~~   A
				// ~~~          with        p1...pi  ~~~   -
				// This corresponds to T[j] having been deleted.
				int score_if_del = D[i][j-1] + indel_cost;
				
				//3: Possibility 3 is:
				// ~~~ optimal alignment of t1...tj    ~~~   -
				// ~~~          with        p1...pi-1  ~~~   T
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
			//done with row i. free up memory used by previous row
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
	 */
	static void recursivelyPrintAlignments(){
		if (i==0 && j == 0){
			alignment.print();
			System.out.println();
			return;
		}
		else {	
			if (A[i][j].contains(EditscriptEvent.DEL /* -> */)){
				alignment.append(T.charAt(j-1), '_');
				j--;
				recursivelyPrintAlignments();
				j++;	//now backtracking...
				alignment.shorten();
			}
			
			if (A[i][j].contains(EditscriptEvent.MUT /* \ */)){
				alignment.append(T.charAt(j-1), P.charAt(i-1));
				i--; j--;
				recursivelyPrintAlignments();
				i++; j++; //now backtracking...
				alignment.shorten();
			}
			
			if (A[i][j].contains(EditscriptEvent.INS /* | */)){
				alignment.append('_', P.charAt(i-1));
				i--;
				recursivelyPrintAlignments();
				i++;	//now backtracking...
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
	 * That is not nice.
	 */
	 //TODO: Let's make it that for every *branch* there is a function call on the stack.