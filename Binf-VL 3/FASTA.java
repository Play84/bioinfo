import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Library with functions to read FASTA format sequence files and perform alignment.
 * Uses ScoringMatrix class
 * v0.1
 * @author kehwan
 *
 */
abstract public class FASTA {
	public int laenge;
	/**
	 * Reads all sequences in the file into memory at once.
	 * Creates an array of strings. Will treat sequences as ascii strings. Warning: assumes upper case with protein!
	 * Index 0 is the first sequence in the file, index 1 is the second, etc. Names are not saved.
	 * @param filename File to read sequences from.
	 * @return String-Array with the sequences.
	 */
	//TODO: Correct fasta parser
	
	static ArrayList<String> readAllSequences(String filename) throws IOException{
		//Read in the first sequence in the fasta file (first command line argument) into buffer and the second into P
		//(Sequences are just ascii strings here)
		//TODO: A correct fasta parser.
		
		BufferedReader source = new BufferedReader(new FileReader(filename));
		ArrayList<String> sequences = new ArrayList<String>();
		String line = source.readLine();
		
		do {
			StringBuilder buffer = new StringBuilder();
			
			//Skip to beginning of the first sequence (or to the end of file, if no sequence)
			while (line != null && (line.length() == 0 || line.charAt(0) != '>'))
				line = source.readLine();
			//now line should contain the comment for a sequence 
				
			//Read until it is over (blank line or >) 
			while((line = source.readLine()) != null && line.length() != 0 && line.charAt(0) != '>')
				buffer.append(line);
			
			//and save the sequence
			if (buffer.length() > 0)
				sequences.add(buffer.toString());
			
		} while (source.ready()); //repeat if there is more in the file
		
		source.close();
		
		return sequences;
	}
	
	/**
	 * Aligns two sequences and returns their distance. O(max(|S1|,|S2|)) space.
	 * @param T ... first sequence to align, as a string (one-letter code amino acid sequence)
	 * @param P ... second sequence
	 * @param S ... substitution matrix
	 * @return
	 */
	static int AlignmentDistance(String T, String P, ScoringMatrix S){
		final int indel_cost = -1;
		
		int m = T.length();
		int n = P.length();
		
		//Calculate scores of the optimal alignments, and the backtracking information.
		int[][] D = new int[n+1][];
				
		//Row 0: alignments of prefixes of the template with ""
		D[0] = new int[m+1];
		for (int j=0;j<=m; j++)
			D[0][j] = j * indel_cost;
			
		//Fill out the matrix row by row. Only save the row currently being filled and the one above. Discard the rest.
		//Fill out the A matrix compeltely though, if alignments are enabled.
		
		for (int i=1;i<=n;i++){	//row i
			D[i] = new int[m+1];
			//Column 0 is for alignments of prefixes of the probe with ""
			D[i][0] = i * indel_cost;
			
			//Fill out row i
			for (int j=1;j<=m;j++){	//cell (i,j)
				
				//Possibility 1 for the optimal alignment of t1...tj with p1...pi is:
				// ~~~ optimal alignment of t1...tj-1~~~   A
				// ~~~          with        p1...pi-1~~~   T
				// This corresponds to T[j] having been mutated into P[i] during the evolution of T into P.
				int score_if_match = D[i-1][j-1] + S.d(P.charAt(i/*string index starts at 0*/-1), T.charAt(j/*string index starts at 0*/-1));
				
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
				
			}
			//done with row i. free up memory used by previous row
			if (i > 0)
				D[i-1] = null;
		}
		//Write out alignment score.
		return D[n][m];
	}

	
	/**
	 * Constructor. Not in use because all methods are class methods
	 * Open FASTA file on Disk.
	 * @param filename
	 *
	FASTAReader(String filename) throws FileNotFoundException{
		source = new BufferedReader(new FileReader(filename));
	}
	*/
}
