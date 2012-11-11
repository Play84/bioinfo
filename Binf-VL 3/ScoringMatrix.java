import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Substitution matrix for use with biological sequence alignment.
 * An instance is a particular Matrix.
 * 
 * Risky cheapskate implementation assumes all seq
 * @author kehwan
 *
 */
public class ScoringMatrix {
	private short matrix [][] = new short[256][256];	//the actual matrix is in here. Characters are encoded as ASCII letters.
	private boolean matrixIsDefinedAt [][] = new boolean[256][256];
	
	/**
	 * Read substitution Matrix from file.
	 * Expected file format: blosum
	 * 
	 * @param fastafile
	 * @throws IOException 
	 */
	ScoringMatrix(String blosumfile_name) throws IOException{
		//Open file
		BufferedReader source = new BufferedReader(new FileReader(blosumfile_name));
		String line;
		//Skip comment lines at the beginning of the file 
		do {
			line = source.readLine();
		} while (line != null && line.charAt(0) == '#');
		//Now we have the first non-comment line of the file in 'line'. (or null if there was only comments)
		String[] spaltennamen = line.split("[ \\t]+");
		//Now we've got a string array such as ["", "A","R","N","D"..."*"] (see blosum.txt from aufgabe 3)
				
		while((line = source.readLine()) != null && ! line.equals("")){
			//Wir haben jetzt eine matrixzeile in line.
			String[] lineAsArray = line.split("[ \\t]+");
			
			char zeile = lineAsArray[0].charAt(0);
			
			for(int spaltennr =1; spaltennr<lineAsArray.length; spaltennr++){
				char spalte = spaltennamen[spaltennr].charAt(0);
				matrix[zeile][spalte] = Short.parseShort(lineAsArray[spaltennr]);
				matrixIsDefinedAt[zeile][spalte] = true;
			}
		}
	}
 
	/**
	 * Scoring function (scoring matrix) of a mutation from amino acid a to amino acid b.
	 * 
	 * @param a: Amino acid a (ascii char, case sensitive, should be upper case)
	 * @param b: Amino acid b (as above)
	 * @return
	 */
	public short d(char a, char b){
		if (matrixIsDefinedAt[a][b]) 
			return matrix[a][b];
		else {
			System.out.println("Error: Scoring matrix did not contain a score for the pair of sequence alphabet elements (" + a + "," + b + ")");
			System.out.println("Wrong capitalization in matrix or input sequence? \n aborting.");
			System.exit(-1);
			return matrix[a][b]; //just so that eclipse does not complain; unreachable
		}
	}
}
