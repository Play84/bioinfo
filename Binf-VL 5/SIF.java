import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * Wraps SIF files
 * @author kehwan
 *
 */
public class SIF {
	private BufferedReader r;
	
	/**
	 * Opens a sif file for reading
	 * @param filename
	 * @throws FileNotFoundException
	 */
	SIF(String filename) throws FileNotFoundException{
		r = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
	}
	
	/**
	 * Closes file. Runs automatically when object is garbage collected.
	 */
	protected void finalize(){
		if (r!=null)
			try {
				r.close();
			} catch (IOException e) {System.out.println("[SIF reader]Warning: finalizer failed to close a file. This may cause the program to use more memory.");}
	}
	/**
	 * Get the next pair of IDs of interacting proteins, as an array of two strings.
	 * Returns null when there are no more interactions in the file.
	 * @return Array with two protein names, or null if no more entries in file.
	 * @throws IOException
	 */
	String[] next() throws IOException{
		String line;
		assert(r != null);
		//Read the next line, expecting format PROTEINNAME pp PROTEINNAME. Skip over any lines not having this format, issuing warnings.
		while ((line = r.readLine()) != null && ! line.matches("^\\s*(\\w|.)+\\s+pp\\s+(\\w|.)+\\s*$"))
			System.out.println("[SIF reader]Warning: the following line did not have the expected format and was skipped:\t"+line);
		
		if (line == null)
			return null;
		else {
			String[] result = line.split("pp");
			result[0] = result[0].trim();
			result[1] = result[1].trim();
			return result;
		}
	}
}
