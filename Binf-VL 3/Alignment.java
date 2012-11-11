
/**
 * Ad-hoc representation of a sequence alignment
 * @author kehwan
 *
 */
public class Alignment {
	private StringBuilder row1 = new StringBuilder();
	private StringBuilder row2 = new StringBuilder();
	
	/**
	 * Attach the given character pair on the left.
	 * @param top
	 * @param bottom
	 */
	public void append(char top, char bottom){
		row1.insert(0, top);
		row2.insert(0, bottom);
	}
	/**
	 * Returns the number of columns.
	 * @return
	 */
	public int length(){
		return row1.length();
	}
	
	/**
	 * Removes leftmost column.
	 */
	public void shorten(){
		row1.deleteCharAt(0);
		row2.deleteCharAt(0);
	}
	
	/**
	 * Prints out the alignment to stdout.
	 */
	public void print(){
		System.out.println(row1);
		System.out.println(row2);
	}
}
