import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Grundlagen der Bioinformatik 2011/12 (Leser, Thomas) HU berlin
 * Aufgabe 3 aus der Uebung
 * 
 * Kommandozeilenprogramm: 
 * Syntax: Treemaker sequences.fasta submatrix.blosum
 * 
 * Bestimmt mit dem Verfahren aus der Vorlesung (Hierarchical clustering) einen Stammbaum zu den Proteinsequenzen im Fasta-File sequences.fasta 
 * @author kehwan
 *
 */
public class TreeMaker {

	/**
	 * 
	 * @param args
	 * @throws IOException 
	 * 
	 */
	public static void main(String[] args) throws IOException {
		//0: Substi-Matrix einlesen
		ScoringMatrix S = new ScoringMatrix(args[1]);
		//1: Alle sequenzen einlesen
		ArrayList<String> sequences = FASTA.readAllSequences(args[0]);
		
		//2: Tabelle der paarweisen Distanzen ausrechen (D(i,j) = Alignment-Abstand von Spezies i zu Spezies j)	
		
		//2.1 hilfsweise legen wir eine liste der speziesnamen an zu denen derzeigt abstände gespeichert sind.
		HashSet<String> spezies = new HashSet<String>();
		for(int i=0; i<sequences.size(); i++)
			spezies.add(Integer.toString(i));
		
		//2.2 Alle Distanzen zwischen sequenzpaaren berechnen
		HashMap<String, Integer> D = new HashMap<String, Integer>();
		for (String spezies1: spezies)
			for (String spezies2: spezies){
				if (! spezies1.equals(spezies2)){
					//speichere den abstand zwischen spezies1 und spezies 2 ab
					D.put(spezies1 + "," + spezies2, FASTA.AlignmentDistance(/*spezies1-sequenz*/sequences.get(Integer.parseInt(spezies1)), /*spezies2-sequenz*/sequences.get(Integer.parseInt(spezies2)), S));
					D.put(spezies2 + "," + spezies1, D.get(spezies1 + "," + spezies2));
				}
			}
		
		/*3: solange die tabelle daten zu mehr als einer spezies entählt:
		 * 		suche das am nächsten verwandte sequenzpaar (a,b)
		 * 		drucke (a,b) aus - zeigt an dass (a,b) zusammengeführt wurde (gemeinsamer vorfahre)
		 * 		loesche die zeilen und spalten zu a und b aus der tabelle und füge einen eintrag für einen gemeinsamen vorfahren ab ein
		 * 		der abstand D(x,ab) der anderen spezies x zum vorfahren ab ist Mittelwert(D(x,a),D(x,b)). Das wird in die tabelle eingetragen
		 */
		while (D.size() > 0){
			//3.1 suche maximum in der tabelle - d.h. die beiden Spezies (a,b) die am nächsten verwandt sind.
			int bestScore = Integer.MIN_VALUE;
			String bestPair = null;
			for(String speziespaar: D.keySet()){
				if (D.get(speziespaar) > bestScore){
					bestPair = speziespaar;
					bestScore = D.get(bestPair);
				}
			}
			// das Paar der 2 maximal verwandten Spezies (zB: "14,5") steht nun in der variable 'bestPair'. 
			
			//3.2 gib (a,b) aus
			System.out.print("(" + bestPair + ") ");
			//System.out.println(": " + bestScore);
			
			//3.3.  Wir updaten die Tabelle so dass *statt* den spezies a und b nun die spezies ab drin steht, mit entsprechendem abstand zu den anderen spezies.
			//3.3.1 Wir nennen den Vorfahren von a und b "ab"
			String a = bestPair.split(",")[0];
			String b = bestPair.split(",")[1];
			String Vorfahre = a + b;
			//3.3.2 Wir fügen zu D die Abstände von allen Spezies x (außer a und b) zur neuen Spezies ab dazu
			// und wir entfernen aus D die Abstände von x zu a bzw. zu b.
			for (String x:spezies){
				if (! x.equals(a) && ! x.equals(b)){
					//abstand von x zu ab hinzufügen
					D.put(x + "," + Vorfahre, /* (Abstand(x,a)+Abstand(x,b))/2 */ (D.get(x + "," + a) + D.get(x + "," + b)) / 2 );
					D.put(Vorfahre + "," + x, D.get(x + "," + Vorfahre));
					
					//abstand von x zu a loeschen
					D.remove(x + "," + a); D.remove(a + "," + x);
					//abstand von x zu b loeschen
					D.remove(x + "," + b); D.remove(b + "," + x);
				}
			}
			
			//Der abstand von a zu b muss auch noch geloescht werden.
			D.remove(a + "," + b); D.remove(b + "," + a);
			
			//3.3.3 wir loeschen a und b aus der liste der speziesnamen und geben den Vorfahren, ab, hinzu.
			spezies.remove(a);
			spezies.remove(b);
			spezies.add(Vorfahre);
		}
		
		//Stammbaum sollte nun ausgegeben worden sein.
	}
}
