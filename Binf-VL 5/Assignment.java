/**
 * Assignment 5, Bioinformatics 2012
 * Uses the excellent 'graphStream' library from http://graphstream-project.org
 * @author kehwan
 *
 */
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.graphstream.algorithm.ConnectedComponents;
import org.graphstream.algorithm.Toolkit;
import org.graphstream.graph.*;
import org.graphstream.graph.implementations.SingleGraph;

public class Assignment {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws ElementNotFoundException 
	 * @throws IdAlreadyInUseException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IdAlreadyInUseException, ElementNotFoundException, IOException, InterruptedException {
		if (args.length <2 ){
			System.out.println("First argument must be name of SIF file to read!");
			System.out.println("Second argument must be a parameter k >= 0 for the determination of the k-cores!");
			System.exit(-1);
		}
		//open the sif file
		SIF s = new SIF(args[0]);
		String[] partner;
		
		//make undirected graph
		Graph g = new SingleGraph("Drosophila");
		g.setStrict(false);
		int edgeID = 1;
		while ( (partner = s.next())  != null){
			g.addNode(partner[0]);
			g.addNode(partner[1]);
			g.addEdge(Integer.toString(edgeID++), partner[0], partner[1], false);
		}
		System.out.println("Done reading file.");
		System.out.println("Proteins:\t" + g.getNodeCount());
		int max_edges = g.getNodeCount() * (g.getNodeCount() - 1)/2; //number of edges in a clique is n*(n-1)/2
		System.out.println("Interactions:\t" + g.getEdgeCount() + " (density " + (100.0*g.getEdgeCount()/max_edges)+"%)");
		
		/*compute cluster coefficient*/

		//there is a function for the cluster coefficient in class Toolkit from the graphStream library, but im afraid we'd get 0pts if we compute it in just 1 line...)
		//we're already avoiding having to implement an adjacency list, etc
		float CC = 0;
		int inner_nodes = 0;	//number of nodes with degree > 1
		for (Node n:g){
			if (n.getDegree() > 1){
				inner_nodes++;
				int links_among_neighbors = 0;
				int max_links = ((n.getDegree())*(n.getDegree() - 1)/2);
				
				//collect neighbors of n
				HashSet<Node> neighbors = new HashSet<Node>();
				for (Edge e : n.getEdgeSet())
					neighbors.add(e.getOpposite(n));
				
				//count links among the neighbors
				for (Node x : (Set<Node>) neighbors.clone()){ //this does exactly what for (Node x : neighbors) ought to do.
					neighbors.remove(x);
					for(Node y: neighbors)
						if (x.hasEdgeBetween(y))
							links_among_neighbors++;
				}
				
				float cc_n = (float) links_among_neighbors / max_links; //clustering coefficient of n
				CC += cc_n; //CC: global clustering coefficient, i.e. average cc_n over all n
				//note: nodes with degree < 2 count as having a clustering coefficient of 0
			}
		}
		CC /= g.getNodeCount();
		System.out.println("Cluster coefficient:\t" + CC);
		System.out.println("Cluster coefficient, computed with library:" + Toolkit.averageClusteringCoefficient(g));
		System.out.println("Cluster coefficient, average over proteins that are not leaves or loners:\t" + (CC * g.getNodeCount()) / inner_nodes);
		System.out.println("Number of proteins that are leaves (1 interaction partner) or loners (0):\t" + (g.getNodeCount() - inner_nodes) + " out of " + g.getNodeCount());
		/*--done with cluster coefficient--*/
		
		//order proteins by degree (nr of interaction partners)
		Node[] nodes = new Node[g.getNodeCount()]; //make array to contain the list of nodes
		g.getNodeSet().toArray(nodes); //put the graph's nodes into it
		Arrays.sort(nodes, new DegreeComparer()); //sort! DegreeComparer defines the total order to sort by, in this case, degree.
		
		System.out.println("The top five proteins by number of interaction partners are:");
		for(int i=0; i < Math.min(5, g.getNodeCount()); i++)
			System.out.println(nodes[i].getId()+"\t"+nodes[i].getDegree());
		
		//Prune graph down using the k-cores method
		int k = Integer.parseInt(args[1]);
		//iterate over all nodes, eliminating nodes that are not part of a k-core, leaving a graph where all connected components are k-cores		
		trimToKcores(g,k,true);
	}
	
	/**
	 * Den gegebenen Graphen auf seine k-cores heruntertrimmen.
	 * Anschließend alles wegschneiden was nicht zum größten k-core gehört.
	 * Ergebnis anzeigen. //FIXME: Zweite und dritte funktion in eigene funktion auslagern...
	 * @param g
	 * @param k
	 * @throws InterruptedException
	 */
	static void trimToKcores(Graph g, int k, boolean visualize_algorithm) throws InterruptedException{
		System.setProperty("gs.ui.renderer", "org.graphstream.ui.j2dviewer.J2DGraphRenderer");
		boolean displaying = false;
		
		/**
		 * Algorithmus zur reduktion eines graphen auf alle seine k-cores, aus der Übung
		 * Schrittweises eliminieren aller knoten die nicht zu irgendeinem k-core gehoeren
		 * Verstehe selbst nicht ganz, wieso das funktioniert
		 */
		int N;
		do {
			//visualize (optional)
			if (visualize_algorithm && !displaying && g.getNodeCount() < 500){
				displaying = true;
				g.display(); //this displays the graph graphically
				Thread.sleep(5000);
			}
			
			//actual algorithm is here
			N = g.getNodeCount();
			Node[] S = new Node[g.getNodeCount()];
			g.getNodeSet().toArray(S); //workaround: for (Node n:g) {...} does not work, because modifying g while iterating over it breaks iteration, and some nodes are not selected for iterating over.
			
			for(Node n: S){				
				int effective_degree = n.getDegree(); //eff-degree: anzahl nachbarn die mindestestens degree k haben
				if (effective_degree < k)
					g.removeNode(n);
				else {
					for (Edge e : n.getEdgeSet())
						if (e.getOpposite(n).getDegree() < k)
							effective_degree--;
					if (effective_degree < k)
						g.removeNode(n);
				}
				if (displaying)
					Thread.sleep(25);
			}
		} while ( g.getNodeCount() < N);
		
		//split graph into ints connected components (=K-cores)
		ConnectedComponents cc_algo = new ConnectedComponents();
		cc_algo.init(g);
		
		//get the largest k-core
		List<Node> kcore;
		kcore = cc_algo.getGiantComponent();
		
		//print out proteins
		System.out.println("The largest "+  k+"-core, with "+kcore.size()+" proteins, consists of:");
			for (Node n: kcore)
				System.out.println(n.getId());
		
		cc_algo.terminate();	
			
		//delete all other nodes from graph	//FIXME: (this is a dirty workaround to mimic the nonexistent operation of 'create subgraph from graph g and node collection kcore'
		Node[] G = new Node[g.getNodeCount()];
		g.getNodeSet().toArray(G);				//FIXME: dirty workaround for doing for (node n : g)
		for (Node n : G)
			if (!kcore.contains(n))
				g.removeNode(n);
		
		//display graph (I.e. largest k-core), if not already displaying
		if (!displaying)
			g.display();

		 	//TODO: configure display to include protein names
			//TODO maybe: make algorithm self-explanatory by improving the video visualization (colors for currently selected neighborhood, ...)
	}
}