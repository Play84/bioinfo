import java.util.Comparator;
import org.graphstream.graph.Node;

/**
 * Comparator, compares GraphStream nodes a and b by degree.
 * Returns negative value if degree(a) < degree(b), zero if degree(a)=degree(b), positive value otherwise.
 * Note: Imposes ordering inconsistent with equals, i.e. the ordering is not a strict ordering. Who'd have thought that.
 * @author kehwan
 *
 */
public class DegreeComparer implements Comparator<Node> {
	public int compare(Node a, Node b) {
		return b.getDegree() - a.getDegree();
	}
}
