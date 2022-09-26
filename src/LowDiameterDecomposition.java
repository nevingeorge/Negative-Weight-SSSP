import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.Set;

public class LowDiameterDecomposition {

	public static void main(String[] args) {

	}
	
	// returns all the vertices in g within a distance of r from source vertex s using Dijkstra's
	public static ArrayList<Integer> volume(Graph g, int s, double r) {
		ArrayList<Integer> output = new ArrayList<Integer>();
		
		Set<Integer> settled = new HashSet<Integer>();
	    PriorityQueue<Node> pq = new PriorityQueue<Node>(g.n, new Node());
		double[] dist = new double[g.n];

        for (int i = 0; i < g.n; i++) {
            dist[i] = Integer.MAX_VALUE;
        }
 
        pq.add(new Node(s, 0));
        dist[s] = 0;
 
        while (settled.size() != g.n) {
            if (pq.isEmpty()) {
                return output;
            }

            int u = pq.remove().node;

            if (settled.contains(u) || dist[u] > r) {
                continue;
            }
            
            output.add(u);
            settled.add(u);
            ArrayList<Integer> neighbors = g.getEdges(u);

            for (int i = 0; i < neighbors.size(); i++) {
                int v = neighbors.get(i);

                if (!settled.contains(v)) {
                    double newDistance = dist[u] + g.getWeight(u, v);

                    if (newDistance < dist[v]) {
                        dist[v] = newDistance;
                    }
                    
                    if (dist[v] <= r) {
                    	pq.add(new Node(v, dist[v]));
                    }
                }
            }
        }

        return output;
	}
	
	// run Dijkstra's using priority queues in O(V + ElogV)
	// guaranteed that g has non-negative edges
	public static double[] Dijkstra(Graph g, int s) {		
		Set<Integer> settled = new HashSet<Integer>();
	    PriorityQueue<Node> pq = new PriorityQueue<Node>(g.n, new Node());
		double[] dist = new double[g.n];

        for (int i = 0; i < g.n; i++) {
            dist[i] = Integer.MAX_VALUE;
        }
 
        pq.add(new Node(s, 0));
        dist[s] = 0;
 
        while (settled.size() != g.n) {
            if (pq.isEmpty()) {
                return dist;
            }

            int u = pq.remove().node;

            if (settled.contains(u)) {
                continue;
            }

            settled.add(u);
            ArrayList<Integer> neighbors = g.getEdges(u);

            for (int i = 0; i < neighbors.size(); i++) {
                int v = neighbors.get(i);

                if (!settled.contains(v)) {
                    double newDistance = dist[u] + g.getWeight(u, v);

                    if (newDistance < dist[v]) {
                        dist[v] = newDistance;
                    }

                    pq.add(new Node(v, dist[v]));
                }
            }
        }
        
        return dist;
	}
}

class Node implements Comparator<Node> {
    public int node;
    public double cost;
 
    public Node() {}
 
    public Node(int node, double cost) {
        this.node = node;
        this.cost = cost;
    }

    @Override public int compare(Node node1, Node node2) {
        if (node1.cost < node2.cost) {
            return -1;
        }
        if (node1.cost > node2.cost) {
            return 1;
        }
 
        return 0;
    }
}