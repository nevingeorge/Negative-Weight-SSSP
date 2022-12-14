import java.util.*;
import org.apache.commons.rng.sampling.distribution.GeometricSampler;
import org.apache.commons.rng.simple.RandomSource;

// input graph G is guaranteed to have non-negative integer edge weights
public class LowDiameterDecomposition {

	public static void main(String[] args) throws Exception {
	}
	
	// calculates the SCCs of g and runs LDD only on SCCs that have large weak diameter
	public static ArrayList<int[]> preLDD(Graph g, int d) throws Exception {
		Random random = new Random();
		if (NegativeWeightSSSP.YES_RANDOM_SEED) {
			random.setSeed(NegativeWeightSSSP.RANDOM_SEED);
		}
		
		// only calculate SCCs CALCULATE_SCC_PROB percent of the time
		if (random.nextDouble() < NegativeWeightSSSP.CALCULATE_SCC_PROB) {
			return LDD(g, d);
		}
		
		ArrayList<ArrayList<Integer>> SCCs = g.SCC();
		ArrayList<int[]> E_sep = new ArrayList<int[]>();
		
		if (SCCs.size() == 1) {
			return LDD(g, d);
		}
		
		for (ArrayList<Integer> SCC : SCCs) {
			if (SCC.size() > 1) {
				Graph SCCSubgraph = new Graph(g.v_max, false);
				SCCSubgraph.addVertices(SCC);
				
				HashSet<Integer> SCCVerts = new HashSet<>(SCC);
				
				for (int v : SCC) {
					ArrayList<Integer> outVertices = new ArrayList<Integer>();
					ArrayList<Integer> weights = new ArrayList<Integer>();
					
					for (int i = 0 ; i < g.adjacencyList[v].length; i++) {
						if (SCCVerts.contains(g.adjacencyList[v][i])) {
							outVertices.add(g.adjacencyList[v][i]);
							weights.add(g.weights[v][i]);
						}
					}
					
					SCCSubgraph.addEdges(v, listToArr(outVertices), listToArr(weights));
				}
				SCCSubgraph.initNullAdjListElts();
				
				int src = SCC.get((int) (random.nextDouble() * SCC.size()));				
				if (hasLargeDiameter(SCCSubgraph, src, d)) {
					E_sep.addAll(LDD(SCCSubgraph, d));
				}
			}
		}
		
		return E_sep;
	}
	
	public static boolean hasLargeDiameter(Graph g, int s, int diameter) {		
		boolean[] settled = new boolean[g.v_max];
		int numSettled = 0;
	    PriorityQueue<Node> pq = new PriorityQueue<Node>(g.v_max, new Node());
		int[] dist = new int[g.v_max];
		for (int i = 0; i < g.v_max; i++) {
            dist[i] = Integer.MAX_VALUE;
        }
        pq.add(new Node(s, 0));
        dist[s] = 0;
 
        while (numSettled != g.n) {
            if (pq.isEmpty()) {
                return false;
            }

            int u = pq.remove().node;

            if (settled[u]) {
                continue;
            }
            
            if (dist[u] > diameter) {
            	return true;
            }

            settled[u] = true;
            numSettled++;
            
            for (int i = 0; i < g.adjacencyList[u].length; i++) {
            	int v = g.adjacencyList[u][i];
            	
                if (!settled[v]) {
                    int newDistance = dist[u] + g.weights[u][i];

                    if (newDistance < dist[v]) {
                        dist[v] = newDistance;
                    }
                    
                    pq.add(new Node(v, dist[v]));
                }
            }
        }
        
        return false;
	}
	
	// OUTPUT: A set of edges e_sep with the following guarantees:
	// ??? each SCC of G\e_sep has weak diameter at most D; that is, if u,v are in the same SCC,
	// then dist_G(u, v) ??? D and dist_G(v, u) ??? D.
	// ??? For every e ??? E, Pr[e ??? e_sep] = O(w(e)??(logn)^2/D +n^???10). These probabilities are not
	// guaranteed to be independent.
	// Each int[] in the output ArrayList has size two and represents an edge (int[0], int[1])
	public static ArrayList<int[]> LDD(Graph g, int d) throws Exception {
		if (NegativeWeightSSSP.PRINT_LDD_SIZE) {
			System.out.println("Size: " + g.n);
		}

		if (g.n <= Math.max(1, NegativeWeightSSSP.LDD_BASE_CASE)) {
			return new ArrayList<int[]>();
		}
		
		Graph g_rev = createGRev(g);
		
		// pick a good, well-connected starting vertex s
		int s = g.vertices.get(0);
		boolean foundGoodS = false;
		for (int v : g.vertices) {
			if (g.adjacencyList[v].length >= 1 || g_rev.adjacencyList[v].length >= 1) {
				s = v;
				foundGoodS = true;
				break;
			}
		}
		
		if (!foundGoodS) {
			return new ArrayList<int[]>();
		}

		int[] condAndi_max = CoreOrLayerRange(g, g_rev, s, d);
		if (condAndi_max[0] == 1) {
			return RandomTrim(g, g_rev, s, d);
		}
		
		int r = (int) Math.ceil(d / (3.0 * Math.log(g.n)));
		int i_min = condAndi_max[1] - r;
		int i_tilda = GeometricSampler.of(RandomSource.MT.create(), calculateGeoProb(g.n, r)).sample();
		int i_rnd = i_min + Math.min(i_tilda, r);
		
		if (condAndi_max[0] == 2) {
			ArrayList<Integer> ball = volume(g, s, i_rnd);
			Graph subGraph = getSubgraph(g, ball, false);
			Graph minusSubGraph = getSubgraph(g, ball, true);
			
			return edgeUnion(layer(g, ball), preLDD(subGraph, d), preLDD(minusSubGraph, d));
		}
		
		if (condAndi_max[0] == 3) {
			ArrayList<Integer> ball = volume(g_rev, s, i_rnd);
			Graph subGraph = getSubgraph(g_rev, ball, false);	
			Graph minusSubGraph = getSubgraph(g_rev, ball, true);
			
			return revEdges(edgeUnion(layer(g_rev, ball), preLDD(subGraph, d), preLDD(minusSubGraph, d)));
		}
		
		throw new Exception("LowDiamDecomposition failed.");
	}
	
	
	public static ArrayList<int[]> revEdges(ArrayList<int[]> edges) {
		ArrayList<int[]> revEdgeSet = new ArrayList<int[]>();
		for (int[] edge : edges) {
			int[] revEdge = new int[2];
			revEdge[0] = edge[1];
			revEdge[1] = edge[0];
			revEdgeSet.add(revEdge);
		}
		
		return revEdgeSet;
	}
	
	
	public static double calculateGeoProb(int n, int r) {
		double prob = Math.pow(Math.log(n), 2) / r;
		return Math.min(1, prob);
	}
	
	
	public static ArrayList<int[]> RandomTrim(Graph g, Graph g_rev, int s, int d) throws Exception {	
		ArrayList<int[]> e_sep = new ArrayList<int[]>();
		
		if (Integer.MAX_VALUE / 4 <= d) {
			return e_sep;
		}

		int[] dist = Dijkstra(g, s);
		int[] dist_rev = Dijkstra(g_rev, s);
		
		int numDistOver4D = 0;
		int over4DVert = -1;
		
		ArrayList<Integer> v_far = new ArrayList<Integer>();
		for (int v : g.vertices) {
			if (Math.max(dist[v], dist_rev[v]) > 2 * d) {
				v_far.add(v);
				
				if (Math.max(dist[v], dist_rev[v]) > 4 * d) {
					numDistOver4D++;
					over4DVert = v;
				}
			}
		}
		
		// base case
		if (numDistOver4D <= 1) {
			// already essentially a SCC of low diameter other than one vertex
			// just add the edges coming out of over4DVert and nothing else
			if (numDistOver4D == 1) {
				for (int v : g.adjacencyList[over4DVert]) {
					int[] edge = {over4DVert, v};
					e_sep.add(edge);
				}
			}			
			
			return e_sep;
		}
		
		ArrayList<Integer> m = new ArrayList<Integer>(); // marked vertices
		int i_max = d;
		int r = (int) Math.ceil(d / (3.0 * Math.log(g.n)));
		int i_min = i_max - r;
		
		int v = diffVertex(v_far, m, g.v_max);
		while (v != -1) {
			int i_rnd = i_min + Math.min(GeometricSampler.of(RandomSource.MT.create(), calculateGeoProb(g.n, r)).sample(), r);
			
			if (dist_rev[v] > 2 * d) {
				Graph gVMinusM = getSubgraph(g, m, true);
				ArrayList<Integer> ball = volume(gVMinusM, v, i_rnd);
				Graph GVMinusMSubGraph = getSubgraph(gVMinusM, ball, false);
				e_sep = edgeUnion(e_sep, layer(gVMinusM, ball), preLDD(GVMinusMSubGraph, d));
				m = vertexUnion(m, ball);
			} else if (dist[v] > 2 * d) {
				Graph gVMinusM_rev = getSubgraph(g_rev, m, true);
				ArrayList<Integer> ball_rev = volume(gVMinusM_rev, v, i_rnd);
				Graph GVMinusMSubGraph_rev = getSubgraph(gVMinusM_rev, ball_rev, false);
				e_sep = edgeUnion(e_sep, revEdges(layer(gVMinusM_rev, ball_rev)), revEdges(preLDD(GVMinusMSubGraph_rev, d)));
				m = vertexUnion(m, ball_rev);
			} else {
				throw new Exception("RandomTrim failed.");
			}
			
			v = diffVertex(v_far, m, g.v_max);
		}
		
		return e_sep;
	}
	
	// returns the subgraph of g containing only the vertices in ball
	// if setMinus is true, the function returns the subgraph of g containing only the vertices outside of the ball
	public static Graph getSubgraph(Graph g, ArrayList<Integer> ball, boolean setMinus) throws Exception {
		boolean[] contains = new boolean[g.v_max];
		for (int i = 0; i < ball.size(); i++) {
			contains[ball.get(i)] = true;
		}
		
		ArrayList<Integer> vert = new ArrayList<Integer>();
		for (int v : g.vertices) {
			if (!setMinus && contains[v]) {
				vert.add(v);
			} else if (setMinus && !contains[v]) {
				vert.add(v);
			}
		}
		
		Graph subGraph = new Graph(g.v_max, false);
		subGraph.addVertices(vert);
		
		for (int u : vert) {
			ArrayList<Integer> edges = new ArrayList<Integer>();
			ArrayList<Integer> weights = new ArrayList<Integer>();
			
			for (int i = 0; i < g.adjacencyList[u].length; i++) {
				int v = g.adjacencyList[u][i];
				
				if (subGraph.containsVertex[v]) {
					edges.add(v);
					weights.add(g.weights[u][i]);
				}
			}
			
			subGraph.addEdges(u, listToArr(edges), listToArr(weights));
		}
		subGraph.initNullAdjListElts();
		
		return subGraph;
	}
	
	// returns the union of two vertex sets
	public static ArrayList<Integer> vertexUnion(ArrayList<Integer> set1, ArrayList<Integer> set2) {
		Set<Integer> set = new HashSet<Integer>();
		addVerticesToSet(set, set1);
		addVerticesToSet(set, set2);
		
		ArrayList<Integer> output = new ArrayList<Integer>();
		for (int v : set) {
			output.add(v);
		}
		
		return output;
	}
	
	
	public static void addVerticesToSet(Set<Integer> set, ArrayList<Integer> vertices) {
		for (int vertex : vertices) {
			set.add(vertex);
		}
	}
	
	// returns the union of three edge sets
	
	public static ArrayList<int[]> edgeUnion(ArrayList<int[]> set1, ArrayList<int[]> set2, ArrayList<int[]> set3) {
		Set<List<Integer>> set = new HashSet<List<Integer>>();
		addEdgesToSet(set, set1);
		addEdgesToSet(set, set2);
		addEdgesToSet(set, set3);
		
		ArrayList<int[]> distinctEdges = new ArrayList<int[]>();
		for (List<Integer> list : set) {
			int[] edge = {list.get(0), list.get(1)};
			distinctEdges.add(edge);
		}
		
		return distinctEdges;
	}
	
	
	public static void addEdgesToSet(Set<List<Integer>> set, ArrayList<int[]> edges) {
		for (int[] edge : edges) {
			List<Integer> list = new ArrayList<Integer>();
			list.add(edge[0]);
			list.add(edge[1]);
			set.add(list);
		}
	}
	
	// returns a vertex in set1 that's not in set2, or -1 if none exists
	
	public static int diffVertex(ArrayList<Integer> set1, ArrayList<Integer> set2, int v_max) {
		boolean[] contains = new boolean[v_max];
		for (int i = 0; i < set2.size(); i++) {
			contains[set2.get(i)] = true;
		}
		
		for (int i = 0; i < set1.size(); i++) {
			if (!contains[set1.get(i)]) {
				return set1.get(i);
			}
		}
		
		return -1;
	}
	
	//	OUTPUT: a pair (Condition,i) where Condition ??? {1,2,3} and i ??? D is a non-negative integer such that
	//	??? if Condition = 1 then n_G(s,i) > 2n and n_G_rev(s,i) > 2n,
	//	- if Condition = 2 then n_G(s, i) ??? 2n and Vol_G(s, i) and Vol_G(s, i ??? ???D/(3lgn)???) are the same canonical range,
	//	??? if Condition = 3 then n_G_rev(s, i) ??? 2n and Vol_G_rev(s, i) and Vol_G_rev(s, i ??? ???D/(3lgn)???) are in the same canonical range.
	//	If Condition ??? {2, 3} then i ??? D/(3lgn).
	// Runs LayerRange on G and G_rev in parallel.
	
	public static int[] CoreOrLayerRange(Graph g, Graph g_rev, int s, int d) throws Exception {
		ArrayList<int[]> farthestDistancesSeen = new ArrayList<int[]>();
		ArrayList<int[]> farthestDistancesSeen_rev = new ArrayList<int[]>();
		double constant = d / (3.0 * Math.log(g.n));
		boolean[] settled = new boolean[g.v_max];
		boolean[] settled_rev = new boolean[g_rev.v_max];
		int numSettled = 0;
		int numSettled_rev = 0;
	    PriorityQueue<Node> pq = new PriorityQueue<Node>(g.v_max, new Node());
	    PriorityQueue<Node> pq_rev = new PriorityQueue<Node>(g.v_max, new Node());
		int[] dist = new int[g.v_max];
		int[] dist_rev = new int[g.v_max];
		init(g, pq, dist, s);
		init(g_rev, pq_rev, dist_rev, s);
        boolean finished = false;
        boolean finished_rev = false;
        int j = -1;
        int j_rev = -1;
 
        while (true) {
        	if (numSettled == g.n) {
        		finished = true;
        	}
        	if (numSettled_rev == g.n) {
        		finished_rev = true;
        	}
        	
        	if (finished && finished_rev) {
        		// case 1
        		if (j == -1 || j_rev == -1) {
        			int[] output = {2, (int) Math.ceil(constant)};
        			return output;
        		}
        		int[] output = {1, Math.max(j, j_rev)};
        		return output;
        	}
        	
        	if (!finished) {
        		int[] result = oneIterationLayerRange(g, pq, settled, numSettled, farthestDistancesSeen, constant, dist, d);
            	if (result != null) {
            		if (result[0] == 1) {
            			j = result[1];
            			finished = true;
            		} else if (result[0] == 2) {
            			// case 2
            			int[] output = {2, result[1]};
                		return output;
            		}
            	}
            	numSettled++;
        	}
        	
        	if (!finished_rev) {
        		int[] result_rev = oneIterationLayerRange(g_rev, pq_rev, settled_rev, numSettled_rev, farthestDistancesSeen_rev, constant, dist_rev, d);
            	if (result_rev != null) {
            		if (result_rev[0] == 1) {
            			j_rev = result_rev[1];
            			finished_rev = true;
            		} else if (result_rev[0] == 2) {
            			// case 3
            			int[] output = {3, result_rev[1]};
                		return output;
            		}
            	}
            	numSettled_rev++;
        	}
        }
	}	
	
	// returns a copy of g but with edges reversed
	@SuppressWarnings("unchecked")
	public static Graph createGRev(Graph g) throws Exception {
		Graph g_rev = new Graph(g.v_max, false);
		g_rev.addVertices(g.vertices);
		
		ArrayList<Integer>[] edges = new ArrayList[g.v_max];
		ArrayList<Integer>[] weights = new ArrayList[g.v_max];
		for (int i = 0; i < g.v_max; i++) {
			edges[i] = new ArrayList<Integer>();
			weights[i] = new ArrayList<Integer>();
		}
		
		for (int v : g.vertices) {
			for (int i = 0; i < g.adjacencyList[v].length; i++) {
				edges[g.adjacencyList[v][i]].add(v);
				weights[g.adjacencyList[v][i]].add(g.weights[v][i]);
			}
		}
		
		for (int i = 0; i < g.v_max; i++) {
			g_rev.addEdges(i, listToArr(edges[i]), listToArr(weights[i]));
		}
		g_rev.initNullAdjListElts();
		
		return g_rev;
	}
	
	public static int[] listToArr(ArrayList<Integer> list) {
		int[] arr = new int[list.size()];
		for (int i = 0; i < list.size(); i++) {
			arr[i] = list.get(i);
		}
		return arr;
	}
	
	//	OUTPUT: a pair (Condition, i) where Condition ??? {1, 2} and i ??? D is a non-negative integer such that
	//	??? if Condition = 1 then n_G(s,i) > 2n,
	//	??? if Condition = 2 then n_G(s, i) ??? 2n and i ??? D/(3lgn); moreover, Vol_G(s, i) and
	//	Vol_G(s, i ??? ???D/(3lgn)???) are in the same canonical range.
	
//	public static int[] LayerRange(Graph g, int s, int d) throws Exception {
//		// variable i in the paper
//		// elts are int[], where int[0] is the farthest distance seen and int[1] is the value of Vol_G(s, int[0])
//		// last int[0] is the farthest distance seen
//		ArrayList<int[]> farthestDistancesSeen = new ArrayList<int[]>();
//		
//		double constant = d / (3.0 * Math.log(g.n));
//		boolean[] settled = new boolean[g.v_max];
//	    PriorityQueue<Node> pq = new PriorityQueue<Node>(g.v_max, new Node());
//		int[] dist = new int[g.v_max];
//		init(g, pq, dist, s);
// 
//		int numSettled = 0;
//        while (numSettled != g.n) {
//        	int[] result = oneIterationLayerRange(g, pq, settled, numSettled, farthestDistancesSeen, constant, dist, d);
//        	
//        	if (result != null) {
//        		return result;
//        	}
//        	
//        	numSettled++;
//        }
//        
//        throw new Exception("Layer range did not find a satisfying i.");
//	}
	
	
	
	public static int[] oneIterationLayerRange(Graph g, PriorityQueue<Node> pq, boolean[] settled, 
			int numSettled, ArrayList<int[]> farthestDistancesSeen, double constant, int[] dist, int d) throws Exception {
		if (pq.isEmpty()) {
			/*
			 * Nothing left to search.
			 * g is disconnected, since n_G(s,i) <= 2n/3.
			 * Will never be the case that i ??? D/(3lgn).
			 * Return i_big = i + ceil(D/(3lgn)).
			 * Guaranteed that i_big will satisfy i_big >= D/(3lgn) and the canonical ranges.
			 */
			int farthestDistanceSeen = farthestDistancesSeen.get(farthestDistancesSeen.size() - 1)[0];
			int i_big = Math.min(d, farthestDistanceSeen + (int) Math.ceil(constant));
			int[] output = {2, i_big};
			return output;
        }

        int u = pq.remove().node;

        if (settled[u]) {
            return null;
        }

        settled[u] = true;
        
        if (farthestDistancesSeen.size() == 0 || dist[u] > farthestDistancesSeen.get(farthestDistancesSeen.size() - 1)[0]) {
        	int[] seenDistance = {dist[u], numSettled + 1};
        	farthestDistancesSeen.add(seenDistance);
        }
        
        int farthestDistanceSeen = farthestDistancesSeen.get(farthestDistancesSeen.size() - 1)[0];
        
        // case 1
        if (numSettled + 1 > 2.0 * g.n / 3.0) {
        	int[] output = {1, farthestDistanceSeen};
        	return output;
        }
        // case 2
        if ((farthestDistanceSeen >= constant) && sameCanonicalRange(farthestDistancesSeen, constant)) {
        	int[] output = {2, farthestDistanceSeen};
        	return output;
        }
        
        updateNeighbors(g, u, settled, pq, dist, d);

        return null;
	}
	
	// Checks whether Vol_G(s, i - ceil[D/(3logn)]) and Vol_G(s, i) are in the same canonical range.
	// Two numbers are in the same canonical range if they lie in the same half-open interval
	// [2^j, 2^{j+1}), where j is a non-negative integer.
	
	public static boolean sameCanonicalRange(ArrayList<int[]> farthestDistancesSeen, double constant) {
		int i = farthestDistancesSeen.get(farthestDistancesSeen.size() - 1)[0];
		int vol1 = farthestDistancesSeen.get(farthestDistancesSeen.size() - 1)[1];
		
		for (int j = farthestDistancesSeen.size() - 2; j >= 0; j--) {
			if (farthestDistancesSeen.get(j)[0] <= i - Math.ceil(constant)) {
				int vol2 = farthestDistancesSeen.get(j)[1];
				
				// check if vol1 and vol2 are in the same canonical range
				if (Math.floor(Math.log(vol1) / Math.log(2)) == Math.floor(Math.log(vol2) / Math.log(2))) {
					return true;
				}
				break;
			}
		}
		return false;
	}
	
	// {(u,v) in E_H | u in V_H(s,r) and v not in V_H(s,r)}
	
	public static ArrayList<int[]> layer(Graph g, ArrayList<Integer> ball) {
		boolean[] contains = new boolean[g.v_max];
		
		for (int i = 0; i < ball.size(); i++) {
			contains[ball.get(i)] = true;
		}
		
		ArrayList<int[]> edges = new ArrayList<int[]>();
		
		for (int u : ball) {
			for (int v : g.adjacencyList[u]) {
				if (!contains[v]) {
					int[] edge = {u, v};
					edges.add(edge);
				}
			}
		}
		
		return edges;
	}
	
	// returns all the vertices in g within a distance of r from source vertex s using Dijkstra's
	
	public static ArrayList<Integer> volume(Graph g, int s, int r) {
		ArrayList<Integer> output = new ArrayList<Integer>();
		
		boolean[] settled = new boolean[g.v_max];
		int numSettled = 0;
	    PriorityQueue<Node> pq = new PriorityQueue<Node>(g.v_max, new Node());
		int[] dist = new int[g.v_max];
		init(g, pq, dist, s);
 
        while (numSettled != g.n) {
            if (pq.isEmpty()) {
                return output;
            }

            int u = pq.remove().node;

            if (settled[u] || dist[u] > r) {
                continue;
            }
            
            output.add(u);
            settled[u] = true;
            numSettled++;
            
            updateNeighbors(g, u, settled, pq, dist, r);
        }

        return output;
	}
	
	// run Dijkstra's using priority queues in O(V + ElogV)
	
	// returns an array containing the distances from s to every vertex in g
	// runs in O(V + ElogV)
	public static int[] Dijkstra(Graph g, int s) {		
		boolean[] settled = new boolean[g.v_max];
		int numSettled = 0;
	    PriorityQueue<Node> pq = new PriorityQueue<Node>(g.v_max, new Node());
		int[] dist = new int[g.v_max];
		init(g, pq, dist, s);
 
        while (numSettled != g.n) {
            if (pq.isEmpty()) {
                return dist;
            }

            int u = pq.remove().node;

            if (settled[u]) {
                continue;
            }

            settled[u] = true;
            numSettled++;
            
            updateNeighbors(g, u, settled, pq, dist, Integer.MAX_VALUE);
        }
        
        return dist;
	}
	
	
	public static void init(Graph g, PriorityQueue<Node> pq, int[] dist, int s) {
		for (int i = 0; i < g.v_max; i++) {
            dist[i] = Integer.MAX_VALUE;
        }
 
        pq.add(new Node(s, 0));
        dist[s] = 0;
	}
	
	
	public static void updateNeighbors(Graph g, int u, boolean[] settled, PriorityQueue<Node> pq, int[] dist, int d) {
        for (int i = 0; i < g.adjacencyList[u].length; i++) {
        	int v = g.adjacencyList[u][i];
        	
            if (!settled[v]) {
                int newDistance = dist[u] + g.weights[u][i];

                if (newDistance < dist[v]) {
                    dist[v] = newDistance;
                }
                
                // only want to process nodes within a distance of d from the source
                if (dist[v] <= d) {
                	pq.add(new Node(v, dist[v]));
                }
            }
        }
	}
}

class Node implements Comparator<Node> {
    public int node;
    public int cost;
 
    public Node() {}
 
    public Node(int node, int cost) {
        this.node = node;
        this.cost = cost;
    }
    
    @Override public boolean equals(Object o) {
    	Node node2 = (Node) o;
    	if ((node == node2.node) && (cost == node2.cost)) {
    		return true;
    	}
    	return false;
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