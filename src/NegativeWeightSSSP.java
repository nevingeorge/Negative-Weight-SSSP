import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.Stack;

public class NegativeWeightSSSP {

	public static void main(String[] args) throws Exception {
//	    Graph g1 = new Graph(5, true);
//	 
//	    g1.addEdge(1, 0, 1);
//	    g1.addEdge(0, 2, 1);
//	    g1.addEdge(2, 1, 1);
//	    g1.addEdge(0, 3, 1);
//	    g1.addEdge(3, 4, 1);
//	    System.out.println("SCC in first graph ");
//	    ArrayList<ArrayList<Integer>> out1 = g1.SCC();
//	    for (ArrayList<Integer> v : out1) {
//	    	for (int u : v) {
//	    		System.out.print(u + " ");
//	    	}
//	    	System.out.println();
//	    }
//	    
////	    ArrayList<int[]> edges = LowDiameterDecomposition.LDD(g1, 0);
////		
////		for (int[] edge : edges) {
////			System.out.println(edge[0] + " " + edge[1]);
////		}
//	 
//	    Graph g2 = new Graph(4, true);
//	    g2.addEdge(0, 1, 1);
//	    g2.addEdge(1, 2, 1);
//	    g2.addEdge(2, 3, 1);
//	    System.out.println("\nSCC in second graph ");
//	    ArrayList<ArrayList<Integer>> out2 = g2.SCC();
//	    for (ArrayList<Integer> v : out2) {
//	    	for (int u : v) {
//	    		System.out.print(u + " ");
//	    	}
//	    	System.out.println();
//	    }
//	    
//	    Graph g4 = new Graph(11, true);
//	    g4.addEdge(0, 1, 1);
//	    g4.addEdge(0, 3, 1);
//	    g4.addEdge(1, 2, 1);
//	    g4.addEdge(1, 4, 1);
//	    g4.addEdge(2, 0, 1);
//	    g4.addEdge(2, 6, 1);
//	    g4.addEdge(3, 2, 1);
//	    g4.addEdge(4, 5, 1);
//	    g4.addEdge(4, 6, 1);
//	    g4.addEdge(5, 6, 1);
//	    g4.addEdge(5, 7, 1);
//	    g4.addEdge(5, 8, 1);
//	    g4.addEdge(5, 9, 1);
//	    g4.addEdge(6, 4, 1);
//	    g4.addEdge(7, 9, 1);
//	    g4.addEdge(8, 9, 1);
//	    g4.addEdge(9, 8, 1);
//	    System.out.println("\nSCC in fourth graph ");
//	    ArrayList<ArrayList<Integer>> out4 = g4.SCC();
//	    for (ArrayList<Integer> v : out4) {
//	    	for (int u : v) {
//	    		System.out.print(u + " ");
//	    	}
//	    	System.out.println();
//	    }
	}
	
	/*
	 * 1. INPUT REQUIREMENTS:
	 * (a) B is positive integer, w is integral, and w(e) ≥ −2B for all e ∈ E
	 * (b) If the graph G does not contain a negative-weight cycle then the input 
	 * 		must satisfy η(GB) ≤ ∆; that is, for every v ∈ V there is a shortest 
	 * 		sv-path in GBs with at most ∆ negative edges
	 * (c) All vertices in G have constant out-degree
	 * 2. OUTPUT: If it terminates, the algorithm returns an integral price function φ 
	 * 	such that wφ(e) ≥ −B for all e ∈ E
	 * 3. RUNNING TIME: If G does not contain a negative-weight cycle, then the 
	 * 	algorithm has epected runtime O(m log3(n) log(∆)). 
	 * 	Remark: If G contains a negative-weight cycle, there is no guarantee 
	 * 	on the runtime, and the algorithm might not even terminate; but if the 
	 * 	algorithm does terminate, it always produces a correct output.
	 */
	public static HashMap<Integer, Integer> ScaleDown(Graph g, int delta, int B) throws Exception {
		// if phi(x) == null, assume that phi(x) = 0
		HashMap<Integer, Integer> phi_2 = new HashMap<Integer, Integer>();
		
		if (delta > 2) { // base case
			double d = delta / 2.0;
			Graph g_B_nneg = createModifiedGB(g, B, true, new HashSet<int[]>(), new HashMap<Integer, Integer>());
			
			// phase 0
			ArrayList<int[]> E_sep = LowDiameterDecomposition.LDD(g_B_nneg, (int) (d * B));
			HashSet<int[]> E_sep_hash = new HashSet<int[]>(E_sep);
			Graph g_B_Esep = createModifiedGB(g, B, false, E_sep_hash, new HashMap<Integer, Integer>());
			ArrayList<ArrayList<Integer>> SCCs = g_B_Esep.SCC();
			
			// phase 1
			HashMap<Integer, Integer> vertexToSCCMap = getVertexToSCCMap(SCCs);
			HashSet<int[]> edgesBetweenSCCs = getEdgesBetweenSCCs(g, vertexToSCCMap);
			Graph H = createModifiedGB(g, 0, false, edgesBetweenSCCs, new HashMap<Integer, Integer>());
			HashMap<Integer, Integer> phi_1 = ScaleDown(H, delta / 2, B);
			
			// phase 2
			Graph g_B_E_sep_phi1 = createModifiedGB(g, B, false, E_sep_hash, phi_1);
			HashMap<Integer, Integer> phi = FixDAGEdges(g_B_E_sep_phi1, SCCs, vertexToSCCMap, edgesBetweenSCCs);
			phi_2 = addPhi(phi_1, phi);
		}
		
		// phase 3
		Graph g_B_phi2 = createModifiedGB(g, B, false, new HashSet<int[]>(), phi_2);
		HashMap<Integer, Integer> phi_prime = ElimNeg(g_B_phi2);
		HashMap<Integer, Integer> phi_3 = addPhi(phi_2, phi_prime);
		
		return phi_3;
	}
	
	/*
	 * Creates G^B_phi = (V, E, w^B_phi), where 
	 * w^B_phi(e) = w(e) + phi(u) - phi(v) if w(e) >= 0, 
	 * and w(e) + B + phi(u) - phi(v) if w(e) < 0.
	 * If nneg == true, w(e) = max{0, w^B(e)}.
	 * Removes all the edges in remEdges.
	 */
	public static Graph createModifiedGB(Graph g, int B, boolean nneg, HashSet<int[]> remEdges, HashMap<Integer, Integer> phi) throws Exception {
		Graph modG = new Graph(g.v_max, false);
		modG.addVertices(g.vertices);
		
		for (int u : g.vertices) {
			for (int v : g.adjacencyList[u]) {
				int[] edge = {u, v};
				if (!remEdges.contains(edge)) {
					double weight = g.weights[u][v];
					
					if (weight < 0) {
						weight += B;
					}
					
					if (nneg) {
						weight = Math.max(0, weight);
					}
					
					if (phi.get(u) != null) {
						weight += phi.get(u);
					}
					if (phi.get(v) != null) {
						weight -= phi.get(v);
					}
					
					modG.addEdge(u, v, weight);
				}
			}
		}
		
		return modG;
	}	
	
	public static HashSet<int[]> getEdgesBetweenSCCs(Graph g, HashMap<Integer, Integer> vertexToSCCMap) {
		HashSet<int[]> edgesBetweenSCCs = new HashSet<int[]>();
		for (int u : g.vertices) {
			for (int v : g.adjacencyList[u]) {
				if (vertexToSCCMap.get(u) != vertexToSCCMap.get(v)) {
					int[] edge = {u, v};
					edgesBetweenSCCs.add(edge);
				}
			}
		}
		
		return edgesBetweenSCCs;
	}
	
	public static HashMap<Integer, Integer> getVertexToSCCMap(ArrayList<ArrayList<Integer>> SCCs) {
		HashMap<Integer, Integer> vertexToSCCMap = new HashMap<Integer, Integer>();
		for (int i = 0; i < SCCs.size(); i++) {
			for (int v : SCCs.get(i)) {
				vertexToSCCMap.put(v, i);
			}
		}
		
		return vertexToSCCMap;
	}

	public static HashMap<Integer, Integer> addPhi(HashMap<Integer, Integer> phi_1, HashMap<Integer, Integer> phi_2) {
		HashMap<Integer, Integer> newPhi = new HashMap<Integer, Integer>();
		
		for (int v : phi_1.keySet()) {
			if (phi_2.containsKey(v)) {
				newPhi.put(v, phi_1.get(v) + phi_2.get(v));
			} else {
				newPhi.put(v, phi_1.get(v));
			}
		}
		
		for (int v : phi_2.keySet()) {
			if (!newPhi.containsKey(v)) {
				// phi_1 does not contain v
				newPhi.put(v, phi_2.get(v));
			}
		}
		
		return newPhi;
	}
	
	public static HashMap<Integer, Integer> FixDAGEdges(Graph g, ArrayList<ArrayList<Integer>> SCCs, HashMap<Integer, Integer> vertexToSCCMap, HashSet<int[]> edgesBetweenSCCs) {
		HashMap<Integer, Integer> topOrdering = topSort(SCCs.size(), createSCCAdjList(SCCs, vertexToSCCMap, edgesBetweenSCCs));
		
		int[] mu = new int[SCCs.size()]; // indices are in topological order (e.g., index 0 corresponds to the first SCC in topological order)
		for (int u : g.vertices) {
			for (int v : g.adjacencyList[u]) {
				int SCCu = vertexToSCCMap.get(u);
				int SCCv = vertexToSCCMap.get(v);
				if ((SCCu != SCCv) && g.weights[u][v] < mu[topOrdering.get(SCCv)]) {
					mu[topOrdering.get(SCCv)] = (int) g.weights[u][v];
				}
			}
		}
		
		HashMap<Integer, Integer> phi = new HashMap<Integer, Integer>();
		int m = 0;
		for (int j = 1; j < SCCs.size(); j++) {
			m += mu[j];
			for (int v : SCCs.get(topOrdering.get(-1 * j - 1))) {
				phi.put(v, m);
			}
		}
		
		return phi;
	}
	
	// returns the adjacency list for the DAG where every SCC is viewed as a single vertex
	public static ArrayList<Integer>[] createSCCAdjList(ArrayList<ArrayList<Integer>> SCCs, HashMap<Integer, Integer> vertexToSCCMap, HashSet<int[]> edgesBetweenSCCs) {
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] SCCAdjList = new ArrayList[SCCs.size()];
		for (int i = 0; i < SCCs.size(); i++) {
			SCCAdjList[i] = new ArrayList<Integer>();
		}
		boolean[][] containsEdge = new boolean[SCCs.size()][SCCs.size()];
		
		for (int[] edge : edgesBetweenSCCs) {
			int u = vertexToSCCMap.get(edge[0]);
			int v = vertexToSCCMap.get(edge[1]);
			if (!containsEdge[u][v]) {
				containsEdge[u][v] = true;
				SCCAdjList[u].add(v);
			}
		}
		
		return SCCAdjList;
	}
	
	
	/*
	 * Input: DAG
	 * Calculates a topological ordering of the n vertices in the DAG
	 * Returns a map of vertex v to its index i in the ordering
	 * The index in the order i is also mapped to the vertex v (bijective mapping);
	 * however, to prevent duplicate keys, instead of using the key i we use the
	 * key -1 * i - 1.
	 */
	public static HashMap<Integer, Integer> topSort(int n, ArrayList<Integer>[] adjList) {
		HashMap<Integer, Integer> topOrdering = new HashMap<Integer, Integer>();
		Stack<Integer> stack = new Stack<Integer>();
		boolean[] visited = new boolean[n];
        
        for (int i = 0; i < n; i++) {
            if (!visited[i]) {
                topSortUtil(i, visited, stack, adjList);
            }
        }
        
        int i = 0;
        while (!stack.empty()) {
        	int v = stack.pop();
        	topOrdering.put(v, i);
        	topOrdering.put(-1 * i - 1, v);
        	i++;
        }
        
        return topOrdering;
	}
	
	public static void topSortUtil(int u, boolean[] visited, Stack<Integer> stack, ArrayList<Integer>[] adjList) {
		visited[u] = true;
		
		for (int v : adjList[u]) {
			if (!visited[v]) {
				topSortUtil(v, visited, stack, adjList);
			}
		}

		stack.push(u);
	}
	
	
	/*
	 * ElimNeg takes as input a graph G = (V,E,w) in which all vertices have constant out-degree. 
	 * The algorithm outputs a price function φ such that w_φ(e) ≥ 0 for all e ∈ E 
	 * and has running time O(log(n)*(n + sum_v∈V ηG(v))); 
	 * note that if G contains a negative-weight cycle then v∈V ηG(v) = ∞ 
	 * so the algorithm will never terminate and hence not produce any output.
	 */
	public static HashMap<Integer, Integer> ElimNeg(Graph g) throws Exception {
		Graph Gs = createGs(g);
		int[] dist = new int[Gs.v_max];
		int s = Gs.v_max - 1;
		dist[s] = 0;
		for (int v = 0; v < s; v++) {
			dist[v] = Integer.MAX_VALUE;
		}
		PriorityQueue<Node> pq = new PriorityQueue<Node>(g.v_max, new Node());
		pq.add(new Node(s, dist[s]));
		HashSet<Integer> marked = new HashSet<Integer>();
		
		while (!pq.isEmpty()) {
			// Dijkstra Phase
			while (!pq.isEmpty()) {
				int v = pq.remove().node;
				marked.add(v);
				
				for (int x : g.adjacencyList[v]) {
					if (g.weights[v][x] >= 0 && (dist[v] + g.weights[v][x] < dist[x])) {
						marked.add(v);
						pq.add(new Node(x, dist[x]));
						dist[x] = dist[v] + (int) g.weights[v][x];
					}
				}
			}
			
			// Bellman-Ford Phase
			for (int v : marked) {
				for (int x : g.adjacencyList[v]) {
					if (g.weights[v][x] < 0 && (dist[v] + g.weights[v][x] < dist[x])) {
						dist[x] = dist[v] + (int) g.weights[v][x];
						pq.add(new Node(x, dist[x]));
					}
				}
				marked.remove(v);
			}
		}
		
		HashMap<Integer, Integer> phi = new HashMap<Integer, Integer>();
		for (int v : g.vertices) {
			phi.put(v, dist[v]);
		}
		return phi;
	}
	
	/*
	 * Returns the graph that is g with an added dummy vertex s 
	 * and edges of weight 0 connecting s to every vertex in G.
	 */
	public static Graph createGs(Graph g) throws Exception {
		int s = g.v_max;
		Graph Gs = new Graph(g.v_max + 1, false);
		Gs.addVertices(g.vertices);
		Gs.addVertex(s);
		
		for (int u : g.vertices) {
			for (int v : g.adjacencyList[u]) {
				Gs.addEdge(u, v, g.weights[u][v]);
			}
			Gs.addEdge(s, u, 0);
		}
		
		return Gs;
	}
}
