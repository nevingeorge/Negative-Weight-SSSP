import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

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
	
//	1. INPUT REQUIREMENTS:
//		(a) B is positive integer, w is integral, and w(e) ≥ −2B for all e ∈ E
//		(b) If the graph G does not contain a negative-weight cycle then the input 
//			must satisfy η(GB) ≤ ∆; that is, for every v ∈ V there is a shortest 
//			sv-path in GBs with at most ∆ negative edges
//		(c) All vertices in G have constant out-degree
//	2. OUTPUT: If it terminates, the algorithm returns an integral price function φ 
//		such that wφ(e) ≥ −B for all e ∈ E
//	3. RUNNING TIME: If G does not contain a negative-weight cycle, then the 
//		algorithm has epected runtime O(m log3(n) log(∆)). 
//		Remark: If G contains a negative-weight cycle, there is no guarantee 
//		on the runtime, and the algorithm might not even terminate; but if the 
//		algorithm does terminate, it always produces a correct output.
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
			Graph H = createModifiedGB(g, 0, false, getEdgesBetweenSCCs(g, SCCs), new HashMap<Integer, Integer>());
			HashMap<Integer, Integer> phi_1 = ScaleDown(H, delta / 2, B);
			
			// phase 2
			Graph g_B_E_sep_phi1 = createModifiedGB(g, B, false, E_sep_hash, phi_1);
			HashMap<Integer, Integer> phi = FixDAGEdges(g_B_E_sep_phi1, SCCs);
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
		Graph modG = new Graph(g.n, false);
		modG.addVertices(g.vertices);
		
		for (int u = 0; u < g.n; u++) {
			if (g.containsVertex[u]) { 
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
		}
		
		return modG;
	}	
	public static HashSet<int[]> getEdgesBetweenSCCs(Graph g, ArrayList<ArrayList<Integer>> SCCs) {
		HashMap<Integer, Integer> vertexToSCCIndex = new HashMap<Integer, Integer>();
		for (int i = 0; i < SCCs.size(); i++) {
			for (int v : SCCs.get(i)) {
				vertexToSCCIndex.put(v, i);
			}
		}
		
		HashSet<int[]> edgesBetweenSCCs = new HashSet<int[]>();
		for (int u = 0; u < g.n; u++) {
			if (g.containsVertex[u]) { 
				for (int v : g.adjacencyList[u]) {
					if (vertexToSCCIndex.get(u) != vertexToSCCIndex.get(v)) {
						int[] edge = {u, v};
						edgesBetweenSCCs.add(edge);
					}
				}
			}
		}
		
		return edgesBetweenSCCs;
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
	
	public static HashMap<Integer, Integer> FixDAGEdges(Graph g, ArrayList<ArrayList<Integer>> SCCs) {
		return null;
	}
	
	public static HashMap<Integer, Integer> ElimNeg(Graph g) {
		return null;
	}
	
	// returns the subgraph of g containing only the vertices in ball
	// if setMinus is true, the function returns the subgraph of g containing only the vertices outside of the ball
	public static Graph getSubgraph(Graph g, ArrayList<Integer> ball, boolean setMinus) throws Exception {
		boolean[] contains = new boolean[g.v_max];
		for (int i = 0; i < ball.size(); i++) {
			contains[ball.get(i)] = true;
		}
		
		ArrayList<Integer> vert = new ArrayList<Integer>();
		for (int i = 0; i < g.v_max; i++) {
			if (g.containsVertex[i]) {
				if (!setMinus && contains[i]) {
					vert.add(i);
				} else if (setMinus && !contains[i]) {
					vert.add(i);
				}
			}
		}
		
		Graph subGraph = new Graph(g.v_max, false);
		subGraph.addVertices(vert);
		
		for (int u : vert) {
			ArrayList<Integer> edgesU = g.adjacencyList[u];
			
			for (int v : edgesU) {
				if (subGraph.containsVertex[v]) {
					subGraph.addEdge(u, v, g.weights[u][v]);
				}
			}
		}
		
		return subGraph;
	}
}
