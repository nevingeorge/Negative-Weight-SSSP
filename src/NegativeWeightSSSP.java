import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.Stack;

public class NegativeWeightSSSP {
	
	public static double startTime; // in ms
	public static final double maxRunTime = 5000; // in seconds

	public static void main(String[] args) throws Exception {
//		BufferedReader f = new BufferedReader(new FileReader("USA-small"));
//		
//		Graph g = new Graph(830, true);
//		
//		String line = f.readLine();
//		while (line != null) {
//			String[] arr = line.split(" ");
//			if (Math.random() < .1) {
//				g.addEdge(Integer.parseInt(arr[1]), Integer.parseInt(arr[2]), -1);
//			} else {
//				g.addEdge(Integer.parseInt(arr[1]), Integer.parseInt(arr[2]), Integer.parseInt(arr[3]));
//			}
//			line = f.readLine();
//		}
//		f.close();
		
//		int[] tree = bitScaling(g, 0, false);
//
//		for (int i = 0; i < 533; i++) {
//			System.out.println("Parent of vertex " + i + ": " + tree[i]);
//		}

		Graph g = new Graph(5, true);
		g.addEdge(0, 1, 2);
		g.addEdge(1, 2, 2);
		g.addEdge(0, 3, 2);
		g.addEdge(3, 4, -2);
		g.addEdge(4, 2, 2);
		int[] tree = bitScaling(g, 0, false);
		
		for (int i = 0; i < 5; i++) {
			System.out.println("Parent of vertex " + i + ": " + tree[i]);
		}
	}
	
	// Runs the bit scaling algorithm of Goldberg and Rao to return a shortest path tree for g_in
	public static int[] bitScaling(Graph g_in, int s, boolean constantOutDegree) throws Exception {
		startTime = System.currentTimeMillis();
		
		Graph g;
		if (constantOutDegree) {
			g = constantOutDegree(g_in);
		} else {
			g = g_in;
		}
		
		int minWeight = Integer.MAX_VALUE;
		for (int u : g.vertices) {
			for (int v : g.adjacencyList[u]) {
				int weight = g.weights.get(Graph.edgeToString(u, v));
				if (weight < minWeight) {
					minWeight = weight;
				}
			}
		}
		
		if (minWeight >= 0) {
			return getShortestPathTree(g_in, s);
		}
		
		int precision = (int) Math.pow(2, (int) logBase2(-1 * minWeight));
		HashMap<Integer, Integer> phi = new HashMap<Integer, Integer>();
		
		while (precision >= 1) {
			Graph gScaledS = new Graph(g.n + 1, true); // scaled graph with added dummy source vertex s
			
			for (int u : g.vertices) {
				for (int v : g.adjacencyList[u]) {
					int phiEdgeValue = 0;
					if (phi.containsKey(u)) {
						phiEdgeValue += phi.get(u);
					} 
					if (phi.containsKey(v)) {
						phiEdgeValue -= phi.get(v);
					}
					
					double weight = g.weights.get(Graph.edgeToString(u, v));
					int roundedWeight = phiEdgeValue + (int) Math.ceil(weight / precision);
					if (roundedWeight < -1) {
						throw new Exception("Bit scaling produced an edge of weight less than -1.");
					}
					
					gScaledS.addEdge(u, v, roundedWeight);
				}
			}
			
			for (int v = 0; v < g.n; v++) {
				gScaledS.addEdge(g.n, v, 0);
			}
			
			int[] tree = SPmain(gScaledS, g.n);
			
			int[] dist = new int[g.n + 1]; // will store the shortest distances from the dummy vertex s to every other vertex
			getDistances(gScaledS, tree, dist, 0, g.n, new boolean[g.n + 1]);
			
			for (int u: g.vertices) {
				// add then multiply by 2
				int phiValue = 0;
				if (phi.containsKey(u)) {
					phiValue = phi.get(u);
				} 
				
				if (precision == 1) {
					phi.put(u, phiValue + dist[u]);
				} else {
					phi.put(u, 2 * (phiValue + dist[u]));
				}
			}
			
			precision /= 2;
		}
		
		for (int u : g.vertices) {
			for (int v : g.adjacencyList[u]) {
				int initWeight = g.weights.get(Graph.edgeToString(u, v));
				g.updateEdgeWeight(u, v, initWeight + phi.get(u) - phi.get(v));
			}
		}
		
		int[] tree = getShortestPathTree(g, s);
		if (constantOutDegree) {
			return treeForInputG(g_in, s, tree);
		}
		return tree;
	}
	
	public static void getDistances(Graph g, int[] tree, int[] dist, int curDis, int curVertex, boolean[] visited) {
		if (!visited[curVertex]) {
			dist[curVertex] = curDis;
			visited[curVertex] = true;
			
			for (int v : g.adjacencyList[curVertex]) {
				if (tree[v] == curVertex) {
					int weight = g.weights.get(Graph.edgeToString(curVertex, v));
					getDistances(g, tree, dist, curDis + weight, v, visited);
				}
			}
		}
	}
	
	/*
	 * Returns the shortest path tree in g from s.
	 * The boolean constantOutDegree is false if we wish to ignore the assumption that g has constant out-degree.
	 */
	public static int[] SPmain(Graph g_in, int s) throws Exception {		
		Graph g = new Graph(g_in.n, true);
		
		for (int u : g_in.vertices) {
			for (int v : g_in.adjacencyList[u]) {
				int initWeight = g_in.weights.get(Graph.edgeToString(u, v));
				g.addEdge(u, v, initWeight * 2 * g_in.n);
			}
		}
		int B = roundPower2(2 * g.n);
		HashMap<Integer, Integer> phi = new HashMap<Integer, Integer>();
		
		for (int i = 1; i <= logBase2(B); i++) {
			Graph g_phi = createModifiedGB(g, 0, false, new HashSet<int[]>(), phi);
			HashMap<Integer, Integer> phi_i = ScaleDown(g_phi, g.n, B / (int) Math.pow(2, i));
			phi = addPhi(phi, phi_i);
		}
		
		// create G^*
		for (int u : g.vertices) {
			for (int v : g.adjacencyList[u]) {
				int initWeight = g.weights.get(Graph.edgeToString(u, v));
				g.updateEdgeWeight(u, v, initWeight + phi.get(u) - phi.get(v) + 1);
			}
		}
		
		return getShortestPathTree(g, s);
	}
	
	/*
	 * Outputs a graph that has constant out degree 2.
	 * For each vertex v (letting OD mean out degree),
	 * if OD(v) = 0, add a vertex v'. make a loop between v and v', as well as self-loops for v and v'.
	 * if OD(v) = 1, add a self-loop.
	 * if OD(v) >= 2, add OD(v) - 1 vertices. Turn the OD(v) vertices into a cycle.
	 * Make each of the OD(v) out edges come out of one of the vertices in the cycle.
	 * Make all of the input edges going into the original vertices.
	 * All the newly added edges have weight 0.
	 */
	public static Graph constantOutDegree(Graph g) throws Exception {
		int numVertices = 0;
		int[] indexOfMainVertex = new int[g.n]; // contains the indices of the original vertices in g
		for (int v = 0; v < g.n; v++) {
			indexOfMainVertex[v] = numVertices;
			int outDegree = g.adjacencyList[v].size();
			if (outDegree == 0) {
				numVertices += 2;
			} else {
				numVertices += outDegree;
			}
		}
		
		Graph gOut = new Graph(numVertices, true);
		
		for (int u = 0; u < g.n; u++) {
			int outDegree = g.adjacencyList[u].size();
			if (outDegree == 0) {
				// self-loops
				gOut.addEdge(indexOfMainVertex[u], indexOfMainVertex[u], 0);
				gOut.addEdge(indexOfMainVertex[u] + 1, indexOfMainVertex[u] + 1, 0);
				
				// loop between v and v'
				gOut.addEdge(indexOfMainVertex[u], indexOfMainVertex[u] + 1, 0);
				gOut.addEdge(indexOfMainVertex[u] + 1, indexOfMainVertex[u], 0);
			} else if (outDegree == 1) {
				int v = g.adjacencyList[u].get(0);
				// original edge
				gOut.addEdge(indexOfMainVertex[u], indexOfMainVertex[v], 0);
				
				// create a self-loop
				gOut.addEdge(indexOfMainVertex[u], indexOfMainVertex[u], 0);
			} else {
				int i = 0;
				for (int v : g.adjacencyList[u]) {
					// edge (u, v) => (u + i, indexOfMainVertex[v]), i in [0, g.adjacencyList[u].size())
					gOut.addEdge(indexOfMainVertex[u] + i, indexOfMainVertex[v], 0);
					
					// create a cycle
					if (i == outDegree - 1) {
						gOut.addEdge(indexOfMainVertex[u] + i, indexOfMainVertex[u], 0);
					} else {
						gOut.addEdge(indexOfMainVertex[u] + i, indexOfMainVertex[u] + i + 1, 0);
					}
					
					i++;
				}
			}
		}
		
		return gOut;
	}
	
	// converts the tree for gOut into a tree for g_in
	public static int[] treeForInputG(Graph g, int s, int[] tree) {
		HashMap<Integer, Integer> newVerttoOldVert = new HashMap<Integer, Integer>();
		int numVertices = 0;
		HashSet<Integer> originalVertices = new HashSet<Integer>();
		for (int v = 0; v < g.n; v++) {
			originalVertices.add(numVertices);
			int outDegree = g.adjacencyList[v].size();
			
			if (outDegree == 0) {
				newVerttoOldVert.put(numVertices, v);
				newVerttoOldVert.put(numVertices + 1, v);
				numVertices += 2;
			} else {
				for (int i = 0; i < outDegree; i++) {
					newVerttoOldVert.put(numVertices + i, v);
				}
				numVertices += outDegree;
			}
		}
		
		int[] treeInputG = new int[g.n];
		for (int v = 0; v < numVertices; v++) {
			if (originalVertices.contains(v)) {
				if (tree[v] == -1) {
					treeInputG[newVerttoOldVert.get(v)] = -1;
				} else {
					treeInputG[newVerttoOldVert.get(v)] = newVerttoOldVert.get(tree[v]);
				}
			}
		}
		
		return treeInputG;
	}
	
	// rounds n up to the nearest power of 2
	public static int roundPower2(int n) {
		return (int) Math.pow(2, Math.ceil(logBase2(n)));
	}
	
	public static double logBase2(int n) {
		return Math.log(n) / Math.log(2);
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
		if (System.currentTimeMillis() - startTime >= maxRunTime * Math.pow(10, 3)) {
			throw new Exception("Exceeded total allotted run time.");
		}
		
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
					int weight = g.weights.get(Graph.edgeToString(u, v));
					
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
				int edgeWeight = g.weights.get(Graph.edgeToString(u, v));
				
				if ((SCCu != SCCv) && edgeWeight < mu[topOrdering.get(SCCv)]) {
					mu[topOrdering.get(SCCv)] = edgeWeight;
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
		PriorityQueue<Node> pq = new PriorityQueue<Node>(Gs.v_max, new Node());
		pq.add(new Node(s, dist[s]));
		HashSet<Integer> marked = new HashSet<Integer>();
		
		while (!pq.isEmpty()) {
			// Dijkstra Phase
			while (!pq.isEmpty()) {
				int v = pq.remove().node;
				marked.add(v);
				
				for (int x : Gs.adjacencyList[v]) {
					int edgeWeight = Gs.weights.get(Graph.edgeToString(v, x));
					
					if (edgeWeight >= 0 && (dist[v] + edgeWeight < dist[x])) {
						marked.add(v);
						pq.add(new Node(x, dist[x]));
						dist[x] = dist[v] + edgeWeight;
					}
				}
			}
			
			// Bellman-Ford Phase
			for (int v : marked) {
				for (int x : Gs.adjacencyList[v]) {
					int edgeWeight = Gs.weights.get(Graph.edgeToString(v, x));
					
					if (edgeWeight < 0 && (dist[v] + edgeWeight < dist[x])) {
						dist[x] = dist[v] + edgeWeight;
						pq.add(new Node(x, dist[x]));
					}
				}
			}
			marked = new HashSet<Integer>();
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
				Gs.addEdge(u, v, g.weights.get(Graph.edgeToString(u, v)));
			}
			Gs.addEdge(s, u, 0);
		}
		
		return Gs;
	}
	
	/*
	 * Let the output be int[] out. 
	 * For vertex i, out[i] = parent of vertex i in the shortest path tree in g from s.
	 * out[i] = -1 indicates that vertex i has no parent in the tree.
	 */
	public static int[] getShortestPathTree(Graph g, int s) {		
		Set<Integer> settled = new HashSet<Integer>();
	    PriorityQueue<Node> pq = new PriorityQueue<Node>(g.v_max, new Node());
		int[] dist = new int[g.v_max];
		int[] tree = new int[g.v_max];
		
		for (int i = 0; i < g.v_max; i++) {
            dist[i] = Integer.MAX_VALUE;
            tree[i] = -1;
        }
 
        pq.add(new Node(s, 0));
        dist[s] = 0;
 
        while (settled.size() != g.n) {
            if (pq.isEmpty()) {
                return tree;
            }

            int u = pq.remove().node;

            if (settled.contains(u)) {
                continue;
            }

            settled.add(u);
            updateTreeNeighbors(g, u, tree, settled, pq, dist);
        }
        
        return tree;
	}
	
	
	public static void updateTreeNeighbors(Graph g, int u, int[] tree, Set<Integer> settled, PriorityQueue<Node> pq, int[] dist) {
		ArrayList<Integer> neighbors = g.adjacencyList[u];

        for (int v : neighbors) {
            if (!settled.contains(v)) {
                int newDistance = dist[u] + g.weights.get(Graph.edgeToString(u, v));

                if (newDistance < dist[v]) {
                    dist[v] = newDistance;
                    tree[v] = u;
                }
                
                pq.add(new Node(v, dist[v]));
            }
        }
	}
}
