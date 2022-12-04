import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.Stack;

public class NegativeWeightSSSP {
	
	public static double startTime; // in ms
	public static final boolean CHECKS = true;

	public static void main(String[] args) throws Exception {
	    String fileName = "USA-small";
		int src = 1;
		
		Graph g_in = readInput(fileName);
		Graph g = getConnectedSubgraph(g_in, src);

		while (g.hasNegCycle() || g.hasNoNegativeEdgeWeights()) {
			g_in = readInput(fileName);
			g = getConnectedSubgraph(g_in, src);
		}
		
		System.out.println("Obtained graph.");
		
		runBellmanFord(g);
		
		int[] tree = bitScaling(g, src);
		
		System.out.println();
		for (int v : g.vertices) {
			System.out.println("Parent of vertex " + v + ": " + tree[v]);
		}
	}
	
	@SuppressWarnings("unchecked")
	public static Graph readInput(String fileName) throws Exception {
		BufferedReader f = new BufferedReader(new FileReader(fileName));
		
		int g_size = Integer.parseInt(f.readLine());
		Graph g = new Graph(g_size, true);
		ArrayList<Integer>[] edges = new ArrayList[g_size];
		ArrayList<Integer>[] weights = new ArrayList[g_size];
		
		for (int i = 0; i < g_size; i++) {
			edges[i] = new ArrayList<Integer>();
			weights[i] = new ArrayList<Integer>();
		}
		
		boolean[][] edge_exists = new boolean[g.v_max][g.v_max];
		
		String line = f.readLine();
		while (line != null) {
			String[] arr = line.split(" ");
			
			int u = Integer.parseInt(arr[1]);
			int v = Integer.parseInt(arr[2]);
			int weight = Integer.parseInt(arr[3]);
			
			if (!edge_exists[u][v]) {
				edges[u].add(v);
				if (Math.random() < .05) {
					// weights[u].add(-1 * weight);
					weights[u].add(-1);
				} else {
					weights[u].add(weight);
				}
				// weights[u].add(weight);
				
				edge_exists[u][v] = true;
			}
			
			line = f.readLine();
		}
		f.close();
		
		for (int i = 0; i < g_size; i++) {
			g.addEdges(i, LowDiameterDecomposition.listToArr(edges[i]), LowDiameterDecomposition.listToArr(weights[i]));
		}
		g.initNullAdjListElts();
		
		return g;
	}
	
	public static Graph getConnectedSubgraph(Graph g, int src) throws Exception {
		boolean[] reachable = new boolean[g.v_max];
		findReachable(g, src, reachable);
		
		Graph subGraph = new Graph(g.v_max, false);
		for (int v = 0; v < g.v_max; v++) {
			if (reachable[v]) {
				subGraph.addVertex(v);
			}
		}
		
		for (int u : subGraph.vertices) {
			ArrayList<Integer> outVertices = new ArrayList<Integer>();
			ArrayList<Integer> weights = new ArrayList<Integer>();
			
			for (int i = 0; i < g.adjacencyList[u].length; i++) {				
				if (reachable[g.adjacencyList[u][i]]) {
					outVertices.add(g.adjacencyList[u][i]);
					weights.add(g.weights[u][i]);
				}
			}
			
			subGraph.addEdges(u, LowDiameterDecomposition.listToArr(outVertices), LowDiameterDecomposition.listToArr(weights));
		}
		subGraph.initNullAdjListElts();
		
		return subGraph;
	}
	
	public static void findReachable(Graph g, int u, boolean[] reachable) {
		if (!reachable[u]) {
			reachable[u] = true;
			
			for (int v : g.adjacencyList[u]) {
				findReachable(g, v, reachable);
			}
		}
	}
	
	public static void runBellmanFord(Graph g) {
		double time = System.currentTimeMillis();
		g.BellmanFord(0);
		double runTime = System.currentTimeMillis() - time;
		double roundedRunTime = ((int) (runTime * 100)) / 100.0;
		System.out.println("Time to complete Bellman-Ford: " + roundedRunTime + " ms.");
	}
	
	// Runs the bit scaling algorithm of Goldberg and Rao to return a shortest path tree for g_in
	public static int[] bitScaling(Graph g, int s) throws Exception {
		startTime = System.currentTimeMillis();

		int minWeight = Integer.MAX_VALUE;
		for (int u : g.vertices) {
			for (int i = 0; i < g.adjacencyList[u].length; i++) {
				int weight = g.weights[u][i];
				if (weight < minWeight) {
					minWeight = weight;
				}
			}
		}
		
		if (minWeight >= 0) {
			System.out.println("Graph has no negative edge weights.");
			return getShortestPathTree(g, s);
		}
		
		int precision = (int) Math.pow(2, (int) logBase2(-1 * minWeight));
		int[] phi = new int[g.v_max];
		
		while (precision >= 1) {
			Graph gScaledS = new Graph(g.v_max + 1, false); // scaled graph with added dummy source vertex s
			gScaledS.addVertices(g.vertices);
			gScaledS.addVertex(g.v_max);
			
			for (int u : g.vertices) {
				int numOutEdges = g.adjacencyList[u].length;
				int[] edges = new int[numOutEdges];
				int[] weights = new int[numOutEdges];
				
				for (int i = 0; i < numOutEdges; i++) {
					int roundedWeight = phi[u] - phi[g.adjacencyList[u][i]] + (int) Math.ceil(g.weights[u][i] / (double) precision);
					if (CHECKS) {
						if (roundedWeight < -1) {
							throw new Exception("Bit scaling produced an edge of weight less than -1.");
						}
					}
					
					edges[i] = g.adjacencyList[u][i];
					weights[i] = roundedWeight;
				}
				
				gScaledS.addEdges(u, edges, weights);
			}
			
			int[] dummyEdges = new int[g.n];
			int[] dummyWeights = new int[g.n];
			for (int i = 0; i < g.n; i++) {	
				dummyEdges[i] = g.vertices.get(i);
				dummyWeights[i] = 0;
			}
			gScaledS.addEdges(g.v_max, dummyEdges, dummyWeights);
			gScaledS.initNullAdjListElts();
			
			int[] tree = SPmain(gScaledS, g.v_max);
			
			int[] dist = new int[g.v_max + 1]; // will store the shortest distances from the dummy vertex s to every other vertex
			getDistances(gScaledS, tree, dist, 0, g.v_max, new boolean[g.v_max + 1]);
			
			if (CHECKS) {
				verifyTree(gScaledS, tree, dist, g.v_max);
				
				if (hasNegativeEdges(gScaledS, dist, 0)) {
					throw new Exception("Bit scaling failed.");
				}
			}
			
			for (int u: g.vertices) {
				// add then multiply by 2
				if (precision == 1) {
					phi[u] += dist[u];
				} else {
					phi[u] = 2 * (phi[u] + dist[u]);
				}
			}
			
			precision /= 2;
		}
		
		for (int u : g.vertices) {
			for (int i = 0; i < g.adjacencyList[u].length; i++) {
				g.weights[u][i] += phi[u] - phi[g.adjacencyList[u][i]];
				
				if (CHECKS) {
					if (g.weights[u][i] < 0) {
						throw new Exception("After applying the phi outputted from running the bit "
								+ "scaling algorithm and SPMain, there exists an edge with negative weight.");
					}
				}
			}
		}
		
		int[] tree = getShortestPathTree(g, s);
		
		double runTime = System.currentTimeMillis() - startTime;
		double roundedRunTime = ((int) (runTime * 100)) / 100.0;
		System.out.println("The program terminated in " + roundedRunTime + " ms.");
		
		return tree;
	}
	
	public static void verifyTree(Graph g, int[] tree, int[] dist, int src) throws Exception {
		// verify that tree is a valid tree
		
		// construct the graph with all the vertices and only the edges in the shortest path tree
		boolean[][] adjList = new boolean[g.v_max][g.v_max];
		for (int u : g.vertices) {
			if (tree[u] != -1) {
				if (!g.containsVertex[tree[u]]) {
					throw new Exception("Graph structure broken.");
				}
				
				adjList[tree[u]][u] = true;
			}
		}
		
		boolean[] visited = new boolean[g.v_max];
		if (containsCycles(g, adjList, src, visited)) {
			throw new Exception("SPMain returned an invalid tree.");
		}
		
		// verify that there aren't any shorter distances
		for (int u : g.vertices) {
			for (int i = 0; i < g.adjacencyList[u].length; i++) {
				int v = g.adjacencyList[u][i];
				
				if (dist[v] > dist[u] + g.weights[u][i]) {
					throw new Exception("SPMain returned a tree that is not a shortest paths tree.");
				}
			}
		}
	}
	
	public static boolean containsCycles(Graph g, boolean[][] adjList, int src, boolean[] visited) {
		if (visited[src]) {
			return true;
		}
		
		visited[src] = true;
		for (int v = 0; v < g.v_max; v++) {
			if (adjList[src][v]) {
				if (containsCycles(g, adjList, v, visited)) {
					return true;
				}
			}
		}
		
		return false;
	}
	
	public static void getDistances(Graph g, int[] tree, int[] dist, int curDis, int curVertex, boolean[] visited) {
		if (!visited[curVertex]) {
			dist[curVertex] = curDis;
			visited[curVertex] = true;
			
			for (int i = 0; i < g.adjacencyList[curVertex].length; i++) {
				int v = g.adjacencyList[curVertex][i];
				if (tree[v] == curVertex) {
					getDistances(g, tree, dist, curDis + g.weights[curVertex][i], v, visited);
				}
			}
		}
	}
	
	/*
	 * Returns the shortest path tree in g from s.
	 */
	public static int[] SPmain(Graph g_in, int s) throws Exception {	
		int scaleFactor = 2 * g_in.n;
		Graph g = getScaledGraph(g_in, scaleFactor);
		int B = roundPower2(scaleFactor);
		int[] phi = new int[g.v_max];
		
		for (int i = 1; i <= logBase2(B); i++) {
			Graph g_phi = createModifiedGB(g, 0, false, null, phi);
			int[] phi_i = ScaleDown(g_phi, g.n, B / (int) Math.pow(2, i));
			
			if (CHECKS) {
				if (hasNegativeEdges(g_phi, phi_i, B / (int) Math.pow(2, i))) {
					throw new Exception("ScaleDown failed.");
				}
			}
			
			phi = addPhi(phi, phi_i);
		}
		
		// create G^*
		for (int u : g.vertices) {
			for (int i = 0; i < g.adjacencyList[u].length; i++) {
				g.weights[u][i] += phi[u] - phi[g.adjacencyList[u][i]] + 1;
				
				if (CHECKS) {
					if (g.weights[u][i] < 0) {
						throw new Exception("After applying the phi outputted from SPMain, "
								+ "there exists an edge with negative weight.");
					}
				}
			}
		}
		
		int[] tree = getShortestPathTree(g, s);
		
		if (CHECKS) {
			if (invalidTree(g, s, tree)) {
				throw new Exception("SPMain get shortest path tree failed.");
			}
		}
		
		return tree;
	}
	
	public static Graph getScaledGraph(Graph g_in, int scaleFactor) throws Exception {
		Graph g = new Graph(g_in.v_max, false);
		g.addVertices(g_in.vertices);
		
		for (int u : g_in.vertices) {
			int[] edges = new int[g_in.adjacencyList[u].length];
			int[] weights = new int[g_in.adjacencyList[u].length];
			
			for (int i = 0; i < g_in.adjacencyList[u].length; i++) {
				edges[i] = g_in.adjacencyList[u][i];
				weights[i] = g_in.weights[u][i] * scaleFactor;
			}
			g.addEdges(u, edges, weights);
		}
		g.initNullAdjListElts();
		
		return g;
	}
	
	public static boolean invalidTree(Graph g, int s, int[] tree) {
		for (int u = 0; u < tree.length; u++) {
			if (g.containsVertex[u]) {
				if ((u != s) && (tree[u] == -1)) {
					return true;
				}
				if ((u == s) && (tree[u] != -1)) {
					return true;
				}
			}
		}
		return false;
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
	public static int[] ScaleDown(Graph g, int delta, int B) throws Exception {
		// if phi(x) == null, assume that phi(x) = 0
		int[] phi_2 = new int[g.v_max];
		
		if (delta > 2) { // base case
			// double d = delta / 2.0;
			// Graph g_B_nneg = createModifiedGB(g, B, true, null, null);
			
			// phase 0
			// ArrayList<int[]> E_sep = LowDiameterDecomposition.LDD(g_B_nneg, (int) (d * B));
			ArrayList<int[]> E_sep = new ArrayList<int[]>();
			
			HashSet<int[]> E_sep_hash = new HashSet<int[]>(E_sep);
			Graph g_B_Esep = createModifiedGB(g, B, false, E_sep_hash, null);
			ArrayList<ArrayList<Integer>> SCCs = g_B_Esep.SCC();
			
			// phase 1
			int[] vertexToSCCMap = getVertexToSCCMap(SCCs, g.v_max);
			HashSet<int[]> edgesBetweenSCCs = getEdgesBetweenSCCs(g, vertexToSCCMap);
			Graph H = createModifiedGB(g, 0, false, edgesBetweenSCCs, null);
			int[] phi_1 = ScaleDown(H, delta / 2, B);
			
			// phase 2
			Graph g_B_E_sep_phi1 = createModifiedGB(g, B, false, E_sep_hash, phi_1);
			int[] phi = FixDAGEdges(g_B_E_sep_phi1, SCCs, vertexToSCCMap, edgesBetweenSCCs);
			phi_2 = addPhi(phi_1, phi);
			
			if (CHECKS) {
				if (hasNegativeEdges(g_B_Esep, phi_2, 0)) {
					throw new Exception("FixDAGEdges failed.");
				}
			}
		}
		
		// phase 3
		Graph g_B_phi2 = createModifiedGB(g, B, false, null, phi_2);
		int[] phi_prime = ElimNeg(g_B_phi2);
		int[] phi_3 = addPhi(phi_2, phi_prime);
		
		if (CHECKS) {
			if (hasNegativeEdges(g_B_phi2, phi_prime, 0)) {
				throw new Exception("ElimNeg failed.");
			}
		}
		
		return phi_3;
	}
	
	public static boolean hasNegativeEdges(Graph g, int[] phi, int B) {
		for (int u : g.vertices) {
			for (int i = 0; i < g.adjacencyList[u].length; i++) {
				if (g.weights[u][i] + phi[u] - phi[g.adjacencyList[u][i]] < -1 * B) {
					return true;
				}
			}
		}
		
		return false;
	}
	
	
	/*
	 * Creates G^B_phi = (V, E, w^B_phi), where 
	 * w^B_phi(e) = w(e) + phi(u) - phi(v) if w(e) >= 0, 
	 * and w(e) + B + phi(u) - phi(v) if w(e) < 0.
	 * If nneg == true, w(e) = max{0, w^B(e)}.
	 * Removes all the edges in remEdges.
	 */
	public static Graph createModifiedGB(Graph g, int B, boolean nneg, HashSet<int[]> remEdges, int[] phi) throws Exception {
		Graph modG = new Graph(g.v_max, false);
		modG.addVertices(g.vertices);
		
		for (int u : g.vertices) {
			int[] edges = new int[g.adjacencyList[u].length];
			int[] weights = new int[g.adjacencyList[u].length];
			
			for (int i = 0; i < g.adjacencyList[u].length; i++) {
				int v = g.adjacencyList[u][i];
				
				int[] edge = {u, v};
				if ((remEdges == null) || !remEdges.contains(edge)) {
					int weight = g.weights[u][i];
					
					if (weight < 0) {
						weight += B;
					}
					
					if (nneg) {
						weight = Math.max(0, weight);
					}
					
					if (phi != null) {
						weight += phi[u] - phi[v];
					}
					
					edges[i] = v;
					weights[i] = weight;
				}
			}
			
			modG.addEdges(u, edges, weights);
		}
		modG.initNullAdjListElts();
		
		return modG;
	}	
	
	public static HashSet<int[]> getEdgesBetweenSCCs(Graph g, int[] vertexToSCCMap) {
		HashSet<int[]> edgesBetweenSCCs = new HashSet<int[]>();
		for (int u : g.vertices) {
			for (int v : g.adjacencyList[u]) {
				if (vertexToSCCMap[u] != vertexToSCCMap[v]) {
					int[] edge = {u, v};
					edgesBetweenSCCs.add(edge);
				}
			}
		}
		
		return edgesBetweenSCCs;
	}
	
	
	public static int[] getVertexToSCCMap(ArrayList<ArrayList<Integer>> SCCs, int numVertices) {
		int[] vertexToSCCMap = new int[numVertices];
		for (int i = 0; i < SCCs.size(); i++) {
			for (int v : SCCs.get(i)) {
				vertexToSCCMap[v] = i;
			}
		}
		
		return vertexToSCCMap;
	}


	public static int[] addPhi(int[] phi_1, int[] phi_2) throws Exception {
		int len = phi_1.length;
		if (len != phi_2. length) {
			throw new Exception("Trying to add phi's of different lengths.");
		}
		
		int[] newPhi = new int[len];
		
		for (int i = 0; i < len; i++) {
			newPhi[i] = phi_1[i] + phi_2[i];
		}
		
		return newPhi;
	}
	
	
	public static int[] FixDAGEdges(Graph g, ArrayList<ArrayList<Integer>> SCCs, int[] vertexToSCCMap, HashSet<int[]> edgesBetweenSCCs) {
		int n = SCCs.size();
		int[] topOrdering = topSort(n, createSCCAdjList(SCCs, vertexToSCCMap, edgesBetweenSCCs));
		
		int[] mu = new int[n]; // indices are in topological order (e.g., index 0 corresponds to the first SCC in topological order)
		for (int u : g.vertices) {
			for (int i = 0; i < g.adjacencyList[u].length; i++) {
				int v = g.adjacencyList[u][i];
				
				int SCCu = vertexToSCCMap[u];
				int SCCv = vertexToSCCMap[v];
				int edgeWeight = g.weights[u][i];
				
				if ((SCCu != SCCv) && edgeWeight < mu[topOrdering[SCCv]]) {
					mu[topOrdering[SCCv]] = edgeWeight;
				}
			}
		}
		
		int[] phi = new int[g.v_max];
		int m = 0;
		for (int j = 1; j < n; j++) {
			m += mu[j];
			for (int v : SCCs.get(topOrdering[j + n])) {
				phi[v] = m;
			}
		}
		
		return phi;
	}
	
	// returns the adjacency list for the DAG where every SCC is viewed as a single vertex
	
	public static ArrayList<Integer>[] createSCCAdjList(ArrayList<ArrayList<Integer>> SCCs, int[] vertexToSCCMap, HashSet<int[]> edgesBetweenSCCs) {
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] SCCAdjList = new ArrayList[SCCs.size()];
		for (int i = 0; i < SCCs.size(); i++) {
			SCCAdjList[i] = new ArrayList<Integer>();
		}
		boolean[][] containsEdge = new boolean[SCCs.size()][SCCs.size()];
		
		for (int[] edge : edgesBetweenSCCs) {
			int u = vertexToSCCMap[edge[0]];
			int v = vertexToSCCMap[edge[1]];
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
	 * key i + n.
	 */
	public static int[] topSort(int n, ArrayList<Integer>[] adjList) {
		int[] topOrdering = new int[2 * n];
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
        	topOrdering[v] = i;
        	topOrdering[i + n] = v;
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
	public static int[] ElimNeg(Graph g) throws Exception {
		Graph Gs = createGs(g);
		int[] dist = new int[Gs.v_max];
		int s = Gs.v_max - 1;
		
		dist[s] = 0;
		for (int v = 0; v < s; v++) {
			dist[v] = Integer.MAX_VALUE;
		}
		
		PriorityQueue<Node> pq = new PriorityQueue<Node>(Gs.v_max, new Node());
		boolean[] inPQ = new boolean[Gs.v_max];
		pq.add(new Node(s, dist[s]));
		inPQ[s] = true;
		
		HashSet<Integer> marked = new HashSet<Integer>();
		
		while (!pq.isEmpty()) {
			// Dijkstra Phase
			while (!pq.isEmpty()) {
				int v = pq.remove().node;
				inPQ[v] = false;
				
				marked.add(v);
				
				for (int i = 0; i < Gs.adjacencyList[v].length; i++) {
					int x = Gs.adjacencyList[v][i];
					int edgeWeight = Gs.weights[v][i];
					
					if (edgeWeight >= 0 && (dist[v] + edgeWeight < dist[x])) {
						marked.add(x);
						
						if (inPQ[x]) {
							pq.remove(new Node(x, dist[x]));
						}
						dist[x] = dist[v] + edgeWeight;
						pq.add(new Node(x, dist[x]));
						inPQ[x] = true;
					}
				}
			}
			
			// Bellman-Ford Phase
			for (int v : marked) {
				for (int i = 0; i < Gs.adjacencyList[v].length; i++) {
					int x = Gs.adjacencyList[v][i];
					int edgeWeight = Gs.weights[v][i];
					
					if (edgeWeight < 0 && (dist[v] + edgeWeight < dist[x])) {
						if (inPQ[x]) {
							pq.remove(new Node(x, dist[x]));
						}
						dist[x] = dist[v] + edgeWeight;
						pq.add(new Node(x, dist[x]));
						inPQ[x] = true;
					}
				}
			}
			marked = new HashSet<Integer>();
		}
		
		int[] phi = new int[g.v_max];
		for (int v : g.vertices) {
			phi[v] = dist[v];
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
			Gs.addEdges(u, g.adjacencyList[u], g.weights[u]);
		}
		
		int[] edges = new int[g.v_max];
		int[] weights = new int[g.v_max];
		for (int i = 0; i < g.v_max; i++) {
			edges[i] = i;
			weights[i] = 0;
		}
		Gs.addEdges(s, edges, weights);
		
		Gs.initNullAdjListElts();
		return Gs;
	}
	
	
	/*
	 * Let the output be int[] out. 
	 * For vertex i, out[i] = parent of vertex i in the shortest path tree in g from s.
	 * out[i] = -1 indicates that vertex i has no parent in the tree.
	 */
	public static int[] getShortestPathTree(Graph g, int s) throws Exception {	
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
	
	
	public static void updateTreeNeighbors(Graph g, int u, int[] tree, Set<Integer> settled, PriorityQueue<Node> pq, int[] dist) throws Exception {
        for (int i = 0; i < g.adjacencyList[u].length; i++) {
        	int v = g.adjacencyList[u][i];
        	
            if (!settled.contains(v)) {
            	int weight = g.weights[u][i];
                int newDistance = dist[u] + weight;

                if (newDistance < dist[v]) {
                    dist[v] = newDistance;
                    tree[v] = u;
                }
                
                pq.add(new Node(v, dist[v]));
            }
        }
	}
}