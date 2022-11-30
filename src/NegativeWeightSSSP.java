import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.Stack;

public class NegativeWeightSSSP {
	
	public static double startTime; // in ms
	public static final double maxRunTime = 5000; // in seconds

	public static void main(String[] args) throws Exception {
		Graph g = readInput("USA-very-small");
		System.out.println(g.hasNegCycle());
		
//		int[] tree = bitScaling(g, 1, false);
//
//		for (int i = 0; i < g_size; i++) {
//			System.out.println("Parent of vertex " + i + ": " + tree[i]);
//		}
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
		
		String line = f.readLine();
		while (line != null) {
			String[] arr = line.split(" ");
			
			edges[Integer.parseInt(arr[1])].add(Integer.parseInt(arr[2]));
//			if (Math.random() < 0) {
//				weights[Integer.parseInt(arr[1])].add(-1);
//			} else {
//				weights[Integer.parseInt(arr[1])].add(Integer.parseInt(arr[3]));
//			}
			weights[Integer.parseInt(arr[1])].add(Integer.parseInt(arr[3]));
			line = f.readLine();
		}
		f.close();
		
		for (int i = 0; i < g_size; i++) {
			g.addEdges(i, LowDiameterDecomposition.listToArr(edges[i]), LowDiameterDecomposition.listToArr(weights[i]));
		}
		g.initNullAdjListElts();
		
		return g;
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
			for (int i = 0; i < g.adjacencyList[u].length; i++) {
				int weight = g.weights[u][i];
				if (weight < minWeight) {
					minWeight = weight;
				}
			}
		}
		
		if (minWeight >= 0) {
			return getShortestPathTree(g_in, s);
		}
		
		int precision = (int) Math.pow(2, (int) logBase2(-1 * minWeight));
		int[] phi = new int[g.v_max];
		
		while (precision >= 1) {
			Graph gScaledS = new Graph(g.n + 1, true); // scaled graph with added dummy source vertex s
			
			for (int u : g.vertices) {
				int numOutEdges = g.adjacencyList[u].length;
				int[] edges = new int[numOutEdges];
				int[] weights = new int[numOutEdges];
				
				for (int i = 0; i < numOutEdges; i++) {
					int roundedWeight = phi[u] - phi[g.adjacencyList[u][i]] + (int) Math.ceil(g.weights[u][i] / precision);
					if (roundedWeight < -1) {
						throw new Exception("Bit scaling produced an edge of weight less than -1.");
					}
					
					edges[i] = g.adjacencyList[u][i];
					weights[i] = roundedWeight;
				}
				
				gScaledS.addEdges(u, edges, weights);
			}
			
			int[] dummyEdges = new int[g.n];
			int[] dummyWeights = new int[g.n];
			for (int v = 0; v < g.n; v++) {	
				dummyEdges[v] = v;
				dummyWeights[v] = 0;
			}
			gScaledS.addEdges(g.n, dummyEdges, dummyWeights);
			gScaledS.initNullAdjListElts();
			
			int[] tree = SPmain(gScaledS, g.n);
			
			int[] dist = new int[g.n + 1]; // will store the shortest distances from the dummy vertex s to every other vertex
			getDistances(gScaledS, tree, dist, 0, g.n, new boolean[g.n + 1]);
			
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
	 * The boolean constantOutDegree is false if we wish to ignore the assumption that g has constant out-degree.
	 */
	public static int[] SPmain(Graph g_in, int s) throws Exception {		
		Graph g = new Graph(g_in.n, true);
		
		for (int u : g_in.vertices) {
			int[] edges = new int[g_in.adjacencyList[u].length];
			int[] weights = new int[g_in.adjacencyList[u].length];
			
			for (int i = 0; i < g_in.adjacencyList[u].length; i++) {
				edges[i] = g_in.adjacencyList[u][i];
				weights[i] = g_in.weights[u][i] * 2 * g_in.n;
			}
			g.addEdges(u, edges, weights);
		}
		g.initNullAdjListElts();
		
		int B = roundPower2(2 * g.n);
		int[] phi = new int[g.v_max];
		
		for (int i = 1; i <= logBase2(B); i++) {
			Graph g_phi = createModifiedGB(g, 0, false, null, phi);
			int[] phi_i = ScaleDown(g_phi, g.n, B / (int) Math.pow(2, i));
			phi = addPhi(phi, phi_i);
		}
		
		// create G^*
		for (int u : g.vertices) {
			for (int i = 0; i < g.adjacencyList[u].length; i++) {
				g.weights[u][i] += phi[u] - phi[g.adjacencyList[u][i]] + 1;
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
			int outDegree = g.adjacencyList[v].length;
			if (outDegree == 0) {
				numVertices += 2;
			} else {
				numVertices += outDegree;
			}
		}
		
		Graph gOut = new Graph(numVertices, true);
		
		for (int u = 0; u < g.n; u++) {
			int outDegree = g.adjacencyList[u].length;
			
			if (outDegree == 0) {
				int[] edges0 = new int[2];
				int[] weights0 = new int[2];
				int[] edges1 = new int[2];
				int[] weights1 = new int[2];
				
				// self-loops
				edges0[0] = indexOfMainVertex[u];
				weights0[0] = 0;
				edges1[0] = indexOfMainVertex[u] + 1;
				weights1[0] = 0;
				
				// loop between v and v'
				edges0[1] = indexOfMainVertex[u] + 1;
				weights0[0] = 0;
				edges1[1] = indexOfMainVertex[u];
				weights1[1] = 0;
				
				gOut.addEdges(indexOfMainVertex[u], edges0, weights0);
				gOut.addEdges(indexOfMainVertex[u] + 1, edges1, weights1);
			} else if (outDegree == 1) {
				int[] edges = new int[2];
				int[] weights = new int[2];
				
				// original edge
				edges[0] = g.adjacencyList[u][0];
				weights[0] = 0;
				// create a self-loop
				edges[1] = indexOfMainVertex[u];
				weights[1] = 0;
				
				gOut.addEdges(indexOfMainVertex[u], edges, weights);
			} else {
				int i = 0;
				for (int v : g.adjacencyList[u]) {
					int[] edges = new int[2];
					int[] weights = new int[2];
					
					// edge (u, v) => (u + i, indexOfMainVertex[v]), i in [0, g.adjacencyList[u].length)
					edges[0] = indexOfMainVertex[v];
					weights[0] = 0;
					
					// create a cycle
					if (i == outDegree - 1) {
						edges[1] = indexOfMainVertex[u];
					} else {
						edges[1] = indexOfMainVertex[u] + i + 1;
					}
					weights[1] = 0;
					
					gOut.addEdges(indexOfMainVertex[u] + i, edges, weights);
					
					i++;
				}
			}
		}
		gOut.initNullAdjListElts();
		
		return gOut;
	}
	
	// converts the tree for gOut into a tree for g_in
	public static int[] treeForInputG(Graph g, int s, int[] tree) {
		int[] newVerttoOldVert = new int[g.v_max];
		int numVertices = 0;
		HashSet<Integer> originalVertices = new HashSet<Integer>();
		for (int v = 0; v < g.n; v++) {
			originalVertices.add(numVertices);
			int outDegree = g.adjacencyList[v].length;
			
			if (outDegree == 0) {
				newVerttoOldVert[numVertices] = v;
				newVerttoOldVert[numVertices + 1] = v;
				numVertices += 2;
			} else {
				for (int i = 0; i < outDegree; i++) {
					newVerttoOldVert[numVertices + i] = v;
				}
				numVertices += outDegree;
			}
		}
		
		int[] treeInputG = new int[g.n];
		for (int v = 0; v < numVertices; v++) {
			if (originalVertices.contains(v)) {
				if (tree[v] == -1) {
					treeInputG[newVerttoOldVert[v]] = -1;
				} else {
					treeInputG[newVerttoOldVert[v]] = newVerttoOldVert[tree[v]];
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
	public static int[] ScaleDown(Graph g, int delta, int B) throws Exception {
		if (System.currentTimeMillis() - startTime >= maxRunTime * Math.pow(10, 3)) {
			throw new Exception("Exceeded total allotted run time.");
		}
		
		// if phi(x) == null, assume that phi(x) = 0
		int[] phi_2 = new int[g.v_max];
		
		if (delta > 2) { // base case
			double d = delta / 2.0;
			Graph g_B_nneg = createModifiedGB(g, B, true, null, null);
			
			// phase 0
			ArrayList<int[]> E_sep = LowDiameterDecomposition.LDD(g_B_nneg, (int) (d * B));
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
		}
		
		// phase 3
		Graph g_B_phi2 = createModifiedGB(g, B, false, null, phi_2);
		int[] phi_prime = ElimNeg(g_B_phi2);
		int[] phi_3 = addPhi(phi_2, phi_prime);
		
		return phi_3;
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
		pq.add(new Node(s, dist[s]));
		HashSet<Integer> marked = new HashSet<Integer>();
		
		while (!pq.isEmpty()) {
			// Dijkstra Phase
			while (!pq.isEmpty()) {
				int v = pq.remove().node;
				marked.add(v);
				
				for (int i = 0; i < Gs.adjacencyList[v].length; i++) {
					int x = Gs.adjacencyList[v][i];
					int edgeWeight = Gs.weights[v][i];
					
					if (edgeWeight >= 0 && (dist[v] + edgeWeight < dist[x])) {
						marked.add(v);
						pq.add(new Node(x, dist[x]));
						dist[x] = dist[v] + edgeWeight;
					}
				}
			}
			
			// Bellman-Ford Phase
			for (int v : marked) {
				for (int i = 0; i < Gs.adjacencyList[v].length; i++) {
					int x = Gs.adjacencyList[v][i];
					int edgeWeight = Gs.weights[v][i];
					
					if (edgeWeight < 0 && (dist[v] + edgeWeight < dist[x])) {
						dist[x] = dist[v] + edgeWeight;
						pq.add(new Node(x, dist[x]));
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
        for (int i = 0; i < g.adjacencyList[u].length; i++) {
        	int v = g.adjacencyList[u][i];
        	
            if (!settled.contains(v)) {
                int newDistance = dist[u] + g.weights[u][i];

                if (newDistance < dist[v]) {
                    dist[v] = newDistance;
                    tree[v] = u;
                }
                
                pq.add(new Node(v, dist[v]));
            }
        }
	}
}
