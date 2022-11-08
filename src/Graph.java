import java.util.ArrayList;
import java.util.Stack;

public class Graph {
	/*
	 *  Upper bound on the greatest index of a vertex. 
	 *  Important for when creating subgraphs, since only a subset 
	 *  of the vertices from 1 to v_max will exist in the subgraph.
	 *  E.g., if v_max = 5, the set of vertices could be {0, 3, 4}.
	 */
	int v_max;
	ArrayList<Integer> vertices; // the specific vertices that are in G
	boolean[] containsVertex;
	
	int n; // number of vertices
	ArrayList<Integer>[] adjacencyList; // for 1 <= i <= v_max, contains an integer Arraylist where an element j indicates an edge (i, j)
	boolean[][] containsEdge;
	int[][] weights;
	
	// Used for calculating SCCs using Tarjan's Algorithm
	int time;
	
	@SuppressWarnings("unchecked")
	public Graph(int numVertices, boolean withAllVertices) {
		v_max = numVertices;
		vertices = new ArrayList<Integer>();
		containsVertex = new boolean[numVertices];
		
		if (withAllVertices) {
			for (int v = 0; v < numVertices; v++) {
				vertices.add(v);
				containsVertex[v] = true;
			}
			n = numVertices;
		} else {
			n = 0;
		}
		
		adjacencyList = new ArrayList[numVertices];
		containsEdge = new boolean[numVertices][numVertices];
		weights = new int[numVertices][numVertices];
		
		for (int i = 0; i < numVertices; i++) {
			ArrayList<Integer> initAdjacencyLists = new ArrayList<Integer>();
			adjacencyList[i] = initAdjacencyLists;
		}
		
		time = 0;
	}
	
	public void addVertices(ArrayList<Integer> verticesToAdd) throws Exception {
		for (int v : verticesToAdd) {
			addVertex(v);
		}
	}
	
	public void addVertex(int v) throws Exception {
		if (0 <= v && v < v_max && !containsVertex[v]) {
			vertices.add(v);
			containsVertex[v] = true;
			n++;
		} else {
			throw new Exception("Add vertex failed.");
		}
	}
	
	public void removeVertex(int v) throws Exception {
		if (0 <= v && v < v_max && containsVertex[v]) {
			vertices.remove(Integer.valueOf(v));
			containsVertex[v] = false;
			n--;
		} else {
			throw new Exception("Remove vertex failed.");
		}
	}
	
	public void addEdge(int v1, int v2, int w) throws Exception {
		if (0 <= v1 && v1 < v_max && 0 <= v2 && v2 < v_max 
				&& containsVertex[v1] && containsVertex[v2] 
						&& !containsEdge[v1][v2]) {
			adjacencyList[v1].add(v2);
			containsEdge[v1][v2] = true;
			weights[v1][v2] = w;
		} else {
			throw new Exception("Add edge failed.");
		}
	}
	
	public void removeEdge(int v1, int v2) throws Exception {
		if (0 <= v1 && v1 < v_max && 0 <= v2 && v2 < v_max 
				&& containsVertex[v1] && containsVertex[v2] 
						&& containsEdge[v1][v2]) {
			adjacencyList[v1].remove(Integer.valueOf(v2));
			containsEdge[v1][v2] = false;
		} else {
			throw new Exception("Remove edge failed.");
		}
	}
	
	public void displayGraph() {
		System.out.println("Displaying the graph.");
		for (int i = 0; i < v_max; i++) {
			if (containsVertex[i]) {
				System.out.print("Vertex " + i + ": ");
				
				for (int j = 0; j < adjacencyList[i].size(); j++) {
					System.out.print(adjacencyList[i].get(j) + " ");
				}
				System.out.println();
			}
		}
	}
	
	// SCC only run on graphs where n == v_max
	// Runs Tarjan's Algorithm for finding SCCs
	ArrayList<ArrayList<Integer>> SCC() {
	    int[] disc = new int[n];
	    int[] low = new int[n];
	    time = 0;
	    for(int i = 0; i < n; i++) {
	        disc[i] = -1;
	        low[i] = -1;
	    }
	     
	    boolean[] stackMember = new boolean[n];
	    Stack<Integer> st = new Stack<Integer>();
	     
	    ArrayList<ArrayList<Integer>> SCCverts = new ArrayList<ArrayList<Integer>>();
	    for(int i = 0; i < n; i++) {
	        if (disc[i] == -1) {
	        	SCCverts.addAll(SCCUtil(i, low, disc, stackMember, st));
	        }
	    }
	    
	    return SCCverts;
	}
		
	public ArrayList<ArrayList<Integer>> SCCUtil(int u, int low[], int disc[], boolean stackMember[], Stack<Integer> st) {
		ArrayList<ArrayList<Integer>> SCCverts = new ArrayList<ArrayList<Integer>>(); 
		
		disc[u] = time;
	    low[u] = time;
	    time += 1;
	    stackMember[u] = true;
	    st.push(u);
	     
	    for (int v : adjacencyList[u]) {
	        if (disc[v] == -1) {
	            SCCverts.addAll(SCCUtil(v, low, disc, stackMember, st));
	            low[u] = Math.min(low[u], low[v]);
	        } else if (stackMember[v] == true) {
	            low[u] = Math.min(low[u], disc[v]);
	        }
	    }
	    
	    int w = -1;
	    if (low[u] == disc[u]) {
	    	ArrayList<Integer> oneSCC = new ArrayList<Integer>();
	        while (w != u) {
	            w = (int) st.pop();
	            oneSCC.add(w);
	            stackMember[w] = false;
	        }
	        SCCverts.add(oneSCC);
	    }
	    
	    return SCCverts;
	}
}
