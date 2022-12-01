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
	int[][] adjacencyList; // length v_max, each element int[i][j] indicates an edge (i, int[i][j])
	int[][] weights; // length v_max, each element int[i][j] is the weight of edge (i, adjacencyList[i][j])
	
	// Used for calculating SCCs using Tarjan's Algorithm
	int time;
	
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
		
		adjacencyList = new int[numVertices][];
		weights = new int[numVertices][];
		
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
	
	public void addEdges(int v, int[] outVertices, int[] newWeights) throws Exception {
		if ((outVertices.length != newWeights.length) || (adjacencyList[v] != null) || (weights[v] != null)) {
			throw new Exception("Add edges failed.");
		}
		
		adjacencyList[v] = new int[outVertices.length];
		weights[v] = new int[outVertices.length];
		
		for (int i = 0; i < outVertices.length; i++) {
			adjacencyList[v][i] = outVertices[i];
			weights[v][i] = newWeights[i];
		}
	}
	
	public void initNullAdjListElts() {
		for (int i = 0; i < adjacencyList.length; i++) {
			if (adjacencyList[i] == null) {
				adjacencyList[i] = new int[0];
				weights[i] = new int[0];
			}
		}
	}
	
	public void displayGraph() {
		System.out.println("Displaying the graph.");
		for (int i = 0; i < v_max; i++) {
			if (containsVertex[i]) {
				System.out.print("Vertex " + i + ": ");
				
				for (int j = 0; j < adjacencyList[i].length; j++) {
					System.out.print(adjacencyList[i][j] + " ");
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
	 
	// Returns true if the graph has a negative weight cycle.
	// Assumes every vertex is being used.
	public boolean hasNegCycle() {
	    boolean[] visited = new boolean[v_max];
	    int dist[] = new int[v_max];
	 
	    for(int i = 0; i < v_max; i++) {
	        if (visited[i] == false) {
	            if (isNegCycleBellmanFord(i, dist)) {
	                return true;
	            }

	            for(int j = 0; j < n; j++) {
	                if (dist[j] != Integer.MAX_VALUE) {
	                    visited[j] = true;
	                }
	            }
	        }
	    }
	    return false;
	}
	
	// Runs Bellman-Ford to determine whether the graph has a negative cycle.
	public boolean isNegCycleBellmanFord(int src, int dist[]) {
	    for(int i = 0; i < v_max; i++) {
	        dist[i] = Integer.MAX_VALUE;
	    }
	    dist[src] = 0;

	    for(int i = 1; i <= v_max - 1; i++) {
	    	for (int u = 0; u < v_max; u++) {
	    		for (int j = 0; j < adjacencyList[u].length; j++) {
	    			int v = adjacencyList[u][j];
	    			if (dist[u] != Integer.MAX_VALUE &&
    	                dist[u] + weights[u][j] < dist[v]) {
    	                dist[v] = dist[u] + weights[u][j];
    	            }
	    		}
	    	}
	    }
	   
	    for (int u = 0; u < v_max; u++) {
    		for (int j = 0; j < adjacencyList[u].length; j++) {
    			if (dist[u] != Integer.MAX_VALUE &&
		            dist[u] + weights[u][j] < dist[adjacencyList[u][j]]) {
		            return true;
		        }
    		}
    	}
	   
	    return false;
	}
	
	public void BellmanFord(int src) {
        int[] dist = new int[v_max];

        for (int i = 0; i < v_max; i++) {
            dist[i] = Integer.MAX_VALUE;
        }
        dist[src] = 0;
 
        for (int i = 1; i < v_max; i++) {
        	for (int u = 0; u < v_max; u++) {
        		for (int j = 0 ; j < adjacencyList[u].length; j++) {
        			int v = adjacencyList[u][j];
        			int weight = weights[u][j];
        			
        			if (dist[u] != Integer.MAX_VALUE
                            && dist[u] + weight < dist[v]) {
                            dist[v] = dist[u] + weight;
        			}
        		}
        	}
        }
    }
}
