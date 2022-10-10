import java.util.ArrayList;

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
	double[][] weights;
	
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
		weights = new double[numVertices][numVertices];
		
		for (int i = 0; i < numVertices; i++) {
			ArrayList<Integer> initAdjacencyLists = new ArrayList<Integer>();
			adjacencyList[i] = initAdjacencyLists;
		}
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
	
	public void addEdge(int v1, int v2, double w) throws Exception {
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
}
