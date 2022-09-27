import java.util.ArrayList;

public class Graph {
	int n; // number of vertices
	ArrayList<Integer>[] adjacencyList; // for 1 <= i <= n, contains an integer Arraylist where an element j indicates an edge (i, j)
	boolean[][] containsEdge;
	double[][] weights;
	
	@SuppressWarnings("unchecked")
	public Graph(int numVertices) {
		n = numVertices;
		adjacencyList = new ArrayList[n];
		containsEdge = new boolean[n][n];
		weights = new double[n][n];
		
		for (int i = 0; i < n; i++) {
			ArrayList<Integer> initAdjacencyLists = new ArrayList<Integer>();
			adjacencyList[i] = initAdjacencyLists;
		}
	}
	
	public void addEdge(int v1, int v2, double w) {
		if (0 <= v1 && v1 <= n && 0 <= v2 && v2 <= n && !containsEdge[v1][v2]) {
			adjacencyList[v1].add(v2);
			containsEdge[v1][v2] = true;
		}
		weights[v1][v2] = w;
	}
	
	public void removeEdge(int v1, int v2) {
		if (0 <= v1 && v1 <= n && 0 <= v2 && v2 <= n && containsEdge[v1][v2]) {
			adjacencyList[v1].remove(Integer.valueOf(v2));
			containsEdge[v1][v2] = false;
		}
	}
	
	public void displayGraph() {
		System.out.println("Displaying the graph.");
		for (int i = 0; i < n; i++) {
			System.out.print("Vertex " + i + ": ");
			
			for (int j = 0; j < adjacencyList[i].size(); j++) {
				System.out.print(adjacencyList[i].get(j) + " ");
			}
			System.out.println();
		}
	}
	
	public boolean hasEdge(int v1, int v2) {
		return containsEdge[v1][v2];
	}
	
	public double getWeight(int v1, int v2) {
		return weights[v1][v2];
	}
	
	public ArrayList<Integer>[] getAdjacencyList() {
		return adjacencyList;
	}
	
	public ArrayList<Integer> getEdges(int v) {
		return adjacencyList[v];
	}
}
