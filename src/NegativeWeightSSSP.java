
public class NegativeWeightSSSP {

	public static void main(String[] args) throws Exception {
	    Graph g1 = new Graph(5, true);
	 
	    g1.addEdge(1, 0, 1);
	    g1.addEdge(0, 2, 1);
	    g1.addEdge(2, 1, 1);
	    g1.addEdge(0, 3, 1);
	    g1.addEdge(3, 4, 1);
	    System.out.println("SSC in first graph ");
	    g1.SCC();
	 
	    Graph g2 = new Graph(4, true);
	    g2.addEdge(0, 1, 1);
	    g2.addEdge(1, 2, 1);
	    g2.addEdge(2, 3, 1);
	    System.out.println("\nSSC in second graph ");
	    g2.SCC();
	}

}
