# Negative Weight SSSP
Code implements the near-linear time algorithm of Aaron Bernstein, Danupon Nanongkai, Christian Wulff-Nilsen that computes the shortest paths from a source vertex in a graph with integral (potentially negative) edge weights (link to paper: https://doi.org/10.48550/arXiv.2203.03456).

### Input file format
The first 10 lines of the input file are parameters for the algorithm._
Line 1: source vertex_
Line 2: 1 to include checks that will throw an exception if something goes wrong with the algorithm, and 0 to ignore the checks._
Line 3: 1 to run low diameter decomposition, 0 to skip low diameter decomposition_
Line 4: 1 to use a random price function to generate negative edge weights, or 2 to use exactly the edge weights given in the input file. A user can also enter a double p between 0 and 1, and with probability p, the algorithm will randomly multiply an edge weight given in the input file by -1._
Line 5: maximum vertex number to read in the input file_
Line 6: 1 to use a random seed, 2 to not use a random seed_
Line 7: If line 6 is 1, give the desired random seed. If line 6 is 2, enter any integer._
Line 8: 1 to display the shortest path tree, 2 to skip displaying the tree_
Line 9: 1 to display the sizes of graphs processed during low diameter decomposition, 2 to avoid displaying the sizes_
Line 10: 1 to make the input graph connected by adding a dummy source vertex and edges of weight 0 between the source vertex and every other vertex in the graph. 2 to skip making the input graph connected._
_
Each of the remaining lines describes an edge in the input graph. Lines must contain the three parameters (vertex 1, vertex 2, integral edge weight), and the entries can be space or comma-separated._

Examples:_
1 3 10_
0,2,14_
2,2,-1_
_
For an example input file, see input_example.txt.