//Import statements
import java.util.*;

/**
 * Outputs the the number of trees in the given list of edges and the total weight of all
 * the minimum-weight spanning trees in the forest, using Prim's algorithm with a heap as 
 * a priority queue.
 * 
 * @author Jeremy Maness
 * @author Gabriel Oliveira
 *
 */
public class prim {

	// Global variables 
	//Adjacency matrix
	static int[][] adjMatrix = null;
	
	//Number of vertices
	static int n = -1;
	
	//Number of edges
	static int m = -1;
	
	//branching factor
	static int brcFactor = -1;
	
	//Scanner to read the input
	private static Scanner scanner;
	
	
	/**
	 * Reads input file and fill up the adjacency matrix.
	 * Calculates branching factor
	 * 
	 * @param args arguments on user input
	 */
	public static void main(String[] args) {
		
		// Reading input
		scanner = new Scanner(System.in);

        // The first number in the stream is the number of vertices
        if (scanner.hasNextInt()) {
        	n = scanner.nextInt();
        }
        System.out.println(n);
        // The second number in the stream is the number of edges
        if (scanner.hasNextInt()) {
        	m = scanner.nextInt();
        }
		System.out.println(m);
        adjMatrix = new int[n][n];
        
        // Calculating branching factor
        // based number of vertices and edges from input
        brcFactor = calculateBranchingFactor(n,m);
        

        while (scanner.hasNext()) {

        	int a = scanner.nextInt();
        	int b = scanner.nextInt();
        	int w = scanner.nextInt();

        	adjMatrix[a][b] = w;
        }

	}
	
	/**
	 * Calculates the branch factor given the number of vertices and the number of edges
	 * @param vertices Number of vertices
	 * @param edges Number of edges
	 * @return branching factor
	 */
	public static int calculateBranchingFactor(int vertices, int edges) {
		return (int) Math.pow(2, Math.ceil(Math.log((float) 6/9) / Math.log(2)));
	}

}
