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
	
	//Array to locate vertices in the heap
	static int[] vertexLocations = null;
	
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
        
        // The second number in the stream is the number of edges
        if (scanner.hasNextInt()) {
        	m = scanner.nextInt();
        }
		
        //Initializing the adjacency matrix
        adjMatrix = new int[n][n];
        
        //Initializing array of vertices locations
        vertexLocations = new int[n];
        
        // Calculating branching factor
        // based number of vertices and edges from input
        brcFactor = calculateBranchingFactor(n,m);
        

        // Reading edges given on the input
        while (scanner.hasNext()) {
        	
        	//First vertex of edge
        	int a = scanner.nextInt();
        	//Second vertex of edge
        	int b = scanner.nextInt();
        	//Weight of the edge
        	int w = scanner.nextInt();

        	//Filling adjacency matrix 
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
	
	
	
	
	
	
	///////////////////////////////////////// Vertex inner class from homework 2 //////////////////////////////////////////
	public class Vertex {
	    private int id = -1;
	    //Log base 2 of the branching factor (used to speed up the calculations of the indexes of the parent/children)
	    private int weight = -1;
	    //counter to track the number of key comparisons
	    private Vertex parent = null;
	    //Flag	    
	    private boolean partOfSpanningTree = false;
	}
	
	
	
	
	

	///////////////////////////////////////// Tree inner class from homework 2 //////////////////////////////////////////
	public class Tree {
		
	}
	
	
	
	
	
	
	
	///////////////////////////////////////// Heap inner class from homework 2 //////////////////////////////////////////
	
	
	/**
	* This class constructs and tests a min-heap in which the nodes can have an arbitrary number of children, as long
	*  as the number is a power of 2.
	* The branching factor (number of children per node) is entered by the user as a command line parameter.
	* The testing works by reading a file with values and keys, which will be entered by the user in the command line
	*  when running the program, and the file will be read using file redirection.
	*/
	public class heap {
	    //array that will keep the values of the heap
	    private ArrayList<Node> array = new ArrayList<>();
	    //size of the heap
	    private int heapSize = 0;
	    //Log base 2 of the branching factor (used to speed up the calculations of the indexes of the parent/children)
	    private int lgBranchingFactor;
	    //counter to track the number of key comparisons
	    private int keyComparisons = 0;

	    /**
	    * Heap constructor.
	    * Initiates the branchingFactor field and the lgBranchingFactor field.
	    * branchingFactor is entered by the user and passed as a parameter from the
	    * main method.
	    */
	    private heap(int brcFactor) {
	        this.lgBranchingFactor = (int) brcFactor;
	    }

	    /**
	    * Inner class for creation of the node.
	    * A Node object represents a node of the heap, and it is composed
	    * by the key and the value.
	    */
	    private class Node {
	        //Key to be compared
	        private final int key;
	        //Value stored on the node
	        private final int value;

	        /**
	        * Constructor
	        * Initiates the key and the value with information from the file
	        * entered by the user.
	        */
	        Node(int key, int value) {
	            this.key = key;
	            this.value = value;
	        }

	        @Override
	        /**
	        * Override the toString of the Node class to format the string output
	        * of a Node
	        */
	        public String toString() {
	            return String.format("%s %s", key, value);
	        }
	    }

	    /**
	     * Inserts the specified node into the heap at the correct position to maintain the min-heap property
	     * 
	     * @param node Node to insert
	     */
	    private void insertValue(Node node) {
	        heapSize++;
	        //Add node to the end of the array
	        array.add(node);
	        
	        //Initiating from the bottom of the heap and going up, makes key comparisons to find the
	        //right place to insert the given node.
	        int i = heapSize - 1;
	        while (i > 0 && keysGreaterThan(array.get(parent(i)), node)) {
	            array.set(i, array.get(parent(i)));
	            i = parent(i);
	        }
	        //Sets the node into the right index after finding the right position.
	        array.set(i, node);
	    }

	    /**
	     * Removes and returns the node with the minimum key in the heap.
	     * 
	     * @return node with the minimum key in the heap
	     */
	    private Node removeMin() {
	        //Throws exception if the heap is empty.
	        if (heapSize < 1) {
	            throw new IllegalStateException("heap underflow");
	        }
	        
	        //Stores the node with the smallest key in the heap (root)
	        Node min = array.get(0);
	        //Replace the root with the last node and decrease the heap size
	        array.set(0, array.get(heapSize - 1));
	        heapSize--;
	        //Calls the Heapify operation to reposition the new root into the 
	        //right position
	        minHeapify(0);
	        //Returns the previous smallest node that was stored
	        return min;
	    }

	    /**
	     * Compares keys of the specified nodes and returns true if the key of node1
	     * is greater than the key of node2.
	     *
	     * This also has a side effect of incrementing the number of key comparisons
	     * during the lifetime of this heap.
	     *
	     * @param node1 Node
	     * @param node2 Node
	     * @return true if the key of node1 is greater than the key of node2
	     */
	    private boolean keysGreaterThan(Node node1, Node node2) {
	        keyComparisons++;
	        return node1.key > node2.key;
	    }

	    /**
	     * Compares keys of the specified nodes and returns true if the key of node1
	     * is less than the key of node2.
	     *
	     * This also has a side effect of incrementing the number of key comparisons
	     * during the lifetime of this heap.
	     *
	     * @param node1 Node
	     * @param node2 Node
	     * @return true if the key of node1 is less than the key of node2
	     */
	    private boolean keysLessThan(Node node1, Node node2) {
	        keyComparisons++;
	        return node1.key < node2.key;
	    }

	    /**
	     * Heapifies this heap starting at the specified index.
	     *
	     * @param i index to start heap-ification.
	     */
	    private void minHeapify(int i) {
	        int smallest = i;

	        for (int j = 1; j <= brcFactor; j++) {
	            int childIndex = child(i, j);

	            if (childIndex < heapSize && keysLessThan(array.get(childIndex), array.get(smallest))) {
	                smallest = childIndex;
	            }
	        }

	        if (smallest != i) {
	            exchange(i, smallest);
	            minHeapify(smallest);
	        }
	    }

	    /**
	     * Returns the index of the parent for the child at the ith index in the heap.
	     * 
	     * @param i index of child
	     * @return the index of the parent
	     */
	    private int parent(int i) {
	        return (i - 1) >> lgBranchingFactor;
	    }

	    /**
	     * Returns the index of the jth child of the ith node
	     *
	     * @param i index of current node
	     * @param j jth child of the ith node
	     * @return the index of the jth child of the ith node
	     */
	    private int child(int i, int j) {
	        return (i << lgBranchingFactor) + j;
	    }

	    /**
	     * Swaps the elements at the specified src and dest indexes of the array
	     *
	     * @param src source index
	     * @param dest destination index
	     */
	    private void exchange(int src, int dest) {
	        Node tmp = array.get(dest);
	        array.set(dest, array.get(src));
	        array.set(src, tmp);
	    }

	    /**
	    * Gets the branching factor entered by the user and create an array out of the values
	    * from the file also entered by the user (read by file redirection).
	    *
	    * Prints error messages if the user entered more than one value for branching factor and
	    * if the entered branching factor is not a power of 2.
	    *
	    * @param args array containing branching factor entered in the command line
	    */
	    public void createHeap(String[] args) {

	    	if (args.length > 1) {
	            System.err.println("usage: heap <branching factor>");
	            System.exit(1);
	        }

	        int branchingFactor = args.length == 0 ? 2 : Integer.parseInt(args[0]);

	        if (!isPow2(branchingFactor)) {
	            System.err.println(String.format("Branching factor %s is not a factor of 2", branchingFactor));
	            System.exit(2);
	        }
	        
	        heap A = new heap(brcFactor);

	        Scanner sc = new Scanner(System.in);
	        while (sc.hasNextLine()) {
	            String line = sc.nextLine();

	            if (line.trim().equals("-1")) {
	                Node min = A.removeMin();
	                System.out.println(min);
	            } else {
	                // input has two integers
	                int[] parts = Arrays.stream(line.split("\\s+"))
	                        .mapToInt(Integer::parseInt)
	                        .toArray();

	                A.insertValue(new Node(parts[0], parts[1]));
	            }
	        }

	        System.out.println("key comparisons: " + A.keyComparisons);
	    }

	    /**
	     * Returns true if the specified integer is a power of 2
	     * (i.e. there exists an integer, k, such that 2^k = val)
	     *
	     * @param val integer value
	     * @return true if the specified integer is a power of 2
	     */
	    private boolean isPow2(int val) {
	        return Integer.bitCount(val) == 1;
	    }
	}

}
