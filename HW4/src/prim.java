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
    static Graph g = null;
    
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
        
        // The second number in the stream is the number of edges
        if (scanner.hasNextInt()) {
            m = scanner.nextInt();
        }
        
        //Initializing the adjacency matrix
        //adjMatrix = new int[n][n];
        
        // Calculating branching factor
        // based number of vertices and edges from input
        brcFactor = calculateBranchingFactor(n,m);
        
        ArrayList<Set<Neighbor>> adjList = new ArrayList<Set<Neighbor>>();
        Set<Vertex> vertices = new HashSet<Vertex>();
        
        // Reading edges given on the input
        while (scanner.hasNext()) {
            //First vertex of edge
            int a = scanner.nextInt();
            //Second vertex of edge
            int b = scanner.nextInt();
            //Weight of the edge
            int w = scanner.nextInt();

            
            vertices.add(new Vertex(a, null, null, false));
            vertices.add(new Vertex(b, null, null, false));
            adjList.get(a).add(new Neighbor(b, w));
        }
        
        g = new Graph(adjList, vertices);

    }
    
    /**
     * Calculates the branch factor given the number of vertices and the number of edges
     * @param vertices Number of vertices
     * @param edges Number of edges
     * @return branching factor
     */
    public static int calculateBranchingFactor(int vertices, int edges) {
        return (int) Math.pow(2, Math.ceil(Math.log((float) edges/vertices) / Math.log(2)));
    }
    
    
    
    
    
    ////////////////////////////////////////////////// Prim's algorithm ///////////////////////////////////////////////////
    public void prim(Graph g, int weight, Vertex root) {
        
        heap priorityQueue = new heap(brcFactor);
        
        root.setDistance(0);
        
        for(Vertex v : g.getVertices()) {
            priorityQueue.insertValue(new prim.heap.Node(v.getDistance(), v.getID()));
        }
        
        while (priorityQueue.heapSize != 0) {
            prim.heap.Node u = priorityQueue.removeMin();
                for (Neighbor v : g.getAdjVertices(u.value)) {
                    //if () {
                        
                    //}
            }
        }
        
        
        
    }
    
    //////////////////////////////////////////////////Graph ///////////////////////////////////////////////////
    static class Graph {
        
        private ArrayList<Set<Neighbor>> adjList = new ArrayList<Set<Neighbor>>();
        private Set<Vertex> vertices = new HashSet<Vertex>();
        
        private Graph(ArrayList<Set<Neighbor>> adjList, Set<Vertex> vertices) {
            setAdjList(adjList);
            setVertices(vertices);
        }
        
        public void setAdjList(ArrayList<Set<Neighbor>> adjList) {
            this.adjList = adjList;
        }
        
        public void setVertices(Set<Vertex> vertices) {
            this.vertices = vertices;
        }
        
        public Set<Neighbor> getAdjVertices(int u) {
            return adjList.get(u);
        }
        
        public Set<Vertex> getVertices(){
            return this.vertices;
        }
        
        public int getEdgeWeight(int u, int v) {
            for (Neighbor n : adjList.get(u)) {
                if (n.getVertexID() == v) {
                    return n.getEdgeWeight();
                }
            }
            return -1;
        }
        
        
    }
    
    ///////////////////////////////////////// Vertex inner class from homework 2 //////////////////////////////////////////
    static class Vertex {
        
        private Integer id = null;
        //Log base 2 of the branching factor (used to speed up the calculations of the indexes of the parent/children)
        private Integer distance = null;
        //Parent of this vertex
        private Vertex parent = null;
        //Flag        
        private boolean partOfSpanningTree = false;
        
        
        private Vertex (Integer id, Integer distance, Vertex parent, boolean partOfSpanningTree) {
            setID(id);
            setDistance(distance);
            setParent(parent);
            setMSPFlag(partOfSpanningTree);
        }
        
        public void setID(int id) {
            this.id = id;
        }
        
        public void setDistance(int distance) {
            this.distance = distance;
        }
        
        public void setParent(Vertex parent) {
            this.parent = parent;
        }
        
        public void setMSPFlag(boolean partOfSpanningTree) {
            this.partOfSpanningTree = partOfSpanningTree;
        }
        
        public int getID() {
            return this.id;
        }
        
        public int getDistance() {
            return this.distance;
        }
        
        public Vertex getParent() {
            return this.parent;
        }
        
        public boolean getMSPFlag() {
            return this.partOfSpanningTree;
        }

        /* (non-Javadoc)
         * @see java.lang.Object#hashCode()
         */
        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + id;
            return result;
        }

        /* (non-Javadoc)
         * @see java.lang.Object#equals(java.lang.Object)
         */
        @Override
        public boolean equals(Object obj) {
            if (this == obj)
                return true;
            if (obj == null)
                return false;
            if (!(obj instanceof Vertex))
                return false;
            Vertex other = (Vertex) obj;
            if (id != other.id)
                return false;
            return true;
        }
        
        
    }
    

    ////////////////////////////////////////////////// Neighbor (Adjacency list nodes) ///////////////////////////////////////////////////
    static class Neighbor {
        int vertexID = -1;
        int edgeWeight = -1;
        
        private Neighbor(int vertexID, int edgeWeight) {
            setVertexID(vertexID);
            setEdgeWeight(edgeWeight);
        }
        
        public void setVertexID(int vertexID) {
            this.vertexID = vertexID;
        }
        
        public void setEdgeWeight(int edgeWeight) {
            this.edgeWeight = edgeWeight;
        }
        
        public int getVertexID() {
            return this.vertexID;
        }
        
        public int getEdgeWeight() {
            return this.edgeWeight;
        }
        
        /* (non-Javadoc)
         * @see java.lang.Object#hashCode()
         */
        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + edgeWeight;
            result = prime * result + vertexID;
            return result;
        }

        /* (non-Javadoc)
         * @see java.lang.Object#equals(java.lang.Object)
         */
        @Override
        public boolean equals(Object obj) {
            if (this == obj)
                return true;
            if (obj == null)
                return false;
            if (!(obj instanceof Neighbor))
                return false;
            Neighbor other = (Neighbor) obj;
            if (edgeWeight != other.edgeWeight)
                return false;
            if (vertexID != other.vertexID)
                return false;
            return true;
        }
        
        
    
    }
    
    

    ///////////////////////////////////////// Tree inner class from homework 2 //////////////////////////////////////////
    public class Tree {
        
    }


    
    
    ///////////////////////////////////////// Heap inner class from homework 2 //////////////////////////////////////////
    
    
    	
	///////////////////////////////////////// Heap inner class from homework 2 //////////////////////////////////////////
	
	
    /**
    * This class constructs and tests a min-heap in which the nodes can have an arbitrary number of children, as long
    *  as the number is a power of 2.
    * The branching factor (number of children per node) is entered by the user as a command line parameter.
    * The testing works by reading a file with values and keys, which will be entered by the user in the command line
    *  when running the program, and the file will be read using file redirection.
    */
    static class heap {
        //array that will keep the values of the heap
        private ArrayList<Node> array = new ArrayList<>();
        //size of the heap
        private int heapSize = 0;
        //Log base 2 of the branching factor (used to speed up the calculations of the indexes of the parent/children)
        private int lgBranchingFactor;
        //counter to track the number of key comparisons
        private int keyComparisons = 0;

        private Integer[] vertexLocations = new Integer[n];

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
        static class Node {
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

            vertexLocations[node.value] = i;
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
            
            vertexLocations[min.value] = null;
            
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
            Node srcNode = array.get(src);
            Node destNode = array.get(dest);

            array.set(dest, srcNode);
            array.set(src, destNode);

            vertexLocations[srcNode.value] = src;
            vertexLocations[srcNode.value] = dest;
        }
	}

   }