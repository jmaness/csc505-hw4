import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.Set;

import static java.lang.Integer.MAX_VALUE;

/**
 * Outputs the the number of trees in the given list of edges and the total weight of all
 * the minimum-weight spanning trees in the forest, using Prim's algorithm with a heap as 
 * a priority queue.
 *
 * @author Jeremy Maness
 * @author Gabriel Oliveira
 */
public class prim {
    private Graph g;
    private int n; // Number of vertices
    private int m; // Number of edges
    private int branchingFactor;

    private prim(Graph g, int n, int m, int branchingFactor) {
        this.g = g;
        this.n = n;
        this.m = m;
        this.branchingFactor = branchingFactor;
    }

    /**
     * Reads input file and fill up the adjacency matrix.
     * Calculates branching factor
     *
     * @param args arguments on user input
     */
    public static void main(String[] args) {

        // Reading input
        Scanner scanner = new Scanner(System.in);

        int n = -1;
        int m = -1;

        // The first number in the stream is the number of vertices
        if (scanner.hasNextInt()) {
            n = scanner.nextInt();
        }

        // The second number in the stream is the number of edges
        if (scanner.hasNextInt()) {
            m = scanner.nextInt();
        }

        /*
         * Calculating branching factor based number of vertices and edges from input
         */
        int branchingFactor = calculateBranchingFactor(n, m);

        List<Set<Neighbor>> adjList = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            adjList.add(null);
        }

        Set<Vertex> vertices = new HashSet<>();

        // Reading edges given on the input
        while (scanner.hasNext()) {
            int a = scanner.nextInt(); // First vertex of edge
            int b = scanner.nextInt(); // Second vertex of edge
            int w = scanner.nextInt(); // Weight of the edge

            vertices.add(new Vertex(a, Integer.MAX_VALUE, null));
            vertices.add(new Vertex(b, Integer.MAX_VALUE, null));

            Set<Neighbor> neighbors = adjList.get(a);

            if (neighbors == null) {
                adjList.set(a, new HashSet<>());
            }

            adjList.get(a).add(new Neighbor(b, w));
        }

        Graph g = new Graph(adjList, vertices);
        new prim(g, n, m, branchingFactor).run();
    }

    /**
     * Calculates the branch factor given the number of vertices and the number of edges
     * @param vertices Number of vertices
     * @param edges Number of edges
     * @return branching factor
     */
    private static int calculateBranchingFactor(int vertices, int edges) {
        return (int) Math.pow(2, Math.ceil(Math.log((float) edges/vertices) / Math.log(2)));
    }

    /**
     * Executes the forest version of Prim's algorithm
     *
     */
    private void run() {
        for (Vertex v : g.getVertices()) {
            if (!v.partOfSpanningTree) {
                mstPrim(g, v);
            }
        }
    }

    /**
     * Prim's algorithm
     *
     * @param g Graph
     * @param root root vertex to start Prim's algorithm
     */
    private void mstPrim(Graph g, Vertex root) {
        heap priorityQueue = new heap(branchingFactor, n);
        root.distance = 0;
        
        Vertex[] vertices = g.getVertices();
        for(Vertex v : vertices) {
            priorityQueue.insertValue(new prim.heap.Node(v.distance, v.id));
        }

        while (priorityQueue.heapSize != 0) {
            prim.heap.Node u = priorityQueue.removeMin();

            Set<Neighbor> neighbors = g.getAdjVertices(u.value);
            if (neighbors != null) {
                for (Neighbor v : g.getAdjVertices(u.value)) {
                    if (priorityQueue.contains(v.vertexId) && v.edgeWeight < vertices[v.vertexId].distance) {
                        vertices[v.vertexId].parent = vertices[u.value];
                        vertices[v.vertexId].distance = v.edgeWeight;
                        priorityQueue.decreaseKey(v.vertexId, v.edgeWeight);
                    }
                }
            }
        }
    }

    /**
     * Graph that uses an adjacency list representation
     *
     */
    static class Graph {
        private ArrayList<Set<Neighbor>> adjList;
        private Vertex[] vertices;
        
        private Graph(List<Set<Neighbor>> adjList, Set<Vertex> vertices) {
            this.adjList = new ArrayList<>(adjList);
            this.vertices = new Vertex[vertices.size()];

            for (Vertex v : vertices) {
                this.vertices[v.id] = v;
            }
        }

        Set<Neighbor> getAdjVertices(int u) {
            return adjList.get(u);
        }
        
        Vertex[] getVertices(){
            return this.vertices;
        }
    }

    /**
     * Vertex that maintains:
     *  - a minimum weight of any edge connecting v to any vertex in the spanning tree
     *  - parent vertex in a spanning tree
     *  - flag to indicate that the vertex is part of a spanning tree
     *
     */
    static class Vertex {

        private Integer id;
        private Integer distance;
        private Vertex parent;
        private boolean partOfSpanningTree;

        private Vertex(Integer id, Integer distance, Vertex parent) {
            this.id = id;
            this.distance = distance;
            this.parent = parent;
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

    /**
     * Neighbor (Adjacency List nodes)
     *
     */
    static class Neighbor {
        private int vertexId;
        private int edgeWeight;

        private Neighbor(int vertexId, int edgeWeight) {
            this.vertexId = vertexId;
            this.edgeWeight = edgeWeight;
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + edgeWeight;
            result = prime * result + vertexId;
            return result;
        }

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
            if (vertexId != other.vertexId)
                return false;
            return true;
        }
    }



    ///////////////////////////////////////// Tree inner class from homework 2 //////////////////////////////////////////
    public class Tree {

    }




    ///////////////////////////////////////// Heap inner class from homework 2 //////////////////////////////////////////

    /**
     * This class constructs and tests a min-heap in which the nodes can have
     * an arbitrary number of children, as long as the number is a power of 2.
     *.
     */
    static class heap {
        private ArrayList<Node> array = new ArrayList<>(); //array that will keep the values of the heap
        private int heapSize = 0; // size of the heap
        private int branchingFactor;

        // Log base 2 of the branching factor (used to speed up the calculations of the indexes of the parent/children)
        private int lgBranchingFactor;

        private Integer[] vertexLocations;

        /**
         * Heap constructor.
         *
         * @param branchingFactor branching factor that is a power of 2 for efficient heap-ification.
         * @param maxVertices maximum number of vertices that can possibly be in the heap simultaneously
         */
        private heap(int branchingFactor, int maxVertices) {
            this.branchingFactor = branchingFactor;
            this.lgBranchingFactor = (int) Math.floor(Math.log(branchingFactor) / Math.log(2));
            this.vertexLocations = new Integer[maxVertices];
        }

        /**
         * Inner class for creation of the node.
         * A Node object represents a node of the heap, and it is composed
         * by the key and the value.
         */
        static class Node {
            private int key; // Key to be compared
            private int value; // Value stored on the node

            /**
             * Node constructor
             *
             * @param key field to sort by in the heap
             * @param value integer to prioritize in the heap
             */
            Node(int key, int value) {
                this.key = key;
                this.value = value;
            }

            /**
             * Override the toString of the Node class to format the string output
             * of a Node
             */
            @Override
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

        boolean contains(int vertexId) {
            return vertexLocations[vertexId] != null;
        }

        void decreaseKey(int vertexId, int key) {
            int nodeId = vertexLocations[vertexId];
            Node node = array.get(nodeId);
            node.key = key;
            minHeapify(nodeId);
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
            return node1.key < node2.key;
        }

        /**
         * Heapifies this heap starting at the specified index.
         *
         * @param i index to start heap-ification.
         */
        private void minHeapify(int i) {
            int smallest = i;

            for (int j = 1; j <= branchingFactor; j++) {
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

            vertexLocations[destNode.value] = src;
            vertexLocations[srcNode.value] = dest;
        }
    }
}
