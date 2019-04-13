/*
 * Import statements
 */
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Scanner;
import java.util.Set;

/**
 * Outputs the the number of trees in the given list of edges and the total weight of all
 * the minimum-weight spanning trees in the forest, using Prim's algorithm with a heap as
 * a priority queue.
 *
 * @author Jeremy Maness
 * @author Gabriel Oliveira
 */
public class prim {
    private Graph g; // Graph to run Prim's algorithm on
    private int n; // Number of vertices
    private int m; // Number of edges
    private int branchingFactor; //Branching factor to use in the heap that will be used as a priority queue

    /**
     * Initializes the object used in the implementation
     * @param g Graph to run Prim's algorithm on
     * @param n Number of vertices in graph g
     * @param m Number of edges in graph g
     * @param branchingFactor // Branching factor based on n and m that will be used in the heap
     */
    private prim(Graph g, int n, int m, int branchingFactor) {
        this.g = g;
        this.n = n;
        this.m = m;
        this.branchingFactor = branchingFactor;
    }

    /**
     * Reads input file.
     * Calculates branching factor.
     * Fills up the adjacency list and sets of vertices from the given input file.
     * Initialized the Graph object and runs the Prim's algorithm
     *
     * @param args arguments on user input
     */
    public static void main(String[] args) {

        /*
         * Reading input
         */
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

        /*
         * Initializing and filling adjacency list and set of vertices from the graph
         * using information from the input file
         */
        List<Set<Neighbor>> adjList = new ArrayList<>(n);
        Set<Vertex> vertices = new HashSet<>();
        
        for (int i = 0; i < n; i++) {
            adjList.add(null);
            vertices.add(new Vertex(i, Integer.MAX_VALUE));
        }

        // Reading edges given on the input
        while (scanner.hasNext()) {
            int a = scanner.nextInt(); // First vertex of edge
            int b = scanner.nextInt(); // Second vertex of edge
            int w = scanner.nextInt(); // Weight of the edge

            vertices.add(new Vertex(a, Integer.MAX_VALUE));
            vertices.add(new Vertex(b, Integer.MAX_VALUE));

            if (adjList.get(a) == null) {
                adjList.set(a, new HashSet<>());
            }

            if (adjList.get(b) == null) {
                adjList.set(b, new HashSet<>());
            }

            adjList.get(a).add(new Neighbor(b, w));
            adjList.get(b).add(new Neighbor(a, w));
        }

        Graph g = new Graph(adjList, vertices); //Initializes graph
        new prim(g, n, m, branchingFactor).run(); //Runs Prim's algorithm
    }

    /**
     * Calculates the branch factor given the number of vertices and the number of edges
     *
     * @param vertices Number of vertices
     * @param edges Number of edges
     * @return branching factor
     */
    private static int calculateBranchingFactor(int vertices, int edges) {
        return Math.max(2, (int) Math.pow(2, Math.ceil(Math.log((float) edges / vertices) / Math.log(2))));
    }

    /**
     * Executes the forest version of Prim's algorithm
     * Calculates total number of trees in the forest and
     * total weight across all trees.
     * Prints output on the console.
     */
    private void run() {
        
        /*
         * Running Prim's algorithms for every vertex that
         * is not part of a spanning tree
         */
        for (Vertex v : g.getVertices()) {
            if (!v.partOfSpanningTree) {
                mstPrim(g, v);
            }
        }

        long numTrees = calculateNumTrees(g); //Calculating number of trees
        long totalWeight = calculateTotalWeight(g); //Calculating total weight of all trees

        System.out.println(String.format("%s %s %s", branchingFactor, numTrees, totalWeight)); //Printing output
    }

    /**
     * Prim's algorithm implementation
     *
     * @param g Graph to run Prim's algorithm on
     * @param root root vertex to start Prim's algorithm
     */
    private void mstPrim(Graph g, Vertex root) {
        
        /*
         * Initializing the priority queue.
         * Setting given vertex as root of the tree.
         */
        MinHeap priorityQueue = new MinHeap(branchingFactor, n);
        root.distance = 0;
        root.partOfSpanningTree = true;

        priorityQueue.insertValue(new Pair(root.distance, root.id)); //Inserting root int the priority queue

        /*
         * Creating array with all the vertices of graph g.
         * This is done to be able to access the vertices in constant time
         */
        Vertex[] vertices = g.getVertices(); 
        
        /*
         * Inserting each vertex that is not part of a spanning tree
         * into the heap
         */
        for(Vertex v : vertices) {
            if (!v.partOfSpanningTree) {
                priorityQueue.insertValue(new Pair(v.distance, v.id));
            }
        }

        /*
         * While there are vertices in the heap, go through the adjacency list of the vertex u with the lowest key
         * and change the parent, the edge weight and decrease the key of vertex v if the current edge eight is greater 
         * than the distance between u and v.
         */
        while (priorityQueue.n != 0) {
            Pair u = priorityQueue.removeMin();
            Set<Neighbor> neighbors = g.getAdjVertices(u.val);
              
            if (neighbors != null) {
                for (Neighbor v : neighbors) {
                    if (priorityQueue.contains(v.vertexId) && v.edgeWeight < vertices[v.vertexId].distance) {
                        vertices[v.vertexId].parent = vertices[u.val];
                        vertices[v.vertexId].distance = v.edgeWeight;
                        vertices[v.vertexId].partOfSpanningTree = true;
                        priorityQueue.decreaseKey(v.vertexId, v.edgeWeight);
                    }
                }
            }
        }
    }

    /**
     * Counts the number of trees in the forest.
     * The runtime complexity is O(n).
     *
     * @param g graph
     * @return number of trees in the forest
     */
    private long calculateNumTrees(Graph g) {
        return Arrays.stream(g.vertices)
                .filter(vertex -> vertex.parent == null)
                .count();
    }

    /**
     * Counts the total weight across all trees in the forest.
     * The runtime complexity is O(n).
     *
     * @param g graph
     * @return total edge weight of the  forest
     */
    private long calculateTotalWeight(Graph g) {
        return Arrays.stream(g.vertices)
                .mapToLong(vertex -> vertex.distance)
                .sum();
    }

    /**
     * Graph that uses an adjacency list representation
     *
     */
    static class Graph {
        private ArrayList<Set<Neighbor>> adjList; //Graph's adjacency list
        private Vertex[] vertices; // Set of all the vertices in the graph

        /**
         * Initializes the grapph with given adjacency list and vertices
         * @param adjList Adjacency list of all vertices in the graph
         * @param vertices Given set of all vertices on the graph
         */
        private Graph(List<Set<Neighbor>> adjList, Set<Vertex> vertices) {
            this.adjList = new ArrayList<>(adjList);
            this.vertices = new Vertex[vertices.size()];

            /*
             * Copying the given set of vertices into an array
             */
            for (Vertex v : vertices) {
                this.vertices[v.id] = v;
            }
        }

        /**
         * Returns adjacency list for a given vertex
         * @param u Vertex to get the adjacency list from
         * @return adjacency list of the given vertex
         */
        Set<Neighbor> getAdjVertices(int u) {
            return adjList.get(u);
        }

        /**
         * Returns an array of Vertex objects containing all the vertices on the graph 
         * @return all the vertices of the graph
         */
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

        private Integer id; //Vertex ID
        private Integer distance; //cost of the edge between this vertex and its parent
        private Vertex parent; //Parent of this vertex in a tree
        private boolean partOfSpanningTree; //Boolean flag to show if this vertex is part of a spanning tree 

        /**
         * Initializes this vertex object
         * @param id ID for the vertex
         * @param distance Distance between this vertex and its parent
         */
        private Vertex(Integer id, Integer distance) {
            this.id = id;
            this.distance = distance;
        }

        /**
         * Compares this vertex with the given object
         * @param o object to be compared
         * @return True if o and this vertex are the same
         */
        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            Vertex vertex = (Vertex) o;
            return Objects.equals(id, vertex.id);
        }

        /**
         * Vertex hash code
         */
        @Override
        public int hashCode() {
            return Objects.hash(id);
        }
    }

    /**
     * Neighbor object
     * It represents a vertex on the adjacency list of another vertex from the graph
     */
    static class Neighbor {
        private int vertexId; //ID of this vertex
        private int edgeWeight; //Weight on the edge between this vertex and its parent

        /**
         * Initializes this object
         * @param vertexId ID for this vertex
         * @param edgeWeight Weight on the edge between this vertex and its parent
         */
        private Neighbor(int vertexId, int edgeWeight) {
            this.vertexId = vertexId;
            this.edgeWeight = edgeWeight;
        }

        /**
         * Neighbor hash code
         */
        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + edgeWeight;
            result = prime * result + vertexId;
            return result;
        }

        /**
         * Compares this neighbor with the given object
         * @param obj object to be compared
         * @return True if obj and this neighbor are the same
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
            if (vertexId != other.vertexId)
                return false;
            return true;
        }
    }


    /**
     * Key/value pair, used in the heap and its interface.
     * @author Dr. Sturgill
     */
    static class Pair {
        public int key; //Key to organize the heap
        public int val; //Value of the heap node

        /**
         * Initializes a pair with key and value
         * @param k Given key
         * @param v Given value
         */
        public Pair( int k, int v ) {
            key = k;
            val = v;
        }
    }

    /**
     * Actual representation of the heap.
     * 
     * @author Dr. Sturgill
     * @author Jeremy Maness
     */
    static class MinHeap {
       
        private long comparisons = 0; // Counter of comparison operations, for comparing performance.
        private int p; // Power of 2 used as the branching factor
        Pair[] tree; // Representation for the heap.
        int n; // Number of elements in the heap.
        int cap; // Capacity of the heap.
        Integer[] vertexLocations; //array of vertices location to link vertices between the heap and the graph

        /**
         * Function to compare keys, so we can also count key comparisons.
         * 
         * @param a Heap node
         * @param b Heap node
         * @return True if the key of node a is less than the key of node b
         */
        private boolean keyLess( Pair a, Pair b ) {
            comparisons += 1;
            return a.key < b.key;
        }

        /**
         * Initializes the heap
         * @param p Power of 2 used as the branching factor
         * @param numVertices Maximum number of vertices that can be stored in the heap
         */
        public MinHeap( int p, int numVertices ) {
            this.p = p;
            cap = 5; //Initial capacity
            n = 0; //Number of vertices
            tree = new Pair [ cap ]; //Representation of the tree
            vertexLocations = new Integer[numVertices];
        }

        /**
         * Remove the minimum value from the heap and executes a heapify operation to reorganize the heap.
         * @return the vertex with minimum key value of the heap
         */
        Pair removeMin() {
            // 
            Pair v = tree[ 0 ];
            tree[ 0 ] = tree[ n - 1 ];
            n -= 1;

            vertexLocations[v.val] = null;
            vertexLocations[tree[0].val] = 0;

            /*
             * We need the branching factor below.
             */
            int branch = 1 << p;

            /*
             * Push this value down until it satisfies the ordering constraint.
             */
            int idx = 0;
            int child = ( idx << p ) + 1;
            while ( child < n ) {
                // Find index of smallest child.
                int m = child;
                int end = child + branch;
                if ( end > n )
                    end = n;
                for ( int i = child + 1; i < end; i++ )
                    if ( keyLess( tree[ i ], tree[ m ] ) )
                        m = i;

                /*
                 * Not happy about this early return.  Would be nice to have it in the condition
                 * on the loop.  Return early if we hit a point where we don't have to swap.
                 */
                if ( ! keyLess( tree[ m ], tree[ idx ] ) )
                    return v;

                /*
                 * Swap the current value with its smallest child
                 */
                Pair temp = tree[ idx ];
                tree[ idx ] = tree[ m ];
                tree[ m ] = temp;

                vertexLocations[tree[idx].val] = idx;
                vertexLocations[tree[m].val] = m;

                /*
                 * Follow the value down into the tree.
                 */
                idx = m;
                child = ( idx << p ) + 1;
            }

            return v;
        }

        /**
         * Insert the given vertex node into the heap and then perfomes a heapify operation to reorganize
         * the heap
         * @param v Vertex node
         */
        void insertValue(Pair v) {
            
            /*
             * If the number of vertices in the graph is larger than the capacity of the heap
             * creates a new tree representation with twice the original capacity, copy all
             * the nodes to the new one and then replace the orginal one.
             */
            if ( n >= cap ) {
                cap *= 2;
                Pair[] t2 = new Pair [ cap ];
                for ( int i = 0; i < n; i++ )
                    t2[ i ] = tree[ i ];
                tree = t2;
            }

            
            /**
             * Put the new value at the end of the heap.
             */
            int idx = n;
            tree[ n ] = v;
            vertexLocations[tree[n].val] = n;

            n++;

            //Move it up in the heap until it's as large as its parent.
            int par = ( idx - 1 ) >> p;
            while ( par >= 0 && keyLess( tree[ idx ], tree[ par ] ) ) {
                // Swap this value with its parent.
                Pair temp = tree[ par ];
                tree[ par ] = tree[ idx ];
                tree[ idx ] = temp;

                vertexLocations[tree[par].val] = par;
                vertexLocations[tree[idx].val] = idx;

                idx = par;
                par = ( idx - 1 ) >> p;
            }
        }

        /**
         *  Return the number of comparisons performed.
         */
        long ccount() {
            return comparisons;
        }

        /**
         * Decreases the key of the given vertex in the heap
         * @param vertexId Given vertex ID
         * @param key Original vertex key value
         */
        void decreaseKey(int vertexId, int key) {
            Integer idx = vertexLocations[vertexId];
            Pair pair = tree[idx];

            if (key > pair.key) {
                throw new RuntimeException("new key is larger than current key");
            }

            pair.key = key;

            int i = idx;
            while (i > 0 && tree[parent(i)].key > pair.key) {
                Pair temp = tree[i];
                tree[i] = tree[parent(i)];
                tree[parent(i)] = temp;

                vertexLocations[tree[parent(i)].val] = parent(i);
                vertexLocations[tree[i].val] = i;

                i = parent(i);
            }
        }

        /**
         * Search the heap for the given vertex
         * @param vertexId ID of the vertex to search
         * @return True if the vertex is in the heap
         */
        boolean contains(int vertexId) {
            return vertexLocations[vertexId] != null;
        }

        /**
         * Returns the parent of the vertex with given index
         * @param idx Vertex index
         * @return Index of the parent
         */
        int parent(int idx) {
            return (idx - 1) >> p;
        }
    }
}
