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
        root.distance = 0;
        root.partOfSpanningTree = true;


        Vertex[] vertices = g.getVertices();
        List<Pair> pairs = new ArrayList<>();

        pairs.add(new Pair(root.distance, root.id));

        /*
         * Inserting each vertex that is not part of a spanning tree
         * into the heap
         */
        for(Vertex v : vertices) {
            if (!v.partOfSpanningTree) {
                pairs.add(new Pair(v.distance, v.id));
            }
        }

        MinHeap priorityQueue = new MinHeap(branchingFactor, n, pairs);

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
    class Pair {
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
    class MinHeap {
        private int p; // Power of 2 used as the branching factor
        private Pair[] tree; // Representation for the heap.
        private int n; // Number of elements in the heap.
        private int cap; // Capacity of the heap.
        private Integer[] vertexLocations; //array of vertices location to link vertices between the heap and the graph

        /**
         * Initializes the heap
         * @param p Power of 2 used as the branching factor
         * @param numVertices Maximum number of vertices that can be stored in the heap
         */
        public MinHeap( int p, int numVertices, List<Pair> pairs ) {
            this.p = p;
            cap = pairs.size(); //Initial capacity
            n = pairs.size(); //Number of vertices
            tree = pairs.toArray(new Pair[0]);
            vertexLocations = new Integer[numVertices];

            for (int i = 0; i < tree.length; i++) {
                vertexLocations[tree[i].val] = i;
            }

            buildMinHeap();
        }

        /**
         * Implements the Build-Min-Heap procedure from "Introduction to Algorithms, Third Edition"
         * (Cormen et al. 2009, p. 157).
         *
         * As shown on p. 157-159, a similar asymptotic analysis shows that this is O(V).
         *
         */
        private void buildMinHeap() {
            n = tree.length;
            for (int i = tree.length / 2; i >= 0; i--) {
                minHeapify(i);
            }
        }

        /**
         * Min heapify procedure that pushes the node at the specified index in
         * the heap down until the heap ordering constraint is satisfied.
         *
         * @param idx index of node in the heap
         */
        private void minHeapify(int idx) {

            /*
             * We need the branching factor below.
             */
            int branch = 1 << p;

            /*
             * Push this value down until it satisfies the ordering constraint.
             */
            int child = ( idx << p ) + 1;
            while ( child < n ) {
                // Find index of smallest child.
                int m = child;
                int end = child + branch;
                if ( end > n )
                    end = n;
                for ( int i = child + 1; i < end; i++ )
                    if (tree[i].key < tree[m].key)
                        m = i;

                /*
                 * Not happy about this early return.  Would be nice to have it in the condition
                 * on the loop.  Return early if we hit a point where we don't have to swap.
                 */
                if (tree[m].key >= tree[idx].key)
                    return;

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
        }

        /**
         * Remove the minimum value from the heap and executes a heapify operation to reorganize the heap.
         *
         * @return the vertex with minimum key value of the heap
         */
        Pair removeMin() {
            Pair v = tree[ 0 ];
            tree[ 0 ] = tree[ n - 1 ];
            n -= 1;

            vertexLocations[v.val] = null;
            vertexLocations[tree[0].val] = 0;

            minHeapify(0);

            return v;
        }

        /**
         * Decreases the key of the given vertex in the heap
         *
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
         *
         * @param vertexId ID of the vertex to search
         * @return True if the vertex is in the heap
         */
        boolean contains(int vertexId) {
            return vertexLocations[vertexId] != null;
        }

        /**
         * Returns the parent of the vertex with given index
         *
         * @param idx Vertex index
         * @return Index of the parent
         */
        int parent(int idx) {
            return (idx - 1) >> p;
        }
    }
}
