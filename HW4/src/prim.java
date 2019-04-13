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

        Graph g = new Graph(adjList, vertices);
        new prim(g, n, m, branchingFactor).run();
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
     *
     */
    private void run() {
        for (Vertex v : g.getVertices()) {
            if (!v.partOfSpanningTree) {
                mstPrim(g, v);
            }
        }

        long numTrees = calculateNumTrees(g);
        long totalWeight = calculateTotalWeight(g);

        System.out.println(String.format("%s %s %s", branchingFactor, numTrees, totalWeight));
    }

    /**
     * Prim's algorithm
     *
     * @param g Graph
     * @param root root vertex to start Prim's algorithm
     */
    private void mstPrim(Graph g, Vertex root) {
        MinHeap priorityQueue = new MinHeap(branchingFactor, n);
        root.distance = 0;
        root.partOfSpanningTree = true;

        priorityQueue.insertValue(new Pair(root.distance, root.id));

        Vertex[] vertices = g.getVertices();
        for(Vertex v : vertices) {
            if (!v.partOfSpanningTree) {
                priorityQueue.insertValue(new Pair(v.distance, v.id));
            }
        }

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

        private Vertex(Integer id, Integer distance) {
            this.id = id;
            this.distance = distance;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            Vertex vertex = (Vertex) o;
            return Objects.equals(id, vertex.id);
        }

        @Override
        public int hashCode() {
            return Objects.hash(id);
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


    // Key/value pair, used in the heap and its interface.
    static class Pair {
        public int key;
        public int val;

        public Pair( int k, int v ) {
            key = k;
            val = v;
        }
    }

    // Actual representation of the heap
    static class MinHeap {
        // Counter of comparison operations, for comparing performance.
        private long comparisons = 0;

        // Power of 2 used as the branchign factor
        private int p;

        // Representation for the heap.
        Pair[] tree;

        // Number of elements in the heap.
        int n;

        // Capacity of the heap.
        int cap;

        Integer[] vertexLocations;

        // Function to compare keys, so we can also count key comparisons.
        private boolean keyLess( Pair a, Pair b ) {
            comparisons += 1;
            return a.key < b.key;
        }

        public MinHeap( int p, int numVertices ) {
            this.p = p;
            cap = 5;
            n = 0;
            tree = new Pair [ cap ];
            vertexLocations = new Integer[numVertices];
        }

        Pair removeMin() {
            // Remove the minimum value and replace it with the last one.
            Pair v = tree[ 0 ];
            tree[ 0 ] = tree[ n - 1 ];
            n -= 1;

            vertexLocations[v.val] = null;
            vertexLocations[tree[0].val] = 0;

            // We need the branching factor below.
            int branch = 1 << p;

            // Push this value down until it satisfies the ordering constraint.
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

                // Not happy about this early return.  Would be nice to ahve it in the condition
                // on the loop.  Return early if we hit a point where we don't have to swap.
                if ( ! keyLess( tree[ m ], tree[ idx ] ) )
                    return v;

                // Swap the current vlaue with its smallest child
                Pair temp = tree[ idx ];
                tree[ idx ] = tree[ m ];
                tree[ m ] = temp;

                vertexLocations[tree[idx].val] = idx;
                vertexLocations[tree[m].val] = m;

                // Follow the value down into the tree.
                idx = m;
                child = ( idx << p ) + 1;
            }

            return v;
        }

        void insertValue( Pair v ) {
            if ( n >= cap ) {
                // Enlarge the heap array and copy everything over.
                cap *= 2;
                Pair[] t2 = new Pair [ cap ];
                for ( int i = 0; i < n; i++ )
                    t2[ i ] = tree[ i ];
                tree = t2;
            }

            // Put the new value at the end of the heap.
            int idx = n;
            tree[ n ] = v;
            vertexLocations[tree[n].val] = n;

            n++;

            // Move it up in the heap until it's as large as its parent.
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

        /** Return the number of comparisons performed. */
        long ccount() {
            return comparisons;
        }

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

        boolean contains(int vertexId) {
            return vertexLocations[vertexId] != null;
        }

        int parent(int idx) {
            return (idx - 1) >> p;
        }
    }
}
