import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A class representing a Graph with vertices and adjacent list to represent graph.
 *
 * @author Neehar
 */

public class MST {
    
	private int vertices; 
	private int edges;
	private HashMap<Integer, Integer>[] adjList;
	private ArrayList<Integer> verticesArray; 
	
	
	/**
	 * @return : ArrayList<Integer>.
	 * @param : null
	 * This method returns a list of Vertices of a Graph
	 */
	
	public ArrayList<Integer> getVerticesArray() {
		return verticesArray;
	}
	
	/**
	 * @return : null
	 * @param : ArrayList<Integer>
	 * This method is used to set verticesArray of graph
	 */
	
	public void setVerticesArray(ArrayList<Integer> verticesArray) {
		this.verticesArray = verticesArray;
	}
	
	/**
	 * @return : int
	 * @param : null
	 * This method returns number of edges of a graph
	 */

	public int getEdges() {
		return edges;
	}
    
	/**
	 * @return : null
	 * @param : int (edges)
	 * This method is used to set Edges for  a Graph
	 */
	
	public void setEdges(int edges) {
		this.edges = edges;
	}
	
	/**
	 * @return : Array of HashMap<Integer, Integer>
	 * @param : null
	 * This method returns an array of HashMaps
	 */
	
	public HashMap<Integer, Integer>[] getAdjList() {
		return adjList;
	}
    

	/**
	 * @return : null
	 * @param : Array of HashMap<Integer, Integer>
	 * This method returns an array of HashMaps
	 */
	
	public void setAdjList(HashMap<Integer, Integer>[] adjList) {
		this.adjList = adjList;
	}
	
	/**
	 * A class representing an Edge of a Graph. This edge will have destination vertex and 
	 * cost associated with edge as its members.
	 *
	 * @author Neehar
	 */
	
	static class Edge {

		private int vertex2;
		private int weight;

		/**
		 * @return the vertex2
		 */
		public int getVertex2() {
			return vertex2;
		}

		/**
		 * @param set
		 *            the vertex2 to edge
		 */
		public void setVertex2(int vertex2) {
			this.vertex2 = vertex2;
		}

		/**
		 * @return the weight
		 */
		public int getWeight() {
			return weight;
		}

		/**
		 * @param set
		 *            the weight to graph
		 */
		public void setWeight(int weight) {
			this.weight = weight;
		}

		/**
		 * @param int vertex, int weight This Constructor sets the weight and
		 *        destination vertex
		 */

		public Edge(int vertex2, int weight) {
			this.vertex2 = vertex2;
			this.weight = weight;
		}

	}
	
	/**
	 * @param null
	 *   This is the default constructor for the graph.
	 */

	public MST(){
		
	}
	
	/**
	 * @param int vertices
	 * @return Graph with entered vertices
	 * This is the parameterized constructor for the graph.
	 * This creates a adjacency list representation of graph with the given 
	 * vertices.
	 */
	
	@SuppressWarnings("unchecked")
	public MST(int vertices) {
		this.vertices= vertices;
		this.edges =0;
		adjList = new HashMap[vertices];
		for (int i=0; i<vertices; i++)
			adjList[i] = (new HashMap<Integer, Integer>());
		
		verticesArray = new ArrayList<Integer>(vertices);
		for (int i=0; i<vertices; i++)
			verticesArray.add(i);
	}
	
	
	/**
	 * @param int vertices, int edges
	 * @return Graph with entered vertices and edges
	 * This is the parameterized constructor for the graph.
	 * This creates a adjacency list representation of graph with the given 
	 * vertices and edges.
	 */
	@SuppressWarnings("unchecked")
	public MST(int vertices, int edges) {
		this.vertices = vertices;
		this.edges =edges;
		adjList = new HashMap[vertices];
		for (int i=0; i<vertices; i++)
			adjList[i] = (new HashMap<Integer, Integer>());
		
		verticesArray = new ArrayList<Integer>(vertices);
		for (int i=0; i<vertices; i++)
			verticesArray.add(i);
	}
	
	/**
	 * @param String
	 * @return returns a graph object with edges read from a file placed in a location.
	 * This is the parameterized constructor for the graph.
	 * This creates a adjacency list representation of graph with the given 
	 * vertices and edges.
	 */
	@SuppressWarnings("unchecked")
	public MST(String fileName) {
		BufferedReader bReader = null;
		try {
			File file = new File(fileName);
			bReader = new BufferedReader(new FileReader(file));
			String firstLine = bReader.readLine();
			this.vertices = Integer.valueOf(firstLine.split(" ")[0]);
			this.edges = Integer.valueOf(firstLine.split(" ")[1]);
			this.verticesArray = new ArrayList<Integer>();
			for(int i = 0;i<vertices;i++){
				verticesArray.add(i);
			}
			adjList = new HashMap[vertices];
			for(int i=0; i < adjList.length; i++){
				adjList[i] = new HashMap<Integer, Integer>();
			}
			String nextLine = null;
			while(null!=(nextLine=bReader.readLine())){
				int vertex1 = Integer.valueOf(nextLine.split(" ")[0]);
				int vertex2 = Integer.valueOf(nextLine.split(" ")[1]);
				int weight = Integer.valueOf(nextLine.split(" ")[2]);
				addEdge(vertex1, vertex2, weight);
			}
			 if(!isConnected(dfs(0, new boolean[vertices]))){
				 System.out.println("Graph is not connected");
				 System.exit(0);
			 }
		} catch (IOException e) {
			System.out.println("Exception : "+e.getMessage());
		} finally{
			try {
				bReader.close();
			} catch (IOException e) {
				System.out.println("Exception : "+e.getMessage());
			}
		}
	}
    
	/**
	 * @param int vertex1, int vertex2, int weight
	 * @return null
	 * Adds an edge to the graph as graph is bidirectional.
	 * */
	public void addEdge(int vertex1, int vertex2, int weight) {
		if(!(vertex1==0 && vertex2==0)){
		adjList[vertex1].put(vertex2,weight);
		adjList[vertex2].put(vertex1,weight);
		}else{
			adjList[0].put(0, 0);
		}
	}
	
	/**
	 * @param int sourceVertex
	 * @return Gets the linked list of a vertex
	 * This returns the HashMap of the vertices with corresponding weights from the source vertex.
	 * 
	 * */
	public HashMap<Integer, Integer> adjList(int v) {
		return adjList[v];
	}
    
	/**
	 * @param null
	 * @return number of vertices of the graph
	 * This returns the number of the vertices of Graph. 
	 * 
	 * */
	 
	public int vertices() {
		return vertices;
	}
    
	/**
	 * @param null
	 * @return null
	 * This method displays the graph. Displays all the destination vertices from a vertex and weight 
	 * associated with edge.  
	 * 
	 * */
	
	public void displayAdjList() {
		for (int i=0; i<vertices; i++) {
			System.out.print(verticesArray.get(i)+" -> ");
			System.out.print(adjList[i].entrySet());
			System.out.println();
		}
	}
	
	/**
	 * @param int sourcevertex, boolean[] visited
	 * @return boolean []
	 * This method does the Depth First Search on a Graph. By doing so 
	 * it fills the boolean array with the value true if a node is visited during 
	 * traversal of the graph.   
	 * 
	 * */
	
	public  boolean[] dfs(int source, boolean[] visited) {
		visited[source] = true;
		for(int i : this.getAdjList()[source].keySet()){
			if(!visited[i]){
				dfs(i,visited);
			}
		}
		return visited;
	}
	
	/**
	 * @param MST WeightedGraph, int sourceVertex
	 * @return Integer []
	 * This method returns the neighbours assiociated with the sourceVertex.
	 * 
	 * */
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public Integer[] getNeighbours(MST weg, int sourceVertex){
		Integer[] neighbourList = null;
		neighbourList = new Integer[weg.getAdjList()[sourceVertex].size()];
		List iList = new ArrayList();
		for(int i:weg.getAdjList()[sourceVertex].keySet()){
			iList.add(i);
		}
		neighbourList = (Integer[]) iList.toArray(neighbourList);
		return neighbourList;
	}
	
	/**
	 * @param MST WeightedGraph, int sourceVertex
	 * @return Integer []
	 * This method returns the weights assiociated with the neighbours of a sourceVertex.
	 * 
	 * */
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public Integer[] getWeights (MST weg, int sourceVertex){
		Integer[] weightsList = null;
		weightsList = new Integer[weg.getAdjList()[sourceVertex].size()];
		List iList = new ArrayList();
		for(int i:weg.getAdjList()[sourceVertex].values()){
			iList.add(i);
		}
		weightsList = (Integer[]) iList.toArray(weightsList);
		return weightsList;
	}
	
	/**
	 * @param int sourceVertex, int destinationVertex
	 * @return int
	 * This method returns the weights assiociated between set of vertices.
	 * 
	 * */
	public int edgeWeight(int vertex1, int vertex2) {
		return adjList[vertex1].get(vertex2);
	}
	
	/**
	 * @param int k
	 * @return int
	 * Generate random number in the range [0 to k]
	 * 
	 * */
	public static int random(int k){
		int min =0;
		int max =k-1;
		return min + (int)(Math.random()*((max-min)+1));
	}
	
	/**
	 * @param int numberof vertices, double density
	 * @return WeightedGraph
	 * Generate random connected graph with random vertices and random weights
	 * till a connected graph is generated. Graph is generated by calling another 
	 * function which does the Graph generation.
	 * 
	 * */
	public static MST generateRandomConnectedGraph(int vertices,
                                                   double density) {
		long starttime = System.currentTimeMillis(); // have to remove
		while (true) {
			MST weightedGraph = generateRandomGraph(vertices, density);
			boolean[] visited = new boolean[weightedGraph.vertices];
			visited = weightedGraph.dfs(0, visited);
			if (isConnected(visited)) {
				long endTime = System.currentTimeMillis();
				return weightedGraph;
			}
		}
	}
	
	/**
	 * @param int numberof vertices, double density
	 * @return WeightedGraph
	 * Generate graph with random vertices and random weights.
	 * If the new edge generated by random function is already present in the graph,
	 * that edge is not added to the graph.
	 * */
	
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public static MST generateRandomGraph(int vertices, double density){
		
		int edges = (int)(density * vertices * (vertices -1) / 200);
		
		MST weightedGraph = new MST(vertices,edges);
		weightedGraph.setAdjList(new HashMap[vertices]);
		for(int i=0; i < weightedGraph.getAdjList().length; i++){
			weightedGraph.getAdjList()[i] = new HashMap();
		}
		int counter = edges;
	
		do {
			int vertex1 = random(vertices);
			int vertex2 = random(vertices);
			int cost = random(1000) + 1;
			if(vertex1 == vertex2){
				continue;
			}
			if(vertex1<vertex2){
				if (!hasEdge(weightedGraph, vertex1,  vertex2)) {
					weightedGraph.addEdge(vertex1, vertex2 ,cost);
					counter--;
					//noOfEdgesGenerated++;
				}
			}else{
                continue;
			}
            }while(counter>0);
		
		return weightedGraph;
	}
	
	/**
	 * @param Boolean[]
	 * @return boolean
	 * This method checks if the graph is connected or not.
	 * */
	
	public static boolean isConnected(boolean[] visited){
		int flag =0;
		for(int i=0; i < visited.length;i++){
			if(visited[i]==true){
				flag++;
			}
		}
		if(flag==visited.length){
			return true;
		}else{
			return false;
		}
	}
	
	/**
	 * @param  MST weightedGraph, int sourceVertex, int destination vertex
	 * @return Boolean 
	 * This method checks if the edge between the given vertices is present or not in the graph.
	 * */
	public static boolean hasEdge(MST weg,int vertex1, int vertex2){
        
		boolean flag = false;
		int v1 = vertex1;
		int v2 = vertex2;
        
        if(weg.getAdjList()[vertex1].isEmpty()){
            return false;
        }else{
            for (int i=0; i < weg.getAdjList()[v1].size();i++){
                if(weg.getAdjList()[v1].containsKey(v2)){
                    flag = true;
                    break;
                }else{
                    flag = false;
                }
            }
        }
		if(flag == true){
			return true;
		}else{
			return false;
		}
	}
	
	/**
	 * @param  MST weightedGraph, int sourceVertex
	 * @return null 
	 * This method displays the minimum spanning edges of the given grap with the minimum cost
	 * . This is implementation using Array.
	 * */
	public static void displayPrimSimple(MST weightedGraph, int source){
		int totSum =0;
		System.out.println("Simple Scheme list of edges");
		int[] mstEdges = primSimple(weightedGraph, source);
		for( int i=1; i <mstEdges.length; i++){
			totSum += weightedGraph.edgeWeight(mstEdges[i], i);
		}
		System.out.println("Total Minimum cost in  Simple Scheme: "+totSum);
		for( int i=1; i <mstEdges.length; i++){
		System.out.println(mstEdges[i] + "-" +i);
		}
	}
	
	/**
	 * @param  MST weightedGraph, int sourceVertex
	 * @return null 
	 * This method displays the minimum spanning edges of the given grap with the minimum cost
	 * . This is implementation using Fibonacci heap.
	 * */
	
	public static void displayPrimFHeap(MST weightedGraph, int source){
		int totSum=0;
		System.out.println("=============================================================================");
		System.out.println("Fibonacci Scheme list of edges ");
		int[] mstEdges = primFHeap(weightedGraph, source);
		for( int i=1; i <mstEdges.length; i++){
			totSum += weightedGraph.edgeWeight(mstEdges[i], i);
		}
		System.out.println("Total Minimum cost in  Fib Scheme : "+totSum);
		for( int i=1; i <mstEdges.length; i++){
			System.out.println(mstEdges[i] + "-" +i);
		}
	}
	
	/**
	 * @param  MST weightedGraph, int sourceVertex
	 * @return int[] 
	 * This method gives the minimum spanning edges of the given grap with the minimum cost
	 * This is implementation using Fibonacci heap.
	 * */
    public static int [] primFHeap (MST weightedGraph, int source) {
        final int [] weights = new int [weightedGraph.vertices];  // shortest WEIGHT to MST
        final int [] MST = new int [weightedGraph.vertices];  // MST Array
        final boolean [] visited = new boolean [weightedGraph.vertices];
        FibonacciHeapNew fHeap = new FibonacciHeapNew();
        FibonacciHeapNew.HeapNode fNode = new FibonacciHeapNew.HeapNode(0, 0);
        fHeap.insert(fNode, 0);
        for (int i=0; i<weights.length; i++) {
            weights[i] = Integer.MAX_VALUE;// replacement for infinity since it is not possible taken max value of Integer class
            visited[i] = false;
        }
        weights[source] = 0;
        MST[0]= -1;
        for (int i=0; i<weights.length-1; i++) {
            FibonacciHeapNew.HeapNode hNode1 = FibonacciHeapNew.extractMin(fHeap);
            if(hNode1==null)
                break;
            int minVertex = hNode1.getValue();
            while(visited[minVertex]==true){
          
                FibonacciHeapNew.HeapNode hNode2 = FibonacciHeapNew.extractMin(fHeap);
                if(hNode2==null)
                    break;
                minVertex = hNode2.getValue();
            }
            visited[minVertex] = true;
            Integer [] neighbours = weightedGraph.getNeighbours(weightedGraph, minVertex); // vertices Array
            Integer [] weightsArray = weightedGraph.getWeights(weightedGraph, minVertex); // weights   
            for (int j=0;j<neighbours.length;j++) {
                
	            final  int neighbour = neighbours[j];
                final int neighbourwt = weightedGraph.edgeWeight(minVertex,neighbour);
                if (visited[neighbour] == false &&  neighbourwt < weights[neighbour]  ) {
                    FibonacciHeapNew.HeapNode hNode = new FibonacciHeapNew.HeapNode(neighbours[j],weights[j]);
                    fHeap.insert(hNode, weightsArray[j]);
	                weights[neighbour] = neighbourwt;
	                MST[neighbour] = minVertex;
                }
            }
        }
        return MST;
    }
	
    
    /**
	 * @param  MST weightedGraph, int sourceVertex
	 * @return int[] 
	 * This method gives the minimum spanning edges of the given grap with the minimum cost
	 * This is implementation using Array.
	 * */
    
	public static int [] primSimple (MST weightedGraph, int source) {
		
        final int [] weights = new int [weightedGraph.vertices];  // shortest WEIGHT to MST
        final int [] MST = new int [weightedGraph.vertices];  // MST Array
        final boolean [] visited = new boolean [weightedGraph.vertices];
        for (int i=0; i<weights.length; i++) {
            weights[i] = Integer.MAX_VALUE;// replacement for infinity since it is not possible taken max value of Integer class
            visited[i] = false;
        }
        weights[source] = 0;
        MST[0]= -1;
        for (int i=0; i<weights.length-1; i++) {
            final  int minVertex = minVertex (weights, visited);
            visited[minVertex] = true;
            Integer [] neighbours = weightedGraph.getNeighbours(weightedGraph, minVertex);
            for (int j=0;j<neighbours.length;j++) {
                final  int neighbour = neighbours[j];
                final int neighbourwt = weightedGraph.edgeWeight(minVertex,neighbour);
                if (visited[neighbour] == false &&  neighbourwt < weights[neighbour] ) {
                    weights[neighbour] = neighbourwt;
                    MST[neighbour] = minVertex;
                }
            }
        }
     
        return MST;
    }
    
	/**
	 * @param  int[] weights, boolean[] visited
	 * @return int 
	 * This method gives the Vertex with the minimum cost associated from a source.
	 * */
    
    private static int minVertex (int [] weights, boolean [] v) {
        int minWt = Integer.MAX_VALUE;
        int minVertex = -1;   // in case of visted vertex
        for (int i=0; i<weights.length; i++) {
            if (v[i]==false && weights[i]<=minWt)
            {
                minVertex=i;
                minWt=weights[i];
            }
        }
        return minVertex;
    }
	
    /**
     * A class representing a Fibonacci heap with HeapNodes.
     *
     * @author Neehar
     */
    
    static  class FibonacciHeapNew {
        
        private HeapNode minPointer = null; // min pointer
        private int listSize = 0; // size
        private static final double oneByLogPhi = 1.0 / Math.log((1.0 + Math.sqrt(5.0)) / 2.0);
        
		/**
		 * @return the minPointer
		 */
		public HeapNode getMinPointer() {
			return minPointer;
		}

		/**
		 * @param minPointer the minPointer to set
		 */
		public void setMinPointer(HeapNode minPointer) {
			this.minPointer = minPointer;
		}

		/**
		 * @return the listSize
		 */
		public int getListSize() {
			return listSize;
		}

		/**
		 * @param listSize the listSize to set
		 */
		public void setListSize(int listSize) {
			this.listSize = listSize;
		} 
        
		 /**
	     * An Inner class representing a Fibonacci heap with many parameters
	     * depicting heap node.
	     *
	     * @author Neehar
	     */
		
        static  class HeapNode { 
        	
            private int     degree ;       // Number of children
            private boolean isMarked; // Whether this node is marked
            private HeapNode right;   // Next value in the list
            private HeapNode left;  // leftious value in the list
            private HeapNode parent; // Parent in the tree, if any.
            private HeapNode child;  // Child node, if any.
            private int value;     // Element being stored here
            private double key; // key
            private int source;
            /**
        	 * @param null
        	 *   This is the default constructor for the HeapNode.
        	 */
            
            public HeapNode(){
                
            }
            
            public HeapNode(int destination, int weight, int source){
            	left = this;
            	right = this;
            	this.value = destination;
            	this.key = weight;
            	this.source = source;
            }
            
            /**
             * @return int
        	 * @param null
        	 *   This is used to get value of a heap node.
        	 */
            
            public int getValue() {
                return value;
            }
            
            /**
             * @return null
        	 * @param int value
        	 *   This is used to set value of a heap node.
        	 */
            public void setValue(int value) {
                this.value = value;
            }
            
            /**
             * @return double
        	 * @param null
        	 *   This is used to get key of a heap node.
        	 */
            public double getKey() {
                return key;
            }
            
            /**
        	 * @param int value, int key
        	 *   This is the overloaded constructor for the HeapNode with given key and value.
        	 */
            
            public HeapNode(int value, int key) {
                left = this;
                right = this;
                this.value = value;
                this.key = key;
            }
            
            
        }
        
        /**
         * @return FibonacciHeap
    	 * @param null
    	 *   This returns a new Fibonacci Heap
    	 */
        public FibonacciHeapNew makeFibHeap(){
            FibonacciHeapNew fHeap = new FibonacciHeapNew();
            fHeap.listSize=0;
            fHeap.minPointer=null;
            return fHeap;
        }
        
        /**
         * @return FibonacciHeapNode 
    	 * @param null
    	 *   This returns a min heap node (pointer) of fibonacci heap. 
    	 */
        public static  HeapNode getMin(FibonacciHeapNew h){
            return h.minPointer;
        }
        
        /**
         * @return FibonacciHeapNode 
    	 * @param HeapNode h, int key 
    	 *   This method inserts node in to the heap and returns heapnode inserted 
    	 */
        public  HeapNode insert(HeapNode hNode, int key){
            hNode.key = key;
            hNode.child = null;
            hNode.parent = null;
            hNode.degree=0;
            hNode.isMarked = false;
            if(null!=minPointer){ // using pointers we traverse through list and insert 
                hNode.left = minPointer;
                hNode.right = minPointer.right;
                minPointer.right = hNode;
                hNode.right.left = hNode;
                if(key < minPointer.key){
                    minPointer= hNode;
                }
            }
            else{
                minPointer = hNode;
            }
            listSize++;
            return hNode;
        }
        
        /**
         * @return FibonacciHeap 
    	 * @param FibonacciHeap h1, Fibonacciheap h2
    	 *   This method  joins two fibonacci heaps and returns a new fibonacci heap 
    	 */
        public static FibonacciHeapNew heapUnion (FibonacciHeapNew heap1, FibonacciHeapNew heap2){
            FibonacciHeapNew unionHeap = new FibonacciHeapNew();
            if(null!=heap1&& null!=heap2){
                unionHeap.minPointer = heap1.minPointer;
                
                if(unionHeap.minPointer!=null){ // concatenation of root list of heap 2 with union heap root list
                    
                    if(null!=heap2.minPointer){ // first do null checks for both heaps and then implement concatenation of lists
                        
                        unionHeap.minPointer.right.left = heap2.minPointer.left;
                        heap2.minPointer.left.right = unionHeap.minPointer.right;
                        unionHeap.minPointer.right = heap2.minPointer;
                        heap2.minPointer.left = unionHeap.minPointer;
                    }
                    
                    if(heap2.minPointer.key < heap1.minPointer.key){
                        unionHeap.minPointer = heap2.minPointer;
                    }
                }else{
                    unionHeap.minPointer = heap2.minPointer;
                }
                unionHeap.listSize = heap1.listSize + heap2.listSize;
            }
            return unionHeap;
        }
        
        /**
         * @return HeapNode 
    	 * @param FibonacciHeap h1
    	 *   This method  remove Min node from Fibonacci heap 
    	 */ 
        
        public static HeapNode extractMin(FibonacciHeapNew H){
            HeapNode toReturn = H.minPointer;
            if(null!=toReturn){
                int noOfChildren = toReturn.degree;
                HeapNode x = toReturn.child;
                HeapNode temp = null; 
                // for each child x of z add to root list of main heap
                // we will do this by first removing x from children list and then add it to the root list
                while(noOfChildren > 0){
                    temp = x.right;
                    x.left.right = x.right; // two steps
                    x.right.left = x.left; // to remove from children list         
                    x.left = H.minPointer;
                    x.right = H.minPointer.right; // all these steps
                    H.minPointer.right = x; // to add x to root list
                    x.right.left =x;
                    
                    // put x parent to null
                    x.parent = null;
                    x = temp;
                    noOfChildren--;
                }
                
                // remove z from root list of main heap
                toReturn.left.right = toReturn.right;
                toReturn.right.left = toReturn.left;
                
                if(toReturn == toReturn.right){
                    H.minPointer =null;
                }else{
                    H.minPointer = toReturn.right;
                    consolidate(H);
                }
                H.listSize--;
            }
            return toReturn;
        }
          
        /**
         * @return null 
    	 * @param FibonacciHeap h1
    	 *   This method  does the pairwise combination of the nodes of same degree
    	 */ 
        
        public static void consolidate(FibonacciHeapNew H){
            int arraySize = ((int) Math.floor(Math.log(H.listSize) * oneByLogPhi)) + 1;
            List<HeapNode> array = new ArrayList<HeapNode>(arraySize);
            for (int i = 0; i < arraySize; i++) {
                array.add(null);
            }
            int noOfRootElements = 0;
            HeapNode x = H.minPointer;
            if (x != null) {
            	noOfRootElements++;
                x = x.right;
                while (x != H.minPointer) {
                	noOfRootElements++;
                    x = x.right;
                }
            }
            // step 1 : For each node in root list do
            while (noOfRootElements > 0) {
                int d = x.degree;
                HeapNode next = x.right;
                for (;;) {
                    HeapNode y = array.get(d);
                    if (y == null) {
                        break;
                    }
                    // step 2 : make one of the nodes a child of the other.
                    if (x.key > y.key) {
                        HeapNode temp = y;
                        y = x;
                        x = temp;
                    }
                    // step 3 : HeapNode y disappears from root list.
                    fibHeapLink(H,y,x);
                    array.set(d, null);
                    d++;
                }
                array.set(d, x);
                x = next;
                noOfRootElements--;
            }
            H.minPointer = null;
            for (int i = 0; i < arraySize; i++) {
                HeapNode y = array.get(i);
                if (y == null) {
                    continue;
                }
                if (H.minPointer != null) {
                    y.left.right = y.right;
                    y.right.left = y.left;
                    y.left = H.minPointer;
                    y.right = H.minPointer.right;
                    H.minPointer.right = y;
                    y.right.left = y;
                    if (y.key < H.minPointer.key) {
                    	H.minPointer = y;
                    }
                } else {
                	H.minPointer = y;
                }
            }
        }
         
        
        /**
         * @return null 
    	 * @param FibonacciHeap h1, HeapNode h1, HeapNode h2
    	 *   This method  does the makes one node a child of other by comparing their keys
    	 */ 
       
        public static void fibHeapLink(FibonacciHeapNew H,HeapNode h1, HeapNode h2){
            // step 1 : remove heap1(y) from rootlist of  heap
            h1.left.right = h1.right;
            h1.right.left = h1.left;
            // step 2 : make h2 (y) a child of h1(x)
            h1.parent = h2;
            // if h2 has a child we have to assign properly
            if(null!=h2.child){
                h1.left = h2.child;
                h1.right = h2.child.right;
                h2.child.right = h1;
                h1.right.left = h1;
            }else{
                h2.child = h1;
                h1.right = h1;
                h1.left = h1;
            }
            h2.degree++;
            // step 3 : y.mark = FALSE
            h1.isMarked = false;
        }
        
        
        /**
         * @return null 
    	 * @param FibonacciHeap h1, HeapNode x, double k
    	 *   This method  decrease the value of a key to the value given
    	 *   and then again does cut and cascading cut of the new node.
    	 */ 
        public static void decreaseKey(FibonacciHeapNew H, HeapNode x, double k){
            if(k > x.key){
                throw new IllegalArgumentException("new key is greater than current key");
            }
            x.key = k;
            HeapNode y = x.parent;
            if((null!=y)&&(x.key<y.key)){
                cut(H,x,y);
                cascadingCut(H,y);
            }
            if(x.key < H.minPointer.key){
                H.minPointer =x;
            }
        }
        
        
        /**
         * @return null 
    	 * @param FibonacciHeap h1, HeapNode x, HeapNode y
    	 *   This method  cuts the given nodes and adds it to the root list of the heap
    	 */ 
        public static void cut(FibonacciHeapNew H,HeapNode x, HeapNode y){
            // step 1 : remove x from childlist of y, decrementing y degree
            x.right.left = x.left;
            x.left.right = x.right;
            y.degree--;
            if(y.degree==0){
                y.child=null;
            }
            if(y.child==x){ // reset y child
                y.child = x.right;
            }
            // step 2 : add x to rootlist of heap
            x.left = H.minPointer;
            x.right = H.minPointer.right;
            H.minPointer.right = x;
            x.right.left = x;
            // step 3 : make x parent to null
            x.parent =null;
            // step 4 : mark x as false
            x.isMarked = false;
        }
        
        
        /**
         * @return null 
    	 * @param FibonacciHeap h1 heapnode y
    	 *   This method  does the cascading cut of a node based on the isMarked value
    	 *   which denotes if a node has lost a child after it has become a child of another node.
    	 */ 
        public static void cascadingCut(FibonacciHeapNew H,HeapNode y){
            // step1 : make z as parent of y
            HeapNode z = y.parent;
            if(null!=z){
                if(y.isMarked==false){
                    y.isMarked=true;
                }else{
                    cut(H,y, z);
                    cascadingCut(H,z);
                }
            }
        }
        
        
        /**
         * @return null 
    	 * @param FibonacciHeap h1, heapnode x
    	 *   This method  deletes a heap node from the heap
    	 */ 
        
        public void delete(FibonacciHeapNew H, HeapNode x){
            decreaseKey(H,x, Double.NEGATIVE_INFINITY);
            extractMin(H);
        }
        
    }

	
    /**
     * @return null 
	 * @param args[]
	 *   This method  is the driver of the execution of a java program. This takes 
	 *   input parameters as string and does the execution.
	 */ 
	
	public static void main(String[] args) {
		String firstParameter =  null;
		String secondParameter = null;
		String thirdParameter =  null;
		if(null!=args[0]){
			firstParameter = args[0];	
		}
		if(firstParameter.equalsIgnoreCase("-r")){
			Integer noOfVertices = 0;
			Double density = 0.0;
			if(null!=args[1]){
				secondParameter = args[1];
				 noOfVertices = Integer.parseInt(secondParameter);
			}
			if(null!=args[2]){
				thirdParameter = args[2];
				 density = Double.parseDouble(thirdParameter);
			}
			MST randomGraph = generateRandomConnectedGraph(noOfVertices, density);
			displayPrimSimple(randomGraph, 0);
			displayPrimFHeap(randomGraph, 0);
			}else if(firstParameter.equalsIgnoreCase("-s")){
			if(null!=args[1]){
				secondParameter = args[1];
			}
			MST generatedGraph = new MST(secondParameter);
			displayPrimSimple(generatedGraph, 0);
			displayPrimFHeap(generatedGraph, 0);
		}
	}
}
