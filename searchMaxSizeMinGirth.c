#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <getopt.h>
#include <stdbool.h>

#include "lib/bitset/bitset.h"


/************************************************
 * Local search algorithm to find graphs
 * of large size and no 3- or 4-cycles.
 *
 * Author: Tibo Van den Eede
 *
 * Some code has been taken from or is based on
 * the following repositories:
 * - GenK2: https://github.com/JarneRenders/GenK2Hypohamiltonian
 * - K2Ham: https://github.com/JarneRenders/K2-Hamiltonian-Graphs/blob/main/readGraph/readGraph6.c
 * - biregGirthGraphs: https://github.com/tiboat/biregGirthGraphs
 * - faultCost: https://github.com/JarneRenders/faultCost/
 ***********************************************/


#define USAGE \
"\nUsage:./searchMaxSizeMinGirth<bitsetSize> -n# [-r -h] [-g# -p# -R# -i# -d#]\n\n"

#define HELPTEXT \
"Searches for a graph of order n with minimum girth g and large size, potentially\n\
starting from a given graph from stdin in graph6 format. It writes graphs\n\
in graph6 format of order n, girth g and larger, increasing size to\n\
stdout that the search algorithm finds. Auxiliary information of\n\
program execution is written to stderr.\n\n\
Underneath are the required arguments.\n\
\n\
    -n#, --numVertices=#\n\
        number of vertices of the potential start graph and graphs\n\
        to search for\n\
\n\
Underneath are the optional arguments.\n\
    -d#, --delEdgeMax=#\n\
        a random number of edges in {1,...,#} is deleted after one iteration.\n\
        The default is 3.\n\
    -g#, --girth=#\n\
        minimum girth of the potential start graph and the obtained graph.\n\
        The default is 5.\n\
    -h, --help\n\
        print this help message\n\
    -i#, --iterations=#\n\
        number of iterations that the local search will perform.\n\
        The default is 80 000 iterations.\n\
    -p#, --prob=#\n\
        probability to add edge between vertices with the largest degree sum\n\
        instead of a random edge. The default probability is 0.5.\n\
    -r, --readGraph\n\
        enables to read a graph in graph6 format via stdin to start the search.\n\
        In case this argument is not set, the search starts from an empty graph.\n\
    -R#, --recencyLimit=#\n\
        do not delete edges that where deleted in the last # iterations.\n\
        The default is 30.\n"



/************************************************
 * Macros
 ***********************************************/

// Macros based on: biregGirthGraphs

// Use this macro when you are not sure which of the two vertices i or j is smaller
#define DIST(i,j) (i < j ? dist[i][j] : dist[j][i])

// Returns the minimum of two values
#define MIN(a,b) (((a)<(b))?(a):(b))
// Returns the maximum of two values
#define MAX(a,b) (((a)>(b))?(a):(b))

// Returns the maximum number of edges in a graph of order numVertices
#define MAXEDGES ((MAXVERTICES*(MAXVERTICES-1)) >> 1)


/************************************************
 * Global variables
 ***********************************************/

// a struct that represents an edge with its two vertices v1 and v2
typedef struct { int v1,v2; } edge;

int n;                 // order of the graph
int startNumEdges;     // number of edges of the starting graph
int bestNumEdges;      // most number of edges found in a graph
int numEdges = 0;      // number of edges of the current graph
int girth = 5;         // the minimum girth that the graph needs to have

int iterDeletedEdge[MAXVERTICES][MAXVERTICES];
                    // iterDeletedEdge[v1][v2] == i  if the edge {v1,v2} was last deleted in iteration i
                    // if {v1,v2} has not been deleted so far, then iterDeletedEdge[v1][v2] == 0

int dist[MAXVERTICES][MAXVERTICES]; // distance matrix of the current graph
int deg[MAXVERTICES];               // deg[v] is the degree of vertex v
edge bestEdges[MAXEDGES];           // list of edges to store the edges with the largest degree sum
edge* toDelEdges;                   // list of edges to delete

bitset possibleEdges[MAXVERTICES];
                    // if v1 is element of the set possibleEdges[v2], then the edge {v1,v2} does not exist
                    // in the current graph and {v1,v2} can be added without creating a cycle
                    // of length smaller than the girth

int numDeletableEdges; // Number of edges that are eligible for deletion (becasue they have not been deleted
                       // too recent.

// Data structures used for efficient updating of distance matrix
int recompute[2*MAXEDGES]; // queue to store pairs of vertices for which their distance should be recomputed
bitset vis;                // a set of visited vertices used for BFS
int q[MAXVERTICES];        // a queue of vertices used for BFS
int qDist[MAXVERTICES];    // a queue of distances corresponding to the queue q in BFS

// Algorithm parameters
bool readInitGraph = false; // true if an input graph in graph6 format should be read from stdin; otherwise the
                            // starting graph if an empty graph
int maxNumIters = 80000;    // number of iterations that the local search should be performed
int recencyLimit = 30;      // edges that where deleted in the last recencyLimit iterations should not be deleted again
double prob = 0.5;          // probability to choose for the largest degree sum strategy when adding an edge
int numEdgesToDelLimit = 3; // maximum number of edges that can be deleted at the end of an iteration



/************************************************
 * Random functions
 ***********************************************/

// Returns a pseudo-random double between lo and hi
double randDoubleBetween(double lo, double hi) {
    return lo + (hi - lo) * (double)rand() / RAND_MAX;
}

// Return a psuedo-random integer in {lo, ..., hi} (so both lo and hi inclusive)
int randIntBetween(int lo, int hi) {
    return lo + rand() % (hi - lo + 1);
    // rand() Returns a pseudo-random integer between 0 and RAND_MAX.
}



/************************************************
 * Printing functions
 ***********************************************/

// Based on: GenK2
// Prints the adjacency list adjList
void printAdjacencyList(int n, bitset adjList[]) {
    for(int i = 0; i < n; i++) {
        fprintf(stderr, "%d:", i);
        forEach(j, adjList[i]) {
            fprintf(stderr, " %d", j);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

// Prints the distance matrix
void printDistMatrix(int n) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            fprintf(stderr, "%d ", dist[i][j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

// Print a bitset
void printBitset(bitset set) {
    forEach(element, set) {
        fprintf(stderr, "%d ", element);
    }
    fprintf(stderr, "\n");
}

// Print the possible edges
void printPossibleEdges(int n) {
    for (int i = 0; i < n; ++i) {
        fprintf(stderr,"%d: ", i);
        printBitset(possibleEdges[i]);
    }
}

// Based on: GenK2
// Function to convert adjacency list to graph6 format and print it
void printGraph6String(int n, bitset adjList[]) {
    // Mainly copied from Gen-K2

    char graphString[8 + n*(n - 1)/2];
    int pointer = 0;

    //  Save number of vertices in the first one, four or 8 bytes.
    if(n <= 62) {
        graphString[pointer++] = (char) n + 63;
    }
    else if(n <= 258047) {
        graphString[pointer++] = 63 + 63;
        for(int i = 2; i >= 0; i--) {
            graphString[pointer++] = (char) ((n >> i*6) & 63) + 63;
        }
    }
    else if(n <= INT_MAX) {
        graphString[pointer++] = 63 + 63;
        graphString[pointer++] = 63 + 63;
        for(int i = 5; i >= 0; i--) {
            graphString[pointer++] = (char) ((n >> i*6) & 63) + 63;
        }
    }
    else {
        fprintf(stderr, "Error: number of vertices too large.\n");
        exit(1);
    }

    // Group upper triangle of adjacency matrix in groups of 6. See B. McKay's
    // graph6 format.
    int counter = 0;
    char charToPrint = 0;
    for(int i = 1; i < n; i++) {
        for(int j = 0; j < i; j++) {
            charToPrint = charToPrint << 1;
            if(contains(adjList[j], i)) {
                charToPrint |= 1;
            }
            if(++counter == 6) {
                graphString[pointer++] = charToPrint + 63;
                charToPrint = 0;
                counter = 0;
            }
        }
    }

    //  Pad final character with 0's.
    if(counter != 0) {
        while(counter < 6) {
            charToPrint = charToPrint << 1;
            if(++counter == 6) {
                graphString[pointer++] = charToPrint + 63;
            }
        }
    }

    //  End with newline and end of string character.
    graphString[pointer++] = '\n';
    graphString[pointer++] = '\0';
    printf("%s", graphString);
}



/************************************************
 * Read graph
 ***********************************************/

// Based on: K2Ham
// Return the previous bit of a char. Unsafe because no defined behaviour if
// character = 0. ctz and clz work with 32 bit numbers.
#define unsafePrev(character, current) (__builtin_ctz(character) - (current) >= 0 ? -1 : (current) -__builtin_clz((character) << (32 - current)) - 1)
#define prev(character,current) ((character) ? unsafePrev(character,current) : -1)

// Based on: K2Ham
bool loadGraph(bitset adjList[], int n, const char * graphString) {
    //	First position after the information relating to the number of vertices.
    int startIndex = 0;
    if (graphString[startIndex] == '>') { // Skip >>graph6<< header.
        startIndex += 10;
    }
    if (n <= 62) {
        startIndex += 1;
    }
    else if (n <= 258048) {
        startIndex += 4;
    }
    else {
        fprintf(stderr,"Error: Program can only handle graphs with %d vertices or fewer.\n", 258048);
        return false;
    }

    //	Taking the remaining characters, subtracting by 63 and concatenating
    //	them represents in binary the concatenation of the upper
    //	triangle of the adjacency matrix of the graph. I.e. the first
    //	bit of this concatenation represents (0,1), the second (0,2), the
    //	third (1,2), the fourth (0,3), etc.
    int currentVertex = 1;
    int sum = 0;
    char finalChar = '\n';
    for (int index = startIndex; graphString[index] != '\n' && (finalChar = graphString[index]) != '\0'; index++) {
        int i;
        for (i = prev(graphString[index] - 63, 6); i != -1; i = prev(graphString[index] - 63, i)) {
            while(5-i+(index-startIndex)*6 - sum >= 0) {
                sum += currentVertex;
                currentVertex++;
            }
            sum -= --currentVertex;
            int neighbour = 5-i+(index - startIndex)*6 - sum;
            add(adjList[currentVertex], neighbour);
			add(adjList[neighbour], currentVertex);
         }
     }
    if(finalChar == '\0') {
        fprintf(stderr,"Error: The g6 string should end with a newline character.\n");
        return false;
    }
    return true;
}

void makeOrReadStartGraph(bitset adjList[]) {
    if (readInitGraph) {
        char *graphString = NULL;
        size_t size;
        if (getline(&graphString, &size, stdin) != -1) {
            if (!loadGraph(adjList, n, graphString)) {
                exit(-1);
            }
        } else {
            fprintf(stdout, "Error reading graph6 string\n");
        }
        free(graphString);
    }
}


/************************************************
 * Graph invariant functions
 ***********************************************/


// Returns true if the graph is connected; otherwise returns false.
bool isConnected(int n, bitset adjList[]) {
    int front = 0, back = 0;
    vis = EMPTY;
    add(vis, 0);
    q[back++] = 0;
    int totVisited = 1;
    while (front < back) {
        int now = q[front++];
        forEach(neigh, adjList[now]) {
            if (!contains(vis,neigh)) {
                add(vis, neigh);
                q[back++] = neigh;
                totVisited++;
            }
        }
    }
    return totVisited == n;
}

// Returns true if there exists a possible edge (that is a pair (v1,v2) such that adding the edge
// {v1,v2} would not create a cycle of length shorter than the girth).
bool edgesCanBeAdded(int n) {
    for (int i = 0; i < n; ++i) {
        if (!isEmpty(possibleEdges[i]))
            return true;
    }
    return false;
}

// Returns true if the difference between currIter and iterDelEdge is smaller than or equal to the recencyLimit
// except if iterDelEdge == 0 or the currIter is smaller than or equal to recencyLimit
bool edgeIsTooRecent(int iterDelEdge, int currIter) {
    if (currIter <= recencyLimit)
        return false;
    return iterDelEdge > 0 && (currIter - iterDelEdge) <= recencyLimit;
}


/************************************************
 * Data structure update functions
 ***********************************************/

// Compute the distance matrix of adjList from scratch
void computeDists(bitset adjList[], int n) {
    bool* visited = (bool*)calloc(n, sizeof(bool));
    if (visited == NULL) {
        fprintf(stderr,"Error allocating memory for visited.\n");
        exit(-1);
    }
    int* queue = (int*)malloc(n * sizeof(int));
    if (queue == NULL) {
        fprintf(stderr,"Error allocating memory for queue.\n");
        exit(-1);
    }
    for (int src = 0; src < n; ++src) {
        for (int i = 0; i < n; ++i) {
            visited[i] = false;
        }
        int front = 0, rear = 0; // Queue for BFS

        // Initialize distances
        for (int i = 0; i < n; i++) dist[src][i] = n;
        dist[src][src] = 0; // Distance to itself is 0

        // Start BFS
        queue[rear++] = src;
        visited[src] = true;

        while (front < rear) {
            int u = queue[front++];
            for (int v = 0; v < n; v++) {
                // If there's an edge and v is not visited
                if (contains(adjList[u], v) && !visited[v]) {
                    dist[src][v] = dist[src][u] + 1; // Update distance
                    visited[v] = true;    // Mark as visited
                    queue[rear++] = v;   // Enqueue
                }
            }
        }
    }
    free(visited);
    free(queue);
}

// Update the distance matrix after numDeletedEdges stored in toDelEdges have been deleted.
void updateDistsAfterDelEdges(bitset adjList[], edge* toDelEdges, int numDeletedEdges) {

    if (!isConnected(n, adjList)) {
        computeDists(adjList, n);
        return;
    }

    int numRecompute = 0;

    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            for (int k = 0; k < numDeletedEdges; k++) {
                if (dist[i][j] == dist[i][toDelEdges[k].v1] + 1 + dist[toDelEdges[k].v2][j] ||
                    dist[i][j] == dist[i][toDelEdges[k].v2] + 1 + dist[toDelEdges[k].v1][j]) {
                    recompute[numRecompute++] = i;
                    recompute[numRecompute++] = j;
                    break;
                }
            }
        }
    }

    int start, end;

    for (int i = 0; i < numRecompute; i+=2) {
        start = recompute[i];
        end = recompute[i+1];

        // Start BFS
        vis = EMPTY;
        int front = 0, rear = 0; // Queue for BFS
        q[rear] = start;
        qDist[rear] = 0;
        rear++;
        add(vis, start);
        bool done = false;
        while (front < rear && !done) {
            int u = q[front];
            int vDist = qDist[front]+1;
            front++;
            forEach(v, adjList[u]) {
                if (v == end) {
                    dist[start][end] = vDist;
                    dist[end][start] = vDist;
                    done = true;
                    break;
                }
                if (!contains(vis,v)) {
                    add(vis, v);
                    q[rear] = v;   // Enqueue
                    qDist[rear] = vDist;
                    rear++;
                }
            }
        }
        if (!done) {
            fprintf(stderr,"Error: BFS not found start:%d to end:%d, but graph is connected.\n", start, end);
            exit(-1);
        }
    }
}

// Compute the legal edges from scratch
void computeLegalEdges(int n) {
    for (int i = 0; i < n; ++i) {
        possibleEdges[i] = EMPTY;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            if (dist[i][j] >= girth-1) {
                add(possibleEdges[i], j);
                add(possibleEdges[j], i);
            }
        }
    }
}



// Based on: biregGirthGraphs
// Add the edge {v1,v2} at iteration iter to the graph represented by adjList[]
void addEdge(bitset adjList[], int n, int v1, int v2, int iter) {

  	if (iter != 0) {
        if (contains(adjList[v1], v2)) {
        	fprintf(stderr,"Error trying to add edge %d %d, but already exists.\n", v1, v2);
        	exit(-1);
    	}
    	add(adjList[v1], v2);
    	add(adjList[v2], v1);
    }

    numEdges++;
    deg[v1]++;
    deg[v2]++;
    dist[v1][v2] = 1;
    dist[v2][v1] = 1;
    removeElement(possibleEdges[v1], v2);
    removeElement(possibleEdges[v2], v1);

    int potDist1, potDist2, minPotDist;
    for(int vi1 = 0; vi1 < n; vi1++) {
        for(int vi2 = vi1+1; vi2 < n; vi2++) {
            potDist1 = DIST(vi1,v1) + 1 + DIST(vi2,v2);
            potDist2 = DIST(vi1,v2) + 1 + DIST(vi2,v1);
            minPotDist = MIN(potDist1,potDist2);
            if(dist[vi1][vi2] > minPotDist) {
                // if previous distance was larger than or equal to girth-1
                // and now it's smaller than girth-1, then we can decrement
                // the amount of legal choices for both vertices
                if(dist[vi1][vi2] >= girth - 1 && minPotDist < girth - 1) {
                    removeElement(possibleEdges[vi1], vi2);
                    removeElement(possibleEdges[vi2], vi1);
                }
                dist[vi1][vi2] = minPotDist;
                dist[vi2][vi1] = minPotDist;
            }
        }
    }
}

// Delete the edge {v1,v2} at iteration currIter. Doesn't update the distance matrix.
void delEdge(bitset adjList[], int v1, int v2, int currIter) {
    removeElement(adjList[v1], v2);
    removeElement(adjList[v2], v1);
    iterDeletedEdge[v1][v2] = currIter;
    iterDeletedEdge[v2][v1] = currIter;
    numEdges--;
    deg[v1]--;
    deg[v2]--;
}

// Sets the variable numDeletableEdges to how much edges can be deleted
// satisfying that they have not been deleted too recent.
void setNumDeletableEdges(bitset adjList[], int n, int currIter) {
    numDeletableEdges = 0;
    for(int i = 0; i < n; i++) {
        forEachAfterIndex(j, adjList[i], i) {
            if (!edgeIsTooRecent(iterDeletedEdge[i][j], currIter)) {
                numDeletableEdges++;
            }
        }
    }
}


/************************************************
 * Allocate (init) and free (end) data structures
 ***********************************************/


bitset* init() {
	srand(time(NULL));

    bitset* adjList = (bitset *)malloc(n * sizeof(bitset));
    if (adjList == NULL) {
        fprintf(stderr,"Error allocating memory for adjList.\n");
        exit(-1);
    }
    toDelEdges = (edge *)malloc(numEdgesToDelLimit * sizeof(edge));
    if (toDelEdges == NULL) {
        fprintf(stderr,"Error allocating memory for toDelEdges.\n");
        exit(-1);
    }

    for (int i = 0; i < n; ++i) {
        possibleEdges[i] = EMPTY;
        adjList[i] = EMPTY;
        for (int j = i+1; j < n; ++j) {
            dist[i][j] = n;
            dist[j][i] = n;
            add(possibleEdges[i], j);
            add(possibleEdges[j], i);
        }
    }

    if (readInitGraph) {
		makeOrReadStartGraph(adjList);
		for (int i = 0; i < n; ++i) {
    		forEachAfterIndex(j, adjList[i], i) {
                addEdge(adjList, n, i, j, 0);
        	}
		}
    }

    fprintf(stderr,"startNumEdges = %d\n", numEdges);

    bestNumEdges = numEdges;

    return adjList;
}

void end(bitset* adjList) {
    free(adjList);
    free(toDelEdges);
}



/************************************************
 * Local search functions
 ***********************************************/


// Choose and add an edge.
void addBestEdge(bitset adjList[], int n, int iter) {
    if (randDoubleBetween(0, 1) < prob) {
        int largestDegSum = 0;
        int numBestEdges = 0;
        // Find largest degree sum
        for (int v1 = 0; v1 < n; ++v1) {
            forEachAfterIndex(v2, possibleEdges[v1], v1) {
                int degreeSum = deg[v1]+deg[v2];
                if (degreeSum == largestDegSum) {
                    numBestEdges++;
                    bestEdges[numBestEdges-1] = (edge) { v1,v2 };
                } else if (degreeSum > largestDegSum) {
                    numBestEdges = 1;
                    bestEdges[numBestEdges-1] = (edge) { v1,v2 };
                    largestDegSum = degreeSum;
                }
            }
        }

        if (numBestEdges != 0) {
            // Add a random edge if there is no "best" edge
            int randNum = randIntBetween(0, numBestEdges-1);
            addEdge(adjList, n, bestEdges[randNum].v1, bestEdges[randNum].v2, iter);
            return;
        }
    }

    int randV1, randV2;
    do {
        randV1 = randIntBetween(0,n-1);
        randV2 = randIntBetween(0,n-1);
    } while (randV1 == randV2 ||
             contains(adjList[randV1], randV2) ||
             !contains(possibleEdges[randV1], randV2));
    addEdge(adjList, n, randV1, randV2, iter);
}

// Choose and delete edges
void deleteEdges(bitset adjList[], int n, int currIter) {
    setNumDeletableEdges(adjList, n, currIter);

    int numToDelEdges = randIntBetween(1, numEdgesToDelLimit);

    if (numDeletableEdges == 0) {
        fprintf(stderr,"There are not no deletable edges. This is caused by setting "
                             "the delEdgeMax and recencyLimit too high.\n");
        fprintf(stderr,"Abort\n");
        // Clean up allocations and exit.
        end(adjList);
        exit(0);
    }
    if (numToDelEdges > numDeletableEdges) {
        numToDelEdges = numDeletableEdges;
    }
    int numDeletedEdges = 0;
    int randV1, randV2;

    // Kind of arbitrary loop limit, but this will probably only break
    // when there is an infinite loop.
    while (numDeletedEdges < numToDelEdges) {
        do {
            randV1 = randIntBetween(0,n-1);
            randV2 = randIntBetween(0,n-1);
        } while (randV1 == randV2 ||
                 !contains(adjList[randV1], randV2) ||
                 edgeIsTooRecent(iterDeletedEdge[randV1][randV2], currIter));
        delEdge(adjList, randV1, randV2, currIter);
        toDelEdges[numDeletedEdges] = (edge) { .v1 = randV1, .v2 = randV2 };
        numDeletedEdges++;
    }
    updateDistsAfterDelEdges(adjList, toDelEdges, numDeletedEdges);
    computeLegalEdges(n);
}



// Performs the local search algorithm
void doLocalSearch() {
    bitset* adjList = init();
    int numIters = 0;

    // Print start graph
    printGraph6String(n, adjList);

    while (true) {
        numIters++;
        if (numIters >= maxNumIters) break;
        while (edgesCanBeAdded(n)) {
            addBestEdge(adjList, n, numIters);
        }
        if (numEdges > bestNumEdges) {
            bestNumEdges = numEdges;
            // fprintf(stderr,"|E|: %d\n", numEdges);
            // printGraph6String(n, adjList);
        }
        printGraph6String(n, adjList);
        deleteEdges(adjList, n, numIters);
    }
    end(adjList);
}


/************************************************
 * Main function
 ***********************************************/

int main(int argc, char* argv[]) {

    // Parsing arguments based on: faultCost

    int opt;
    while (true) {
        int option_index = 0;
        static struct option long_options[] =
        {
            {"prob", required_argument, NULL, 'p'},
            {"numVertices", required_argument, NULL, 'n'},
            {"readGraph", no_argument, NULL, 'r'},
            {"girth", required_argument, NULL, 'g'},
            {"recencyLimit", required_argument, NULL, 'R'},
            {"delEdgeMax", required_argument, NULL, 'd'},
            {"iterations", required_argument, NULL, 'i'},
            {"help", no_argument, NULL, 'h'}
        };

        opt = getopt_long(argc, argv, "d:g:hi:n:p:rR:", long_options, &option_index);
        if (opt == -1) break;
        switch(opt) {
            case 'p':
                prob = (double) strtof(optarg, NULL);
                break;
            case 'n':
                n = (int) strtol(optarg, NULL, 10);
                break;
            case 'g':
                girth = (int) strtol(optarg, NULL, 10);
                break;
            case 'r':
                readInitGraph = true;
                break;
            case 'R':
                recencyLimit = (int) strtol(optarg, NULL, 10);
                break;
            case 'd':
                numEdgesToDelLimit = (int) strtol(optarg, NULL, 10);
                break;
            case 'i':
                maxNumIters = (int) strtol(optarg, NULL, 10);
                break;
            case 'h':
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr, "%s", HELPTEXT);
                return 0;
            case '?':
                fprintf(stderr,"Error: Unknown option: %c\n", optopt);
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr,"Use %s --help for more detailed instructions.\n", argv[0]);
                return 1;
        }
    }

    // Check that required arguments are set
    if (n == -1) {
        fprintf(stderr,"Error: Missing required argument.\n");
        fprintf(stderr, "%s\n", USAGE);
        exit(-1);
    }

    fprintf(stderr,"Local search starting with:\n");
    fprintf(stderr,"read start graph = %s\n", readInitGraph ? "true" : "false");
    fprintf(stderr,"girth = %d\n", girth);
    fprintf(stderr,"probability = %f\n", prob);
    fprintf(stderr,"recency limit = %d\n", recencyLimit);
    fprintf(stderr,"max edge to del = %d\n", numEdgesToDelLimit);
    fprintf(stderr,"number of iterations = %d\n", maxNumIters);

    clock_t start = clock();

    doLocalSearch();

    clock_t end = clock();

    fprintf(stderr,"Local search finished\n");
    fprintf(stderr,"bestNumEdges = %d\n", bestNumEdges);
    fprintf(stderr,"time = %fs\n", (double)(end - start) / CLOCKS_PER_SEC);

    return 0;
}
