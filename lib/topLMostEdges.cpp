#include <bits/stdc++.h>


/************************************************
 * This program reads graph6 strings via stdin
 * and returns the top L graphs (as graph6 strings)
 * of largest size via stdout.
 * L is given as an argument to the executable;
 * e.g. ./topLMostEdges 10
 * returns the top 10 graphs of largest size.
 *
 * Author: Jorik Jooken
 ***********************************************/

// Unsafe because no defined behaviour if character = 0. ctz and clz work with 32 bit numbers.
#define unsafePrev(character, current) (__builtin_ctz(character) - current >= 0 ? -1 : current -__builtin_clz((character) << (32 - current)) - 1)

#define prev(character,current) (character ? unsafePrev(character,current) : -1)

using namespace std;

const int nb_bits=4096;
#define MAXVERTICES 4096


int n;
vector< vector<int> > graph;


int getNumberOfVertices(string graphString) 
{
	if(graphString.size() == 0){
        printf("Error: String is empty.\n");
        abort();
    }
    else if((graphString[0] < 63 || graphString[0] > 126) && graphString[0] != '>') {
    	printf("Error: Invalid start of graphstring.\n");
    	abort();
    }

	int index = 0;
	if (graphString[index] == '>') { // Skip >>graph6<< header.
		index += 10;
	}

	if(graphString[index] < 126) { // 0 <= n <= 62
		return (int) graphString[index] - 63;
	}

	else if(graphString[++index] < 126) { 
		int number = 0;
		for(int i = 2; i >= 0; i--) {
			number |= (graphString[index++] - 63) << i*6;
		}
		return number;
	}

	else if (graphString[++index] < 126) {
		int number = 0;
		for (int i = 5; i >= 0; i--) {
			number |= (graphString[index++] - 63) << i*6;
		}
		return number;
	}

	else {
		printf("Error: Format only works for graphs up to 68719476735 vertices.\n");
		abort();
	}
}

void loadGraph(string graphString, int numberOfVertices) {
    vector<int> emp;
    graph.assign(n,emp);
	int startIndex = 0;
	if (graphString[startIndex] == '>') { // Skip >>graph6<< header.
		startIndex += 10;
	}
	if (numberOfVertices <= 62) {
		startIndex += 1;
	}
	else if (numberOfVertices <= MAXVERTICES) {
		startIndex += 4;
	}
	else {
		printf("Error: Program can only handle graphs with %d vertices or fewer.\n",MAXVERTICES);
		abort();
	}

	int currentVertex = 1;
	int sum = 0; 
	for (int index = startIndex; index<graphString.size(); index++) {
		int i;
		for (i = prev(graphString[index] - 63, 6); i != -1; i = prev(graphString[index] - 63, i)) {
			while(5-i+(index-startIndex)*6 - sum >= 0) {
				sum += currentVertex;
				currentVertex++;
			}
			sum -= --currentVertex;
			int neighbour = 5-i+(index - startIndex)*6 - sum;
            graph[currentVertex].push_back(neighbour);
            graph[neighbour].push_back(currentVertex);
		}
	}
}

int main(int argc, char *argv[])
{
    ios::sync_with_stdio(false);
    cin.tie(0);
    int K=atoi(argv[1]);
    long long nb_graphs_read_from_input=0;
    string line;
    priority_queue< pair<int, string>, vector< pair<int, string> >, greater< pair<int, string> > > pq;
    while(getline(cin,line))
    {
        //line+="\n";
        nb_graphs_read_from_input++;
        n = getNumberOfVertices(line);
        loadGraph(line,n);
        int nbEdges=0;
        for(int i=0; i<n; i++)
        {
            nbEdges+=graph[i].size();
        }
        nbEdges/=2;
        pq.push(make_pair(nbEdges,line));
        if(pq.size()>K) pq.pop();
    }
    while(!pq.empty())
    {
        pair<int, string> pa=pq.top();
        pq.pop();
        cout << pa.second << endl;
    }
    return 0;
}
