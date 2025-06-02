# searchMaxSizeMinGirth

This repository contains code and data related to the manuscript "J. Goedgebeur, J. Jooken, G. Joret and T. Van den Eede. Improved lower bounds on the maximum size of graphs with
girth 5, manuscript". First, there is a local search method, written in C,
that finds graphs of a specific order $n$, girth at least $5$ and large size. Secondly we have two Python scripts (`upRun.py` and `downRun.py`) that
perform the local search on the range of orders $\{50, \ldots, 203\}$, using graphs of one order less or more as start graph.

## Data

### Start graphs

The directory `bestGraphs/startGraphs` contains the graphs which are used as a start graph for the local search.

- All graphs up to order 64 are the graphs found by Narjess Afzaly and Brendan McKay at https://users.cecs.anu.edu.au/~bdm/data/extremal.html, more specifically downloadable in sparse6 format via [this link](https://users.cecs.anu.edu.au/~bdm/data/extremal/c34_all.tar.gz).
- The graph of order 80 was downloaded from House of Graphs via [this link](https://houseofgraphs.org/graphs/30635).
- The graphs of order 96 and 156 by Leif K. Jørgensen come from [this link](https://people.math.aau.dk/%7Eleif/research/girth5/156) and [this link](https://people.math.aau.dk/%7Eleif/research/girth5/96.html), respectively.
- The graphs of order 124, 126, 154 and 203 were made available by Geoffrey Exoo and are available via [this link](https://web.archive.org/web/20140306232202/http://ginger.indstate.edu:80/ge/CAGES/g10.5.124a), [this link](https://web.archive.org/web/20100726172606/http://ginger.indstate.edu/ge/CAGES/g10.5.126), [this link](https://web.archive.org/web/20150116023828/http://ginger.indstate.edu:80/ge/CAGES/g11.5.154) and [this link](https://web.archive.org/web/20150116023828/http://ginger.indstate.edu:80/ge/CAGES/g12.5.203), respectively 


### Results
We performed two iterations of the up and down run each, in which the top 150 graphs are kept for each order.
The directory `resultGraphs` contains for each iteration of a down and an up run, the graphs obtained for each order. More precisely, we have the zipped folders
- `resultGraphs/iter1Up_150.zip`,
- `resultGraphs/iter1Down_150.zip`, 
- `resultGraphs/iter2Up_150.zip`, and
- `resultGraphs/iter2Down_150.zip`.

When unzipped, all contain graphs of large size obtained after performing an iteration of the up or down run, saved per order:
`n<n>.g6` contains graphs (in graph6 format) of large size, girth at least 5 and order `<n>`.

The graphs in these zipped files which improve upper bounds on $n(\lbrace r, m\rbrace ; 5)$ (= the minimum order of a graph of degree set
$$\lbrace r,m\rbrace$$, with $$r < m$$, and girth $$5$$) are additionally stored in `biregular.zip`. When unzipped, the 
file `r<r>_m<m>_g5_n<n>.g6` stores the found graphs with degree set \{ `<r>`,`<m>` \}, girth 5 and order `<n>`.


## Code


The core local search algorithm is implemented in `searchMaxSizeMinGirth.c`. This algorithm searches for graphs of a given order, given minimum girth and large size.
On top of that, we have two Python scripts (`upRun.py` and `downRun.py`) that
perform the local search on the range of orders $\{50, \ldots, 203\}$, using graphs of one order less or more as start graph.

These code files make use of other code available at
- `lib/nauty` contains some files from the [nauty](https://pallini.di.uniroma1.it/) library version 2.8.9,
- `lib/bitset` contains a bitset implementation, taken over from [this repo](https://github.com/JarneRenders/K2-Hamiltonian-Graphs/blob/main/) and [this repo](https://github.com/JorikJooken/vertexGirthRegularGraphs/tree/master/Code),
- `lib/topLMostEdges.cpp` contains code that returns the top L graphs of largest size, and
- `common.py` contains some common code used by both `upRun.py` and `downRun.py`.


### searchMaxSizeMinGirth.c

#### Compilation

Only `searchMaxSizeMinGirth.c` needs to be compiled, which can be done with the following command.

```
make
```
This creates the executables `searchMaxSizeMinGirth64`, `searchMaxSizeMinGirth128`, `searchMaxSizeMinGirth192` and `searchMaxSizeMinGirth256`, 
where the ending number denotes the size of the bitset used.
The executable `searchMaxSizeMinGirth<x>` supports the local search on graphs up to and including `<x>` vertices. It is desired
to choose the executable that has the smallest `<x>` that is greater than or equal to the order of the graph, for optimal 
efficiency.

#### Execution

Usage: `./searchMaxSizeMinGirth<x> -n# [-r -h] [-g# -p# -R# -i# -d#]` with `<x>` the bitset size.

Searches for a graph of order n with minimum girth g and large size, potentially
starting from a given graph from stdin in graph6 format. It writes graphs
in graph6 format of order n, girth g and larger, increasing size to
stdout that the search algorithm finds. Auxiliary information of
program execution is written to stderr.

Underneath are the required arguments.

    -n#, --numVertices=#
        number of vertices of the potential start graph and graphs
        to search for

Underneath are the optional arguments.

    -d#, --delEdgeMax=#
        a random number of edges in {1,...,#} is deleted after one iteration.
        The default is 3.
    -g#, --girth=#
        minimum girth of the potential start graph and the obtained graph.
        The default is 5.
    -h, --help
        print this help message
    -i#, --iterations=#
        number of iterations that the local search will perform.
        The default is 80 000 iterations.
    -p#, --prob=#
        probability to add edge between vertices with the largest degree sum
        instead of a random edge. The default probability is 0.5.
    -r, --readGraph
        enables to read a graph in graph6 format via stdin to start the search.
        In case this argument is not set, the search starts from an empty graph.
    -R#, --recencyLimit=#
        do not delete edges that where deleted in the last # iterations.
        The default is 30.

### upRun.py and donwnRun.py

#### Compilation

Compile `searchMaxSizeMinGirth.c` (see above).

Compile `lib/topLMostEdges.cpp` with the following commands
```
cd lib
g++ -std=c++11 -O3 topLMostEdges.cpp -o topLMostEdges
```

#### Execution

You should always first perform an up run, before executing a down run.

```
python upRun.py <iter> <l> <max_processes>
```
```
python downRun.py <iter> <l> <max_processes>
```

- `<iter>`: iteration number (i.e., 1, 2, ...)
- `<l>`: maximum number of graphs to save for each order.
- `<max_processes>`: maximum number of executions of `searchMaxSizeMinGirth.c` to perform simultaneously. This can be useful in order to not overload the computer system you are running this program on.

The top `<l>` graphs obtained for order `<n>` are saved in `bestGraph/iter<iter>Up_<l>/n<n>.g6` for an up run and in `bestGraph/iter<iter>Down_<l>/n<n>.g6` for a down run.
Also, immediate graphs will be saved. More specifically, an up run will save graphs in the following created directories
- `upRun_iter<iter>_<l>/ptg/n<n>.g6`: the graphs obtained after adding an isolated vertex to each graph of order `n-1` from the same directory (so `upRun_iter<iter>_<l>/ptg/n<n-1>.g6`)
- `upRun_iter<iter>_<l>/ls/out_ls_<n>_lineNum<i>_p<p>_recency<recency>_numDel<numDelEdge>_iters<numIters>.g6`: the graphs obtained from running `searchMaxSizeMinGirth.c` on the `<i>`th graph from `upRun_iter<iter>_<l>/ptg/n<n>.g6` where `<p>`, `<recency>`, `<numDelEdge>`, `<numIters>` are the values of the parameters to `searchMaxSizeMinGirth.c`  

Likewise, a down run will save graphs in the following directories
- `downRun_iter<iter>_<l>/ptg/n<n>.g6`: the top `<l>` graphs of largest size obtained after removing a vertex of each graph of order `n+1` from the same directory (so `upRun_iter<iter>_<l>/ptg/n<n+1>.g6`)
- `downRun_iter<iter>_<l>/ls/out_ls_<n>_lineNum<i>_p<p>_recency<recency>_numDel<numDelEdge>_iters<numIters>.g6`: the graphs obtained from running `searchMaxSizeMinGirth.c` on the `<i>`th graph from `downRun_iter<iter>_<l>/ptg/n<n>.g6` where `<p>`, `<recency>`, `<numDelEdge>`, `<numIters>` are the values of the parameters to `searchMaxSizeMinGirth.c`

#### Example 
When performing two iterations (of an up and down run each) while keeping the top 100 graphs of largest size and running 
at most 50 instances of the search simultaneously, one can run the following commands (in the given order).
```
python upRun.py 1 100 50
python downRun.py 1 100 50
python upRun.py 2 100 50
python downRun.py 2 100 50
```

## Authors

- Jorik Jooken, jorik [dot] jooken [at] kuleuven [dot] be

- Tibo Van den Eede, tibo [dot] vandeneede [at] kuleuven [dot] be
