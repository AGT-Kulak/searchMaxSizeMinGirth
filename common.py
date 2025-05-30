import networkx as nx
import os

def getNumLinesInFile(file_path):
    count = 0
    with open(file_path) as fp:
        for line in fp:
            if line.strip():
                count += 1
    return count

def getMostEdges(file_paths):
    maxEdges = -1
    for file_path in file_paths:
        Gs = nx.read_graph6(file_path)

        if not isinstance(Gs, list):
            Gs = [Gs]

        for G in Gs:
            numEdges = nx.number_of_edges(G)
            if numEdges > maxEdges:
                maxEdges = numEdges

    return maxEdges


def getBitsetSize(n):
    bitsetSize = -1
    if n <= 64:
        bitsetSize = 64
    elif n <= 128:
        bitsetSize = 128
    elif n <= 192:
        bitsetSize = 192
    elif n <= 256:
        bitsetSize = 256
    else:
        print(f"Error: n={n} is too large, only support up to n=256.")
        exit(-1)
    return bitsetSize




startGraphs = dict()

def readStartGraphs(dir_path):
    for filename in os.listdir(dir_path):
        f = os.path.join(dir_path, filename)

        readGraphs = nx.read_graph6(f)
        G = None
        if isinstance(readGraphs, list):
            G = readGraphs[0]
        else:
            G = readGraphs
        
        startGraphs[nx.number_of_nodes(G)] = f
    print(startGraphs)

def isStartGraphOfOrder(order):
    return order in startGraphs.keys()

def startGraphFileOfOrder(order):
    return startGraphs[order]