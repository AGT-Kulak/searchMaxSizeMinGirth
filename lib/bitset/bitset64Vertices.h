
// Downloaded from: https://github.com/JarneRenders/K2-Hamiltonian-Graphs/blob/main/bitset64Vertices.h

// FOR GRAPHS UP TO 64 VERTICES
#ifndef BITSET_MACROS
#define BITSET_MACROS

#include <stdint.h>

//  Bitset macros, we assume nodes are labeled 0,1,2,...
typedef uint64_t bitset;

//  Returns an empty bitset. 
#define EMPTY (bitset) 0LL

//  Returns a bitset containing only node.
#define singleton(node) ((bitset) 1LL << (node))

//  Returns the union of set1 and set2.
#define union(set1, set2) ((set1) | (set2))

//  Returns the intersection of set1 and set2.
#define intersection(set1, set2) ((set1) & (set2))

//  Adds node to set.
#define add(set, node) ((set) = union((set),singleton(node)))

//  Returns set1\set2 (set difference).
#define difference(set1, set2) ((set1) & ~(set2))

//  Removes node from set.
#define removeElement(set, node) ((set) = difference((set), singleton(node)))

//  Check if set is empty.
#define isEmpty(set) ((set) == 0)

//  Returns the size of the set.
#define size(set) (__builtin_popcountll(set))

//	Check if set1 equals set2.
#define equals(set1, set2) ((set1) == (set2))

//	Loops over all elements of the set. 
#define forEach(element, set) for (int element = next((set), -1); (element) != -1; (element) = next((set), (element)))

//	Loops over all elements of the set starting from start (not included). 
#define forEachAfterIndex(element, set, start) for (int element = next((set), (start)); (element) != -1; (element) = next((set), (element)))

//  Returns -1 if the set is empty. Otherwise it executes unsafeNext(set, current).
#define next(set, current)  (isEmpty(set) ? -1 : unsafeNext((set), (current)) ) //the builtin clz and ctz compiler functions have unexpected behavior at zero.

// Returns the first element after current in the set. If there are no more elements after current, it returns -1. 
// Unsafe because the behaviour when the set is empty is not guaranteed.
#define unsafeNext(set, current) (__builtin_clzll(set) + (current) >= 63 ? -1 : (current) + __builtin_ctzll((set) >> ((current) + 1)) + 1)

//  Checks whether node is an element of set.
#define contains(set, node) (!isEmpty(intersection((set), singleton(node))))

//	Take the complement of a set in a universe consisting of sizeOfUniverse elements. 
//	E.g.: complement of {0,2} if there are only 4 elements is {1,3} and not 64-bit complement. 
#define complement(set, sizeOfUniverse) (~(set) << (64-(sizeOfUniverse)) >> (64-(sizeOfUniverse)))

#endif