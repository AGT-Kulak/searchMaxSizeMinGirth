
// Downloaded from: https://github.com/JarneRenders/K2-Hamiltonian-Graphs/blob/main/bitset128VerticesArray.h

// FOR GRAPHS UP TO 128 VERTICES
#ifndef BITSET_MACROS
#define BITSET_MACROS

#include <stdint.h>

//  Bitset macros, we assume nodes are labeled 0,1,2,...
typedef struct bitset {uint64_t parts[2];} bitset;

//  Returns an empty bitset. 
#define EMPTY (bitset) {(uint64_t) 0LL, (uint64_t) 0LL}

//  Returns a bitset containing only node.
#define singleton(node) ((node)/64) == 0 ? (bitset) {(uint64_t) 1LL << (node), (uint64_t) 0LL} : (bitset) {(uint64_t) 0LL, (uint64_t) 1LL << (node)%64}

//  Returns the union of set1 and set2.
#define union(set1, set2) (bitset) {((set1).parts[0] | (set2).parts[0]), ((set1).parts[1] | (set2).parts[1])}

//  Returns the intersection of set1 and set2.
#define intersection(set1, set2) (bitset) {((set1).parts[0] & (set2).parts[0]), ((set1).parts[1] & (set2).parts[1])}

//  Adds node to set.
#define add(set, node) ((set) = union((set),singleton(node)))

//  Returns set1\set2 (set difference).
#define difference(set1, set2) (bitset) {((set1).parts[0] & ~(set2).parts[0]), ((set1).parts[1] & ~(set2).parts[1])}


//  Removes node from set.
#define removeElement(set, node) ((set) = difference((set), singleton(node)))

//  Check if set is empty.
#define isEmpty(set) !((set).parts[0] | (set).parts[1]) 

//  Returns the size of the set.
#define size(set) (__builtin_popcountll((set).parts[0])+__builtin_popcountll((set).parts[1]))

//	Check is set1 equals set2.
#define equals(set1, set2) (((set1).parts[0] == (set2).parts[0]) && ((set1).parts[1] == (set2).parts[1]))

#define unsafeNextLessThan63(set, current) (__builtin_clzll((set).parts[0]) + (current) >= 63 ? safeNextGreaterThanOrEq63(set, 63) : (current) + __builtin_ctzll((set).parts[0] >> ((current) + 1)) + 1)
#define unsafeNextGreaterThanOrEq63(set, current) (__builtin_clzll((set).parts[1]) + ((current) - 64) >= 63 ? -1 : (current) + __builtin_ctzll((set).parts[1] >> ((current) - 64 + 1)) + 1)
#define safeNextGreaterThanOrEq63(set, current) ((set).parts[1] ? unsafeNextGreaterThanOrEq63(set, current) : -1)
#define safeNextLessThan63(set, current) ((set).parts[0] ? unsafeNextLessThan63(set, current) : safeNextGreaterThanOrEq63(set, 63))
#define next(set, current)  ((current) < 63 ? safeNextLessThan63((set), (current)) : safeNextGreaterThanOrEq63((set), (current)))

//  Checks whether node is an element of set.
#define contains(set, node) (!isEmpty(intersection((set), singleton(node))))

//	Loops over all elements of the set. 
#define forEach(element, set) for (int element = next((set), -1); (element) != -1; (element) = next((set), (element)))

//	Loops over all elements of the set starting from start (not included). 
#define forEachAfterIndex(element, set, start) for (int element = next((set), (start)); (element) != -1; (element) = next((set), (element)))

//	Take the complement of a set in a universe consisting of sizeOfUniverse elements. 
//	E.g.: complement of {0,2} if there are only 4 elements is {1,3} and not 64-bit complement. 
#define complement(set, sizeOfUniverse) (sizeOfUniverse <= 64 ? (bitset) {~(set).parts[0] << (64-(sizeOfUniverse)) >> (64-(sizeOfUniverse)), (set).parts[1]} : (bitset) {~(set).parts[0], ~(set).parts[1] << (64 - ((sizeOfUniverse) - 64)) >> (64-((sizeOfUniverse) - 64) )})

#endif