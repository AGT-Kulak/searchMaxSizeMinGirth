
// Based on: https://github.com/JorikJooken/vertexGirthRegularGraphs/blob/master/Code/bitset.h

#ifndef BITSETCHOOSER
#define BITSETCHOOSER

#ifdef USE_64_BIT
	#include "bitset64Vertices.h"
	#define MAXVERTICES 64

#elif defined(USE_128_BIT)
	#include "bitset128VerticesArray.h"
	#define MAXVERTICES 128

#elif defined(USE_192_BIT)
	#include "bitset192VerticesArray.h"
	#define MAXVERTICES 192

#elif defined(USE_256_BIT)
	#include "bitset256VerticesArray.h"
	#define MAXVERTICES 256

#endif

#endif
