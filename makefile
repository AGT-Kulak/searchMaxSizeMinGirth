compiler=gcc
flags= -O4 -march=native

all: searchMaxSizeMinGirth64 searchMaxSizeMinGirth128 searchMaxSizeMinGirth192 searchMaxSizeMinGirth256


searchMaxSizeMinGirth64: searchMaxSizeMinGirth.c lib/bitset
	$(compiler) -DUSE_64_BIT  -o searchMaxSizeMinGirth64  searchMaxSizeMinGirth.c $(flags)

searchMaxSizeMinGirth128: searchMaxSizeMinGirth.c lib/bitset
	$(compiler) -DUSE_128_BIT -o searchMaxSizeMinGirth128 searchMaxSizeMinGirth.c $(flags)

searchMaxSizeMinGirth192: searchMaxSizeMinGirth.c lib/bitset
	$(compiler) -DUSE_192_BIT -o searchMaxSizeMinGirth192 searchMaxSizeMinGirth.c $(flags)

searchMaxSizeMinGirth256: searchMaxSizeMinGirth.c lib/bitset
	$(compiler) -DUSE_256_BIT -o searchMaxSizeMinGirth256 searchMaxSizeMinGirth.c $(flags)


clean:
	rm -f searchMaxSizeMinGirth64 searchMaxSizeMinGirth128 searchMaxSizeMinGirth192 searchMaxSizeMinGirth256