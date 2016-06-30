#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

typedef float btrack_chunk_t;

#ifndef M_PI
#define M_PI 3.14159265359
#endif

#define BTRACK_MIN(x, y) ((x) < (y)) ? (x) : (y)
#define BTRACK_MAX(x, y) ((x) > (y)) ? (x) : (y)

#define BTRACK_STRINGIFY(x) BTRACK_STRINGIFY2(x)
#define BTRACK_STRINGIFY2(x) #x
#define BTRACK_ASSERT(x) if(!(x)){fprintf(stderr,"Error, assertion failed: " __FILE__ " line %d: " BTRACK_STRINGIFY(x) "\n", __LINE__); exit(EXIT_FAILURE);}
