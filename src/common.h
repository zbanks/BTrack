#include <math.h>
#include <stdlib.h>

#define true 1
#define false 0

#ifndef M_PI
#define M_PI 3.14159265359
#endif

#define MIN(x, y) ((x) < (y)) ? (x) : (y)
#define MAX(x, y) ((x) > (y)) ? (x) : (y)

#define ASSERT(x) if(!x){fprintf(stderr,"Error, assertion failed: " __FILE__ " line %d.\n", __LINE__); exit(EXIT_FAILURE);}
