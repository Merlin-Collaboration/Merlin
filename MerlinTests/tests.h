#undef NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>

#define assert_close(x1, x2, d) assert( fabs(x1 - x2) < d )

#define assert_throws(func_call, exception)\
try{(func_call); fprintf(stderr, "Assertion "#func_call" throws "#exception" failed: Nothing thrown\n"); abort();}\
catch(exception &e){}\
catch(...){fprintf(stderr, "Assertion "#func_call" throws "#exception" failed: Other exception thrown\n");}

