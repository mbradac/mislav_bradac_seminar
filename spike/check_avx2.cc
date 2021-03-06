#include <cstdio>
extern "C" {
#include <immintrin.h>
}

#ifdef __AVX2__
const bool kHasAvx2 = true;
#else
const bool kHasAvx2 = false;
#endif

#ifdef __SSE4_1__
const bool kHasSse4_1 = true;
#else
const bool kHasSse4_1 = false;
#endif

#ifdef __SSE4_2__
const bool kHasSse4_2 = true;
#else
const bool kHasSse4_2 = false;
#endif

int main() {
  printf("Has AVX2: %d\n", kHasAvx2);
  printf("Has SSE4.1: %d\n", kHasSse4_1);
  printf("Has SSE4.2: %d\n", kHasSse4_2);
  return 0;
}
