#include <cstdio>
extern "C" {
#include <immintrin.h>  // AVX2 and lower
}

#ifdef __AVX2__
const bool kHasAvx2 = true;
#else
const bool kHasAvx2 = false;
#endif

int main() {
  printf("Has AVX2: %d\n", kHasAvx2);
  return 0;
}
