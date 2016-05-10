#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cassert>
extern "C" {
#include <immintrin.h>
}

const int kNumStrings = 10000;
const int kMaxLength = 999999;

void Generate(char *a, char *b, int length) {
  for (int i = 0; i < length; ++i) {
    b[i] = a[i] = rand() % 255;
  }
}

bool CompareString(const char *a, const char *b, int length) {
  for (int i = 0; i < length; ++i) {
    if (*a != *b) return false;
    ++a;
    ++b;
  }
  return true;
}

bool CompareStringSimd(const char *a, const char *b, int length) {
  int i;
  for (i = 0; i < length && (i & 0xf); ++i) {
    if (*a != *b) return false;
    ++a;
    ++b;
  }
  for (; i + 15 < length; i += 16) {
    __m128i av = _mm_load_si128((const __m128i *)a);
    __m128i bv = _mm_load_si128((const __m128i *)b);
    __m128i cmp = _mm_cmpeq_epi32(av, bv);
    int mask = _mm_movemask_epi8(cmp);
    if (mask != 0xFFFF) return false;
    a += 16;
    b += 16;
  }
  for (; i < length; ++i) {
    if (*a != *b) return false;
    ++a;
    ++b;
  }
  return true;
}

int main() {
#ifndef __SSE2__
  fprintf(stderr, "No support for SSE2\n");
  return 1;
#endif
  srand(101);
  static char a[kMaxLength];
  static char b[kMaxLength];
  double simd_time = 0.0;
  double no_simd_time = 0.0;

  int difference_position = -1;
  Generate(a, b, kMaxLength);
  for (int i = 0; i < kNumStrings; ++i) {
    bool same = i < kNumStrings / 2;
    if (!same) {
      if (difference_position != -1) {
        b[difference_position] = a[difference_position];
      }
      difference_position = rand() % kMaxLength;
      b[difference_position] = a[difference_position] + 1;
    }
    {
      clock_t start_time = clock();
      bool simd_res = CompareStringSimd(a, b, kMaxLength);
      simd_time += (double)(clock() - start_time) / CLOCKS_PER_SEC;
      assert(simd_res == same);
    }
    {
      clock_t start_time = clock();
      bool res = CompareString(a, b, kMaxLength);
      no_simd_time += (double)(clock() - start_time) / CLOCKS_PER_SEC;
      assert(res == same);
    }
  }
  printf(
      "Compared %d strings with length %d\nSimd time: %.2lf\nNo simd time: "
      "%.2lf\n",
      kNumStrings, kMaxLength, simd_time, no_simd_time);
  return 0;
}
