#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstring>
#include <cstdint>
extern "C" {
#include <immintrin.h>
}

const int kNumCases = 10000;
const int kNumVectors = 10000;
const double eps = 1e-5;

float RandFloat() { return (float)rand() / RAND_MAX; }

inline void Mult(const float a[4], const float b[4][4], float c[4]) {
  c[0] = a[0] * b[0][0] + a[1] * b[1][0] + a[2] * b[2][0] + a[3] * b[3][0];
  c[1] = a[0] * b[0][1] + a[1] * b[1][1] + a[2] * b[2][1] + a[3] * b[3][1];
  c[2] = a[0] * b[0][2] + a[1] * b[1][2] + a[2] * b[2][2] + a[3] * b[3][2];
  c[3] = a[0] * b[0][3] + a[1] * b[1][3] + a[2] * b[2][3] + a[3] * b[3][3];
}

inline __m128 MultSimd(const __m128 a, const __m128 b[4]) {
  __m128 a0 = _mm_shuffle_ps(a, a, 0x00);
  __m128 a1 = _mm_shuffle_ps(a, a, 0x55);
  __m128 a2 = _mm_shuffle_ps(a, a, 0xaa);
  __m128 a3 = _mm_shuffle_ps(a, a, 0xff);
  __m128 a0b0 = _mm_mul_ps(a0, b[0]);
  __m128 a1b1 = _mm_mul_ps(a1, b[1]);
  __m128 a2b2 = _mm_mul_ps(a2, b[2]);
  __m128 a3b3 = _mm_mul_ps(a3, b[3]);
  __m128 a01b01 = _mm_add_ps(a0b0, a1b1);
  __m128 a23b23 = _mm_add_ps(a2b2, a3b3);
  return _mm_add_ps(a01b01, a23b23);
}

int main() {
#ifndef __SSE__
  fprintf(stderr, "No support for SSE2\n");
  return 1;
#endif
  srand(101);
  clock_t simd_time = 0;
  clock_t no_simd_time = 0;

  static float a[kNumVectors][4];
  float b[4][4];
  static float c[kNumVectors][4];

  for (int k = 0; k < kNumCases; ++k) {
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        b[i][j] = RandFloat();
      }
    }
    for (int i = 0; i < kNumVectors; ++i) {
      for (int j = 0; j < 4; ++j) {
        a[i][j] = RandFloat();
      }
    }
    memset(c, 0, sizeof c);

    static __m128 av[kNumVectors];
    __m128 bv[4];
    static __m128 cv[kNumVectors];
    for (int i = 0; i < kNumVectors; ++i) {
      av[i] = _mm_set_ps(a[i][3], a[i][2], a[i][1], a[i][0]);
    }
    for (int i = 0; i < 4; ++i) {
      bv[i] = _mm_set_ps(b[i][3], b[i][2], b[i][1], b[i][0]);
    }

    clock_t start_time = clock();
    for (int i = 0; i < kNumVectors; ++i) {
      cv[i] = MultSimd(av[i], bv);
    }
    clock_t end_time = clock();
    simd_time += end_time - start_time;

    start_time = clock();
    for (int i = 0; i < kNumVectors; ++i) {
      Mult(a[i], b, c[i]);
    }
    end_time = clock();
    no_simd_time += end_time - start_time;
    for (int i = 0; i < kNumVectors; ++i) {
      float _simd_c[8];
      float *simd_c = (float *)((intptr_t)(_simd_c + 4) & ~(0xf));
      _mm_store_ps(simd_c, cv[i]);
      for (int j = 0; j < 4; ++j) {
        assert(abs(simd_c[j] - c[i][j]) < eps);
      }
    }
  }
  printf(
      "Multplied %d matrices.\nSimd time: %.2lf\nNo simd time: "
      "%.2lf\n",
      kNumCases * kNumVectors, (double)simd_time / CLOCKS_PER_SEC,
      (double)no_simd_time / CLOCKS_PER_SEC);
  return 0;
}
