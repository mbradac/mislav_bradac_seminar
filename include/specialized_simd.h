#ifndef MISLAV_BRADAC_SEMINAR_SIMD_INTRINSICS_H_
#define MISLAV_BRADAC_SEMINAR_SIMD_INTRINSICS_H_

#include <limits>
extern "C" {
#include <immintrin.h>
}
#include "search.h"

namespace SmithWatermanSIMD {

#ifdef __AVX2__

const int kRegisterSize = 32;

inline __m256i _mm256_set1(char x) { return _mm256_set1_epi8(x); }
inline __m256i _mm256_set1(short x) { return _mm256_set1_epi16(x); }
inline __m256i _mm256_set1(int x) { return _mm256_set1_epi32(x); }

template<typename T>
inline __m256i SearchNormal(const __m256i **query_score, int query_length,
                            __m256i *f, __m256i *h, __m256i q,
                            __m256i r, __m256i results, const __m256i *z);
template<>
inline __m256i SearchNormal<short> (const __m256i **query_score,
                                    int query_length,__m256i *f, __m256i *h,
                                    __m256i q, __m256i r, __m256i results,
                                    const __m256i *z) {
  return SearchNormal16(query_score, query_length, f, h, q, r, results, z);
}
template<>
inline __m256i SearchNormal<char> (const __m256i **query_score,
                                   int query_length,__m256i *f, __m256i *h,
                                   __m256i q, __m256i r, __m256i results,
                                   const __m256i *z) {
  return SearchNormal8(query_score, query_length, f, h, q, r, results, z);
}

template<typename T>
inline __m256i SearchMasked(const __m256i **query_score, int query_length,
                            __m256i *f, __m256i *h, __m256i q, __m256i r,
                            __m256i results, const __m256i *z, __m256i mask);
template<>
inline __m256i SearchMasked<short> (const __m256i **query_score,
                                    int query_length,__m256i *f, __m256i *h,
                                    __m256i q, __m256i r, __m256i results,
                                    const __m256i *z, __m256i mask) {
  return SearchMasked16(query_score, query_length,
                        f, h, q, r, results, z, mask);
}
template<>
inline __m256i SearchMasked<char> (const __m256i **query_score,
                                   int query_length,__m256i *f, __m256i *h,
                                   __m256i q, __m256i r, __m256i results,
                                   const __m256i *z, __m256i mask) {
  return SearchMasked8(query_score, query_length,
                      f, h, q, r, results, z, mask);
}

inline int extract(char result) {
  return (unsigned char)result ^ 0x80;
}
inline int extract(short result) {
  return (unsigned short)result ^ 0x8000;
}

inline bool did_overflow(char result) {
  return std::numeric_limits<unsigned char>::max() <= (unsigned char)result;
}
inline bool did_overflow(short result) {
  return std::numeric_limits<unsigned short>::max() <= (unsigned short)result;
}

template<typename T>
struct Masks {
};

template<>
struct Masks<char> {
  const static char kZero = 0x80;
};

template<>
struct Masks<short> {
  const static short kZero = 0x8000;
};

template<>
struct Masks<int> {
  const static int kZero = 0x80000000;
};

#endif  // __AVX2__

}
#endif  // MISLAV_BRADAC_SEMINAR_SIMD_INTRINSICS_H_
