#ifndef MISLAV_BRADAC_SEMINAR_SEARCH_H_
#define MISLAV_BRADAC_SEMINAR_SEARCH_H_

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus
#include <immintrin.h>
#ifdef __cplusplus
}
#endif  // __cplusplus

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus
__m256i SearchNormal16(const __m256i **query_score, int query_length,
                       __m256i *f, __m256i *h, __m256i q,
                       __m256i r, __m256i results, const __m256i *z);
#ifdef __cplusplus
}
#endif  // __cplusplus

#endif  // MISLAV_BRADAC_SEMINAR_SEARCH_H_
