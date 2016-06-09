#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "smith_waterman.h"
#include "search.h"
extern "C" {
#include <immintrin.h>
}

namespace SmithWatermanSIMD {

#ifdef __AVX2__

__m256i SearchMasked16(const char *query, int query_length,
                       const __m256i *score, __m256i *fa, __m256i *ha,
                       __m256i q, __m256i r, __m256i results, __m256i mask) {
  __m256i e = _mm256_setzero_si256();
  __m256i z = _mm256_setzero_si256();
  __m256i h = _mm256_setzero_si256();
  results = _mm256_and_si256(results, mask);
  for (int j = 0; j < query_length; ++j) {
    __m256i prev_h = _mm256_and_si256(ha[j + 1], mask);
    __m256i f = _mm256_and_si256(fa[j + 1], mask);
    h = _mm256_adds_epi16(h, score[(int)query[j]]);
    h = _mm256_max_epi16(h, z);
    h = _mm256_max_epi16(h, e);
    h = _mm256_max_epi16(h, f);
    ha[j + 1] = h;
    results = _mm256_max_epi16(results, h);
    e = _mm256_max_epi16(_mm256_subs_epi16(e, r),
                         _mm256_subs_epi16(h, q));
    f = _mm256_max_epi16(_mm256_subs_epi16(f, r),
                         _mm256_subs_epi16(h, q));
    fa[j + 1] = f;
    h = prev_h;
  }
  return results;
}

std::vector<short> SmithWaterman(Sequence query, std::vector<Sequence> database,
                                 const ScoreMatrix &matrix, int q, int r) {
  if ((int)database.size() == 0) {
    return std::vector<short>();
  }
  int i;
  int done = 0;
  int idx[16] = {0};
  int who[16];
  memset(who, -1, sizeof who);
  short unpacked_results[16] __attribute__((aligned(256 / 8)));
  short unpacked_mask[16] __attribute__((aligned(256 / 8)));
  short score_cell[16] __attribute__((aligned(256 / 8)));
  int query_length = query.sequence.size();
  __m256i *h, *f, *score;
  assert(posix_memalign((void **)&h, 256 / 8, 32 * (query_length + 1)) == 0);
  assert(posix_memalign((void **)&f, 256 / 8, 32 * (query_length + 1)) == 0);
  assert(posix_memalign((void **)&score, 256 / 8,
                        32 * matrix.alphabet_size()) == 0);
  h[0] = _mm256_setzero_si256();
  __m256i results = _mm256_setzero_si256();
  __m256i mask = _mm256_setzero_si256();
  __m256i qv = _mm256_set1_epi16(q);
  __m256i rv = _mm256_set1_epi16(r);
  TranslateSequence(&query, matrix);
  std::vector<short> all_results(database.size());
  bool masked = true;

  for (i = 0; i < 16 && i < (int)database.size(); ++i) {
    idx[i] = 0;
    who[i] = i;
    TranslateSequence(&database[i], matrix);
  }

  while (true) {
    for (int j = 0; j < matrix.alphabet_size(); ++j) {
      for (int k = 0; k < 16; ++k) {
        if (who[k] != -1) {
          score_cell[k] = matrix.get_score(
              j, database[who[k]].sequence[idx[k]]);
        }
      }
      score[j] = _mm256_load_si256((const __m256i *)score_cell);
    }
    if (masked) {
      results = SearchMasked16(query.sequence.c_str(), query_length, score,
                               f, h, qv, rv, results, mask);
    } else {
      results = SearchNormal16(query.sequence.c_str(), query_length, score,
                               f, h, qv, rv, results);
    }
    masked = false;
    _mm256_store_si256((__m256i *)unpacked_results, results);
    for (int j = 0; j < 16; ++j) {
      if (who[j] != -1 && ++idx[j] == (int)database[who[j]].sequence.size()) {
        all_results[who[j]] = unpacked_results[j];
        ++done;
        unpacked_mask[j] = 0;
        if (i != (int)database.size()) {
          TranslateSequence(&database[i], matrix);
          idx[j] = 0;
          who[j] = i++;
          masked = true;
        } else {
          who[j] = -1;
        }
      } else {
        unpacked_mask[j] = ~0;
      }
    }
    mask = _mm256_load_si256((__m256i *)unpacked_mask);
    if (done == (int)database.size()) break;
  }

  free(f);
  free(h);
  free(score);
  return all_results;
}

#else

std::vector<short> SmithWaterman(Sequence query, std::vector<Sequence> database,
                                 const ScoreMatrix &matrix, int q, int r) {
  (void)query;
  (void)database;
  (void)matrix;
  (void)q;
  (void)r;
  fprintf(stderr, "No support for AVX2.\n");
  assert(false);
  return std::vector<short>();
}

#endif
}
