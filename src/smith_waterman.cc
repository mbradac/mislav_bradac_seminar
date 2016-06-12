#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "smith_waterman.h"
#include "search.h"
extern "C" {
#include <immintrin.h>
}
#include <ctime>

namespace SmithWatermanSIMD {

#ifdef __AVX2__

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
                        4 * 32 * matrix.alphabet_size()) == 0);
  const __m256i **query_score =
      (const __m256i **)malloc(sizeof(score) * query_length);
  assert(query_score != nullptr);
  __m256i z = _mm256_set1_epi16(0x8000);
  h[0] = z;
  __m256i results = z;
  __m256i mask = z;
  __m256i qv = _mm256_set1_epi16(q);
  __m256i rv = _mm256_set1_epi16(r);
  TranslateSequence(&query, matrix, 1);
  for (int j = 0; j < query_length; ++j) {
    query_score[j] = &score[4 * (int)query.sequence[j]];
  }
  std::vector<short> all_results(database.size());
  bool masked = true;

  for (i = 0; i < 16 && i < (int)database.size(); ++i) {
    idx[i] = 0;
    who[i] = i;
    TranslateSequence(&database[i], matrix, 4);
  }

  while (true) {
    int database_chars[16];
    for (int l = 0; l < 4; ++l) {
      for (int j = 0; j < 16; ++j) {
        database_chars[j] =
            who[j] == -1 ? 0 : database[who[j]].sequence[idx[j] + l];
      }
      for (int j = 0; j < matrix.alphabet_size(); ++j) {
        const int *matrix_row = matrix.get_matrix_row(j);
        for (int k = 0; k < 16; ++k) {
          score_cell[k] = matrix_row[database_chars[k]];
        }
        score[4 * j + l] = _mm256_load_si256((const __m256i *)score_cell);
      }
    }
    if (masked) {
      results = SearchMasked16(query_score, query_length, f,
                               h, qv, rv, results, &z, mask);
    } else {
      results = SearchNormal16(query_score, query_length, f,
                               h, qv, rv, results, &z);
    }
    masked = false;
    _mm256_store_si256((__m256i *)unpacked_results, results);
    for (int j = 0; j < 16; ++j) {
      idx[j] += 4;
      if (who[j] != -1 && idx[j] == (int)database[who[j]].sequence.size()) {
        all_results[who[j]] = unpacked_results[j] & 0x7FFF;
        ++done;
        unpacked_mask[j] = 0x8000;
        if (i != (int)database.size()) {
          TranslateSequence(&database[i], matrix, 4);
          idx[j] = 0;
          who[j] = i++;
          masked = true;
        } else {
          who[j] = -1;
        }
      } else {
        unpacked_mask[j] = 0;
      }
    }
    mask = _mm256_load_si256((__m256i *)unpacked_mask);
    if (done == (int)database.size()) break;
  }

  free(f);
  free(h);
  free(score);
  free(query_score);
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
