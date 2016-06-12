#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "smith_waterman.h"
#include "specialized_simd.h"
extern "C" {
#include <immintrin.h>
}

namespace SmithWatermanSIMD {

#ifdef __AVX2__

const int kBucketSize = 1024;

template <typename T>
int Search(const Sequence &query, const Sequence *database, int n,
           const ScoreMatrix &matrix, int q, int r, int *all_results) {
  if (n <= 0) {
    return 0;
  }

  // Initialize.
  const int kNumParallel = kRegisterSize / sizeof(T);
  int i;
  int done = 0;
  int query_length = query.sequence.size();
  int idx[kNumParallel] = {0};
  int who[kNumParallel];
  memset(who, -1, sizeof who);
  __m256i z = _mm256_set1(Masks<T>::kZero);
  T unpacked_results[kNumParallel] __attribute__((aligned(kRegisterSize)));
  T unpacked_mask[kNumParallel] __attribute__((aligned(kRegisterSize)));
  T score_cell[kNumParallel] __attribute__((aligned(kRegisterSize)));
  __m256i *h, *f, *score;
  assert(posix_memalign((void **)&h, kRegisterSize,
                        sizeof(*h) * (query_length + 1)) == 0);
  assert(posix_memalign((void **)&f, kRegisterSize,
                        sizeof(*f) * (query_length + 1)) == 0);
  assert(posix_memalign(
             (void **)&score, kRegisterSize,
             kTimesUnrolled * sizeof(*score) * matrix.alphabet_size()) == 0);
  const __m256i **query_score =
      (const __m256i **)malloc(sizeof(score) * query_length);
  assert(query_score != nullptr);
  h[0] = z;
  __m256i results = z;
  __m256i mask = z;
  __m256i qv = _mm256_set1((T)q);  // Cast is necessary to call right overload.
  __m256i rv = _mm256_set1((T)r);  // Cast is necessary to call right overload.

  for (int j = 0; j < query_length; ++j) {
    query_score[j] = &score[kTimesUnrolled * (int)query.sequence[j]];
  }
  bool masked = true;
  for (i = 0; i < kNumParallel && i < n; ++i) {
    idx[i] = 0;
    who[i] = i;
  }

  // Solve.
  while (true) {
    // Preapare query_score.
    int database_chars[kNumParallel];
    for (int l = 0; l < kTimesUnrolled; ++l) {
      for (int j = 0; j < kNumParallel; ++j) {
        database_chars[j] =
            who[j] == -1 ? 0 : database[who[j]].sequence[idx[j] + l];
      }
      for (int j = 0; j < matrix.alphabet_size(); ++j) {
        const int *matrix_row = matrix.get_matrix_row(j);
        for (int k = 0; k < kNumParallel; ++k) {
          score_cell[k] = matrix_row[database_chars[k]];
        }
        score[kTimesUnrolled * j + l] =
            _mm256_load_si256((const __m256i *)score_cell);
      }
    }

    // Search.
    if (masked) {
      results = SearchMasked<T>(query_score, query_length, f, h, qv, rv,
                                results, &z, mask);
    } else {
      results =
          SearchNormal<T>(query_score, query_length, f, h, qv, rv, results, &z);
    }

    // Update results, load new sequences if necessary.
    masked = false;
    _mm256_store_si256((__m256i *)unpacked_results, results);
    for (int j = 0; j < kNumParallel; ++j) {
      idx[j] += kTimesUnrolled;
      if (who[j] != -1 && idx[j] == (int)database[who[j]].sequence.size()) {
        all_results[who[j]] = extract(unpacked_results[j]);
        if (did_overflow((T)all_results[who[j]])) {  // Cast is necessary.
          return 1;
        }
        ++done;
        unpacked_mask[j] = Masks<T>::kZero;
        if (i != n) {
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
    if (done == n) break;
  }

  // Deallocate resources and return result.
  free(f);
  free(h);
  free(score);
  free(query_score);
  return 0;
}

int SmithWaterman(Sequence query, std::vector<Sequence> database,
                  const ScoreMatrix &matrix, int q, int r,
                  ScoreRange score_range, int *results) {
  TranslateSequence(&query, matrix, 1);
  for (Sequence &sequence : database) {
    TranslateSequence(&sequence, matrix, kTimesUnrolled);
  }

  if (score_range == kChar) {
    int ret = Search<char>(query, database.data(), (int)database.size(), matrix,
                           q, r, results);
    return ret;
  }
  if (score_range == kShort) {
    return Search<short>(query, database.data(), (int)database.size(), matrix,
                         q, r, results);
  }
  if (score_range == kDynamic) {
    for (int i = 0; i < (int)database.size(); i += kBucketSize) {
      int n = std::min((int)database.size() - i, kBucketSize);
      int err = Search<char>(query, database.data() + i, n, matrix, q, r,
                             results + i);
      if (err) {
        err = Search<short>(query, database.data() + i, n, matrix, q, r,
                            results + i);
        if (err) {
          return err;
        }
      }
    }
  }
  return 0;
}

#else

int SmithWaterman(Sequence query, std::vector<Sequence> database,
                  const ScoreMatrix &matrix, int q, int r,
                  ScoreRange score_range, int *results) {
  (void)query;
  (void)database;
  (void)matrix;
  (void)q;
  (void)r;
  (void)score_range;
  (void)results;
  fprintf(stderr, "No support for AVX2.\n");
  assert(false);
  return 1;
}

#endif  // __AVX2__
}
