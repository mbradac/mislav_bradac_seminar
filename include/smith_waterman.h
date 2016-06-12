#ifndef MISLAV_BRADAC_SEMINAR_SMITH_WATERMAN_H_
#define MISLAV_BRADAC_SEMINAR_SMITH_WATERMAN_H_

#include <vector>
#include "sequence.h"

namespace SmithWatermanSIMD {

enum ScoreRange {
  kChar, kShort, kDynamic
};

int SmithWaterman(Sequence query, std::vector<Sequence> database,
                  const ScoreMatrix &matrix, int q, int r,
                  ScoreRange score_range, int *results);
}

#endif  // MISLAV_BRADAC_SEMINAR_SMITH_WATERMAN_H_
