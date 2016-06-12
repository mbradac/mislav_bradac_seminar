#ifndef MISLAV_BRADAC_SEMINAR_SMITH_WATERMAN_H_
#define MISLAV_BRADAC_SEMINAR_SMITH_WATERMAN_H_

#include <vector>
#include "sequence.h"

namespace SmithWatermanSIMD {

std::vector<int> SmithWaterman(Sequence query, std::vector<Sequence> database,
                               const ScoreMatrix &matrix, int q, int r);
}

#endif  // MISLAV_BRADAC_SEMINAR_SMITH_WATERMAN_H_
