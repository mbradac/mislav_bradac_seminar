#include "sequence.h"

#include <cassert>
#include <string>
#include <vector>

using namespace SmithWatermanSIMD;

std::string matrix_in = "  b d  a \n   1   3  2 \n  1 3 4 \n 6 4 2";
std::vector<std::vector<int>> out = {{1, 3, 2}, {1, 3, 4}, {6, 4, 2}};

int main() {
  ScoreMatrix matrix;
  if (matrix.Init(matrix_in) != 0) {
    printf("WRONG!!!\n");
    return 0;
  }
  bool ok = true;
  for (int i = 0; i < ScoreMatrix::kMatrixSize; ++i) {
    for (int j = 0; j < ScoreMatrix::kMatrixSize; ++j) {
      int a = i < 3 && j < 3 ? out[i][j] : -1;
      ok &= a == matrix.get_score(i, j);
    }
  }
  if (ok) {
    printf("OK\n");
  } else {
    printf("WRONG\n");
  }
  return 0;
}
