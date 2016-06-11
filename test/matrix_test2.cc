#include "sequence.h"

#include <cassert>
#include <string>
#include <vector>

using namespace SmithWatermanSIMD;

std::string matrix_in = "  b d  a \n   1   3  2 \n  1 3 4 \n 6 4 2";
std::vector<std::string> sequence_in = {"abdabd", "abcabc"};
std::vector<std::vector<char>> sequence_out = {{2, 0, 1, 2, 0, 1, 3, 3},
                                               {2, 0, 'c', 'a', 'b', 'c'}};
std::vector<int> out = {0, -1};

int main() {
  ScoreMatrix matrix;
  if (matrix.Init(matrix_in) != 0) {
    printf("WRONG!!!\n");
    return 0;
  }
  for (int i = 0; i < (int)sequence_in.size(); ++i) {
    Sequence sequence;
    sequence.sequence = sequence_in[i];
    if (TranslateSequence(&sequence, matrix, 4) != out[i]) {
      printf("WRONG\n");
    } else if (sequence.sequence.size() != sequence_out[i].size()) {
      printf("WRONG\n");
    } else {
      bool same = true;
      for (int j = 0; j < (int)sequence.sequence.size(); ++j) {
        same &= sequence.sequence[j] == sequence_out[i][j];
      }
      if (!same) {
        printf("WRONG\n");
      } else {
        printf("OK\n");
      }
    }
  }
  return 0;
}
