#include "sequence.h"

#include <cassert>
#include <string>
#include <vector>

using namespace SmithWatermanSIMD;

std::vector<std::string> in = {"", "ab 1", "  a  b c\n 1 2  3 \n 4 5 6",
                               "  a b c\n 1 2  3 \n 4 5 6\n 7 8",
                               "  b d  a \n   1   3  2 \n  1 3 4 \n 6 4 2"};

std::vector<int> out = {-1, -2, -3, -3, 0};

int main() {
  ScoreMatrix matrix;
  for (int i = 0; i < (int)in.size(); ++i) {
    if (matrix.Init(in[i]) != out[i]) {
      printf("WRONG\n");
    } else {
      printf("OK\n");
    }
  }
  return 0;
}
