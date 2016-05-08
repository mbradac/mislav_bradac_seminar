#include "sequence.h"

#include <algorithm>
#include <iostream>
#include <istream>
#include <iterator>
#include <string>

using namespace SmithWatermanSIMD;

int main() {
  std::cin >> std::noskipws;
  std::istream_iterator<char> it(std::cin);
  std::istream_iterator<char> end;
  std::string results(it, end);
  std::vector<Sequence> sequences;
  int error = ParseFasta(results, &sequences);
  if (error) {
    printf("Error: %d!\n", error);
    return 0;
  }
  printf("Number of sequences: %zu\n", sequences.size());
  int sum = 0;
  int max = 0;
  for (const auto &sequence : sequences) {
    sum += (int)sequence.sequence.size();
    max = std::max(max, (int)sequence.sequence.size());
  }
  printf("Largest sequence: %d\n", max);
  printf("Average sequence: %lf\n", sum / (double)sequences.size());
  return 0;
}
