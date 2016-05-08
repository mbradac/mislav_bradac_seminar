#include "sequence.h"

#include <cassert>
#include <string>
#include <vector>

using namespace SmithWatermanSIMD;

std::vector<std::string> in = {
    ">abc def ghi\ndefsafa\nsafdsafda\n>sda\nsafsdaafdsafa\n>\nasfafa", "fdsa",
    ">fdsafda dsafdsaf", ">fdsa fdsaf\nfdsafdsa\n>fdsafa dafa dsa",
    ">fdsa fdsa\nsfasdaf\nfdsafsda\n\n>sfa fdsa\nfdasfda"};

std::vector<int> err = {0, -3, -2, -2, -1};

std::vector<std::vector<Sequence>> out = {
    {{"abc", "def ghi", "defsafasafdsafda"},
     {"sda", "", "safsdaafdsafa"},
     {"", "", "asfafa"}},
    {},
    {},
    {{"fdsa", "fdsaf", "fdsafdsa"}},
    {}};

int main() {
  for (int i = 0; i < (int)in.size(); ++i) {
    std::vector<Sequence> sequences;
    int error = ParseFasta(in[i], &sequences);
    bool same = true;
    if (error == err[i] && sequences.size() == out[i].size()) {
      for (int j = 0; j < (int)sequences.size(); ++j) {
        same &= sequences[j].identifier == out[i][j].identifier &&
                sequences[j].description == out[i][j].description &&
                sequences[j].sequence == out[i][j].sequence;
      }
    } else {
      same = false;
    }
    if (same) {
      printf("OK\n");
    } else {
      printf("WRONG\n");
    }
  }
  return 0;
}
