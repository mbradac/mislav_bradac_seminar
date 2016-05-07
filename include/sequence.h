#include <string>
#include <vector>

namespace SmithWatermanSIMD {

struct Sequence {
  std::string identifier;
  std::string description;
  std::string sequence;
};

int ParseFasta(const std::string &fasta, std::vector<Sequence> *sequences);
}
