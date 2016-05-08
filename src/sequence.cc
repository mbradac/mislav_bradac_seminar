#include "sequence.h"

#include <sstream>

namespace SmithWatermanSIMD {

int ParseFasta(const std::string &fasta, std::vector<Sequence> *sequences) {
  std::istringstream fasta_stream(fasta);
  std::string line;
  Sequence sequence;
  bool is_first = true;

  while (getline(fasta_stream, line)) {
    if (line.size() == 0U) return -1;
    if (line[0] == '>') {
      if (!is_first) {
        if (sequence.sequence == "") return -2;
        sequences->push_back(sequence);
        sequence.identifier = sequence.description = sequence.sequence = "";
      }
      std::istringstream line_stream(line.substr(1));
      line_stream >> sequence.identifier;
      getline(line_stream >> std::ws, sequence.description);
    } else {
      if (is_first == true) return -3;
      sequence.sequence += line;
    }
    is_first = false;
  }
  if (sequence.sequence == "") return -2;
  sequences->push_back(sequence);
  return 0;
}
}
