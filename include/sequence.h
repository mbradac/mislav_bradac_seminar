#ifndef MISLAV_BRADAC_SEMINAR_SEQUENCE_H_
#define MISLAV_BRADAC_SEMINAR_SEQUENCE_H_

#include <string>
#include <vector>

namespace SmithWatermanSIMD {

struct Sequence {
  std::string identifier;
  std::string description;
  std::string sequence;
};

int ParseFasta(const std::string &fasta, std::vector<Sequence> *sequences);

class ScoreMatrix {
 public:
  ScoreMatrix() {}
  ScoreMatrix &operator=(const ScoreMatrix &) = delete;
  ScoreMatrix(const ScoreMatrix &) = delete;
  int Init(const std::string &score_matrix);
  inline int get_position(int c) const { return position_of_letter_[c]; }
  inline int get_score(int i, int j) const { return matrix_[i][j]; }
  inline int alphabet_size() const { return alphabet_size_; }

 private:
  int alphabet_size_;
  int matrix_[256][256];
  int position_of_letter_[256];
};

int TranslateSequence(Sequence *sequence, const ScoreMatrix &matrix);
}

#endif  // MISLAC_BRADAC_SEMINAR_SEQUENCE_H_
