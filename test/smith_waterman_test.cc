#include <cassert>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "sequence.h"
#include "smith_waterman.h"

using namespace SmithWatermanSIMD;

std::string ReadFile(const char *path) {
  std::ifstream file;
  file.open(path);
  std::stringstream file_stream;
  file_stream << file.rdbuf();
  file.close();
  return file_stream.str();
}

int main() {
  std::vector<std::string> query_paths = {
    "test_data/query/P19930.fasta",
    "test_data/query/P18080.fasta",
    "test_data/query/test2.fasta"};
  std::vector<std::string> database_paths = {
    "test_data/database/uniprot_sprot196.fasta",
    "test_data/database/uniprot_sprot12071.fasta",
    "test_data/database/uniprot_sprot12071.fasta"};
  std::vector<std::string> matrix_paths = {
    "test_data/matrix/blosum50.mat",
    "test_data/matrix/blosum50.mat",
    "test_data/matrix/blosum50.mat"};
  std::vector<std::string> results_paths = {
    "test_data/results/P19930_sprot196_blosum50_r1q3",
    "test_data/results/P18080_sprot12071_blosum50_r1q3",
    "test_data/results/test2_sprot12071_blosum50_r1q3"};
  std::vector<ScoreRange> score_ranges = {kShort, kDynamic, kChar};
  std::vector<int> out = {0, 0, 0};
  std::vector<int> rs = {1, 1, 1};
  std::vector<int> qs = {3, 3, 3};

  for (int t = 0; t < (int)query_paths.size(); ++t) {
    std::vector<Sequence> query_vector;
    std::string query_string = ReadFile(query_paths[t].c_str());
    assert(ParseFasta(query_string, &query_vector) == 0);
    assert(query_vector.size() == 1U);
    std::vector<Sequence> database;
    std::string database_string = ReadFile(database_paths[t].c_str());
    assert(ParseFasta(database_string, &database) == 0);
    std::string matrix_string = ReadFile(matrix_paths[t].c_str());
    ScoreMatrix matrix;
    assert(matrix.Init(matrix_string) == 0);

    std::vector<int> results(database.size());
    clock_t start_time = clock();
    if (SmithWaterman(query_vector[0], database, matrix, qs[t], rs[t],
                      score_ranges[t], results.data()) != out[t]) {
      printf("WRONG!\n");
      return 1;
    }
    clock_t end_time = clock();
    std::vector<int> real_results;
    std::ifstream results_file;
    results_file.open(results_paths[t].c_str());
    int result;
    while (results_file >> result) {
      real_results.push_back(result);
    }
    results_file.close();
    assert(results.size() == real_results.size());
    for (int i = 0; i < (int)results.size(); ++i) {
      if (results[i] != real_results[i]) {
        printf("WRONG!\n");
        return 1;
      }
    }
    printf("OK %.3lf\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);
  }
  return 0;
}
