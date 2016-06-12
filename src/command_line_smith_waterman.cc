#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <ctime>
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

void PrintUsage() {
  printf("-q FILE     FILE is query file\n");
  printf("-d FILE     FILE is database file\n");
  printf("-m FILE     FILE is matrix file\n");
  printf("-r N        N is gap extend penalty        [default: 1]\n");
  printf("-e N        N is gap open_extend penalty   [default: 3]\n");
  printf("-x N        N is either 0, 1 or 2          [default: 0\n");
  printf("                 0 - dynamic score range        [default: 0]\n");
  printf("                 1 - 8bit score range\n");
  printf("                 2 - 16bit score range\n");
  printf("-s          silent - don't print results\n");
}

void PrintErrorUsageAndExit(const char *message) {
  printf("%s\n", message);
  printf("\n");
  PrintUsage();
  exit(1);
}

int main(int argc, char **argv) {
  // Initialize.
  char *query_path = nullptr;
  char *database_path = nullptr;
  char *matrix_path = nullptr;
  int r = 1;
  int q = 3;
  ScoreRange score_range = kDynamic;
  bool silent = false;

  // Parse arguments.
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-q") == 0) {
      if (i + 1 == argc) {
        PrintErrorUsageAndExit("File path expected after -q");
      }
      query_path = argv[++i];
    } else if (strcmp(argv[i], "-d") == 0) {
      if (i + 1 == argc) {
        PrintErrorUsageAndExit("Database path expected after -d");
      }
      database_path = argv[++i];
    } else if (strcmp(argv[i], "-m") == 0) {
      if (i + 1 == argc) {
        PrintErrorUsageAndExit("Matrix path expected after -m");
      }
      matrix_path = argv[++i];
    } else if (strcmp(argv[i], "-r") == 0) {
      if (i + 1 == argc) {
        PrintErrorUsageAndExit("Value expected after -r");
      }
      r = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-e") == 0) {
      if (i + 1 == argc) {
        PrintErrorUsageAndExit("Value expected after -e");
      }
      q = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-x") == 0) {
      if (i + 1 == argc) {
        PrintErrorUsageAndExit("Value expected after -x");
      }
      int x = atoi(argv[++i]);
      if (x == 0) {
        score_range = kDynamic;
      } else if (x == 1) {
        score_range = kChar;
      } else if (x == 2) {
        score_range = kShort;
      } else {
        PrintErrorUsageAndExit("Invalid score range value");
      }
    } else if (strcmp(argv[i], "-s") == 0) {
      silent = true;
    } else {
      PrintErrorUsageAndExit("Unknown flag!");
    }
  }

  // Check for missing arguments.
  if (query_path == nullptr) {
    PrintErrorUsageAndExit("Query path not specified");
  }
  if (database_path == nullptr) {
    PrintErrorUsageAndExit("Database path not specified");
  }
  if (matrix_path == nullptr) {
    PrintErrorUsageAndExit("Matrix path not specified");
  }

  // Read, check and parse all requires files.
  std::vector<Sequence> query_vector;
  std::string query_string = ReadFile(query_path);
  if (ParseFasta(query_string, &query_vector) != 0) {
    PrintErrorUsageAndExit("Invalid query file or path");
  }
  if (query_vector.size() != 1U) {
    PrintErrorUsageAndExit("Query file must contain one sequence");
  }
  std::vector<Sequence> database;
  std::string database_string = ReadFile(database_path);
  if (ParseFasta(database_string, &database) != 0) {
    PrintErrorUsageAndExit("Invalid database file or path");
  }
  std::string matrix_string = ReadFile(matrix_path);
  ScoreMatrix matrix;
  if (matrix.Init(matrix_string) != 0) {
    PrintErrorUsageAndExit("Invalid matrix file.\n"
                           "For matrix file example look at test_data/matrix");
  }

  // Solve.
  std::vector<int> results(database.size());
  clock_t start_time = clock();
  int err = SmithWaterman(query_vector[0], database, matrix, q,
                          r, score_range, results.data());
  clock_t end_time = clock();

  // Print.
  if (err) {
    printf("Score range not big enough\n");
    return 0;
  }
  if (silent == false) {
    for (int result : results) {
      printf("%d\n", result);
    }
  }
  fprintf(stderr, "%.3lf\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);
  return 0;
}
