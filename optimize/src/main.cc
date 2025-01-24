#include <iostream>
#include <string>
#include <signal.h>

#include "ConstraintFile.H"
#include "LocalSearch.h"
#include "SpecificationFile.h"

using namespace std;

int main(int argc, char const *argv[]) {
  if (argc < 7) {
    return 1;
  }
  string modelFile(argv[1]);
  string constrFile(argv[2]);
  unsigned long long maxTime = atoi(argv[3]);
  int seed = atoi(argv[4]);
  string initCAFile(argv[5]);
  string outputFile(argv[6]);

  int __forced_greedy_percent = 0;
  int use_weight = 1;
  int use_score2 = 0;
  int use_potential_score2 = 0;

  int use_Reset = 0;
  int use_graded_tabu = 0;
  int try_every = 1;
  int _change = 1;
  int score_every_cell = 0;

  int use_withdraw = 1;
  int use_random_step = 1;
  
  int thread_num = 32;
  int choices_num = 100;
  if (argc >= 8) {
    use_weight = atoi(argv[7]);
    if (argc >= 9) {
      thread_num = atoi(argv[8]);
      if (argc >= 10)
        choices_num = atoi(argv[9]);
    }
  }

  SpecificationFile specificationFile(modelFile);
  ConstraintFile constraintFile(constrFile);
  localSearch(specificationFile, constrFile, maxTime, seed, initCAFile,
  __forced_greedy_percent, use_weight, use_score2, use_potential_score2,
  
  use_Reset, use_graded_tabu, try_every, _change, score_every_cell,
  
  use_withdraw, use_random_step,
  
  thread_num, choices_num, outputFile);
  return 0;
}
