#ifndef LOCALSEARCH_OPTIMIZER_INCLUDE_H
#define LOCALSEARCH_OPTIMIZER_INCLUDE_H
#endif

#include <cmath>
#include <fstream>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <thread>

// #include "./core/Solver.h"
// #include "../minisat_ext/BlackBoxSolver.h"
// #include "../minisat_ext/Ext.h"

#include "ConstraintFile.H"
#include "Coverage.h"
// #include "OptionTupleSet.h"
#include "LineVarTupleSet.h"
#include "SAT.h"
#include "Tabu.h"
#include "TupleSet.h"
#include "Valid_check.h"
#include "mersenne.h"

class CoveringArray {
public:
  CoveringArray(const SpecificationFile &specificationFile,
                const ConstraintFile &constraintFile, unsigned long long maxT,
                int seed, int __forced_greedy_percent_, int use_weight_,
                int use_score2, int use_potential_score2_,
                
                const int use_Reset_, 
                const int use_graded_tabu_, 
                const int try_every_, 
                const int _change_, 
                const int score_every_cell_,
                
                const int use_withdraw_,
                const int use_random_step_,
                
                const int thread_num_,
                const int choices_num_,
                const std::string outputFile_);

  void greedyConstraintInitialize();
  void greedyConstraintInitialize2();
  // void actsInitialize(const std::string file_name);
  void arrayInitialize(const std::string file_name);


  void simplifyInitialize();

  void optimize();

private:
  
  std::string outputFile;

  std::vector<std::pair<int, int> > coverByLineIndex;
  int thread_num = 32;
  void arrayInitialize_thread(const std::vector<unsigned> &begin_columns, u_int64_t nums);

  Valid::Validater validater;
  SATSolver satSolver;
  std::vector<bool> option_constrained_indicator;
  Mersenne mersenne;
  const SpecificationFile &specificationFile;
  std::vector<std::vector<unsigned>> bestArray; // = array;
  std::vector<std::vector<unsigned>> array;
  Coverage coverage;
  TupleSet uncoveredTuples;
  std::set<unsigned> varInUncovertuples;
  LineVarTupleSet oneCoveredTuples, twoCoveredTuples;
  Tabu<Entry> entryTabu;
  Tabu<Entry> lineTabu;

  unsigned long long maxTime;
  unsigned long long step;
  clock_t clock_start;

  int sampling_num = 100;
  int __forced_greedy_percent = 0;
  int use_weight = 0, use_score2 = 0;
  int use_potential_score2 = 0;
  std::size_t potential_set_size = 1000000;

  int use_Reset = 1;
  // int reset_times = 2000; // FastCA
  int reset_times = 5000; // FastCA2

  int add_weight_times = 0;
  void Reset();
  void force_tuple(const unsigned lineIndex, const unsigned encode);

  int use_graded_tabu = 1;
  int try_every = 1;
  int _change = 1;
  int score_every_cell = 1;

  int use_withdraw = 1;
  int choices_num = 100;
  void WithdrawTuple (const unsigned encode);

  void random_step ();

  int use_random_step = 0;


  std::vector<InputClause> constr;
  int only_use_init_weight = 0;
  int use_init_weight = 1;

  void adaptive_adjustment();
  void cover(const unsigned encode, const unsigned oldLineIndex);
  void uncover(const unsigned encode, const unsigned oldLineIndex);
  // produce one row at least cover one uncovered tuple.
  // Producing the row without update coverage
  void produceSatRow(std::vector<unsigned> &newLine, const unsigned encode);
  // greedily produce one row at least cover one uncovered tuple.
  // producing the row AND updating coverage
  void mostGreedySatRow(const unsigned lineIndex, const unsigned encode);
  // void mostGreedySatRow(std::vector<unsigned> &newLine, const unsigned encode);
  void mostGreedySatRow2(std::vector<unsigned> &newLine, const unsigned encode);
  void replaceRow(const unsigned lineIndex, const unsigned encode);
  void replaceRow2(const unsigned lineIndex, const unsigned encode);
  void removeUselessRows();
  void removeUselessRows2();
  void removeOneRow();
  void removeOneRowGreedy();
  void removeOneRowRandom();
  long long varScoreOfRow(const unsigned var, const unsigned lineIndex);
  long long varScoreOfRow2(const unsigned var, const unsigned lineIndex);
  std::pair<long long, long long> varScoreOfRow3(const unsigned var, const unsigned lineIndex);
  void replace(const unsigned var, const unsigned lineIndex);

  long long multiVarRow(const std::vector<unsigned> &sortedMultiVars,
                        const unsigned lineIndex, const bool change = false);
  long long multiVarScoreOfRow(const std::vector<unsigned> &sortedMultiVars,
                               const unsigned lineIndex);
  std::pair<long long, long long> multiVarScoreOfRow2(const std::vector<unsigned> &sortedMultiVars,
                                const unsigned lineIndex);
  void multiVarReplace(const std::vector<unsigned> &sortedMultiVars,
                       const unsigned lineIndex);
  
  void force_patching(const unsigned encode);

  void tabugw();
  void tmpPrint();
  bool verify(const std::vector<std::vector<unsigned>> &resultArray);
  bool checkCovered(unsigned encode);


  void check_score2(int Flag);
  // std::vector<std::pair<int, int> > coverByLineIndex;

  // choose s1: -1,   choose s2: 1, same : 0
  int cmp(std::pair<long long, long long> s1, std::pair<long long, long long> s2) {

    long long gap = 0;
    if (s1.first == s2.first) {
      if (s1.second - s2.second <= gap && s1.second - s2.second >= -gap)
        return 0;
      return s1.second < s2.second ? -1 : 1;
    }

    else if(s1.first == s2.first + 1 && s1.second + 10 < s2.second)
      return -1;
    else if(s1.first == s2.first - 1 && s1.second > s2.second + 10)
      return 1;

    return s1.first < s2.first ? -1 : 1;
  }

  void potential_score2_init();

#ifndef NDEBUG
  void print();
#endif
  void t();
};
