#include "LocalSearch.h"
#include "TupleSet.h"

// #include "ActsSolver.h"
#include "CoveringArray.h"
#include <unistd.h>

void localSearch(const SpecificationFile &specificationFile,
                 const ConstraintFile &constraintFile,
                 const unsigned long long maxTime, int seed,
                 const std::string &initCAFile,
                 const int __forced_greedy_percent,
                 const int use_weight,
                 const int use_score2,
                 const int use_potential_score2,
                 const int use_Reset, 
                 const int use_graded_tabu, 
                 const int try_every, 
                 const int _change, 
                 const int score_every_cell,
                 const int use_withdraw,
                 const int use_random_step,
                 
                 const int thread_num,
                 const int choices_num,
                 const std::string outputFile) {

  CoveringArray c(specificationFile, constraintFile, maxTime, seed,
    __forced_greedy_percent, use_weight, use_score2, use_potential_score2,
    
    use_Reset, use_graded_tabu, try_every, _change, score_every_cell,
    
    use_withdraw, use_random_step,
    
    thread_num, choices_num, outputFile);
  // c.greedyConstraintInitialize2();
  // c.greedyConstraintInitialize();
  c.arrayInitialize(initCAFile);

  c.simplifyInitialize();

  // ActsSolver ActsSolver;
  // char filename[L_tmpnam];
  // if (!tmpnam(filename)) {
  //   std::cerr << "tmp file name error" << std::endl;
  //   abort();
  // }
  // std::string acts_res_filename = filename;
  // acts_res_filename += std::to_string(getpid());
  // ActsSolver.solve(specificationFile, constraintFile, acts_res_filename);
  // c.actsInitialize(acts_res_filename);
  // std::string cmd = (std::string) "rm " + acts_res_filename;
  // if (system(cmd.c_str()) != 0) {
  //   std::cerr << "can't remove acts result file" << std::endl;
  //   exit(0);
  // };
  //	return ;
  c.optimize();
}
