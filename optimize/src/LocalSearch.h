#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H

#include "ConstraintFile.H"
#include "SpecificationFile.h"

void localSearch(const SpecificationFile &specificationFile,
                 const ConstraintFile &constrFile,
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
                 const std::string outputFile);

#endif /* end of include guard: LOCALSEARCH_H */
