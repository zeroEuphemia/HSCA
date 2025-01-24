#ifndef SamplingCA_BlackBox_Solver_h
#define SamplingCA_BlackBox_Solver_h

#include <vector>
#include "Solver.h"

using std::vector;

namespace CDCLSolver {

class Solver{
public:
    Solver();
    virtual ~Solver();
    void read_clauses(int nvar, const vector<vector<int> >& clauses);
    void add_assumption(int var, int truth_value);
    void add_assumption(int lit);
    void clear_assumptions();
    bool solve();
    void get_solution(vector<int>& tc);
    void setPolarity(int vid, bool bit);

    void _addClause(const vector<int>& cls1, const vector<bool>& cls2);
    void _removeClause_last();
    void addVar(int nvar);

protected:
    void* internal_solver;
    vector<int> assumptions;
};

}

#endif