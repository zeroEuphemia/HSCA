#include "BlackBoxSolver.h"
#include <iostream>
using namespace CDCLSolver;

Solver::Solver(){
    Minisat::Solver* sv = new Minisat::Solver(true, 1);
    internal_solver = sv;
    assumptions.clear();
}

Solver::~Solver(){ 
    Minisat::Solver* sv = (Minisat::Solver*)internal_solver; 
    delete sv; 
}

void Solver::read_clauses(int nvar, const vector<vector<int> >& clauses){
    // puts("TAT");
    Minisat::Solver* sv = (Minisat::Solver*)internal_solver; 
    Minisat::vec<Minisat::Lit> lits;
    for (const vector<int>& cl: clauses){
        lits.clear();
        for (int rawvar: cl){
            int var = abs(rawvar) - 1;
            while (var >= sv->nVars()) sv->newVar();
            lits.push(Minisat::mkLit(var, rawvar < 0));
        }
        sv->addClause_(lits);
    }

    while (nvar > sv->nVars()){
        sv->newVar();
    }

    sv->user_pol.growTo(sv->nVars(), (Minisat::lbool((uint8_t)2)));
}

void Solver::addVar(int nvar) {    
    Minisat::Solver* sv = (Minisat::Solver*)internal_solver;
    bool changed = false;
    while (nvar >= sv->nVars())
        sv->newVar(), changed = true;
    if (changed)
        sv->user_pol.growTo(sv->nVars(), (Minisat::lbool((uint8_t)2)));
}

void Solver::_addClause(const vector<int>& cls1, const vector<bool>& cls2) {
    Minisat::Solver* sv = (Minisat::Solver*)internal_solver; 
    Minisat::vec<Minisat::Lit> lits;
    lits.clear();

    int nvar = 0, num = cls1.size();
    for (int i = 0; i < num; i ++) {
        // std::cout << cls1[i] << " " << cls2[i] << "\n"; 
        lits.push(Minisat::mkLit(cls1[i], cls2[i]));
        nvar = std::max(nvar, cls1[i]);
    }
    while (nvar >= sv->nVars())
        sv->newVar();
    sv->addClause_(lits);
    sv->user_pol.growTo(sv->nVars(), (Minisat::lbool((uint8_t)2)));
}

void Solver::_removeClause_last() {
    Minisat::Solver* sv = (Minisat::Solver*)internal_solver; 
    // sv->removeClause(sv->clauses[sv->clauses.size() - 1]);
    sv->_removeClause_last();
}

void Solver::add_assumption(int var, int truth_value){
    assumptions.emplace_back(truth_value ? var + 1: -var - 1);
}

void Solver::add_assumption(int lit){
    assumptions.emplace_back(lit);
}

void Solver::clear_assumptions(){
    assumptions.clear();
}

bool Solver::solve(){
    Minisat::Solver* sv = (Minisat::Solver*)internal_solver; 
    if (!sv->simplify()) return false;

    Minisat::vec<Minisat::Lit> assu;
    for (int x: assumptions){
        assu.push(Minisat::mkLit(abs(x) - 1, x < 0));
    }

    bool res = sv->solve(assu);
    return res;
}

void Solver::get_solution(vector<int>& tc){
    Minisat::Solver* sv = (Minisat::Solver*)internal_solver; 
    tc.clear();
    for (int i = 0; i < sv->nVars(); i++){
        if (sv->model[i] == (Minisat::lbool((uint8_t)0))) tc.emplace_back(1);
        else if (sv->model[i] == (Minisat::lbool((uint8_t)1))) tc.emplace_back(0);
        else tc.emplace_back(-1);
    }
}

void Solver:: setPolarity(int vid, bool bit) {
    Minisat::Solver* sv = (Minisat::Solver*)internal_solver; 
    if(bit)
        sv->setPolarity(vid, (Minisat::lbool((uint8_t)1)));
    else
        sv->setPolarity(vid, (Minisat::lbool((uint8_t)0)));
}
