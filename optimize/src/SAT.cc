// Copyright 2008, 2009 Brady J. Garvin

// This file is part of Covering Arrays by Simulated Annealing (CASA).

// CASA is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// CASA is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with CASA.  If not, see <http://www.gnu.org/licenses/>.

// Other SAT solvers can be substituted by altering SAT.H and SAT.C.

//#include "sat/SAT.H"

#include "SAT.h"
#include "Options.h"
#include <iostream>
InputClause::InputClause() : maxVariable(-1) {}

// jkunlin
InputClause::InputClause(const std::vector<InputTerm> &terms)
    : maxVariable(-1) {
  unsigned size = terms.size();
  for (unsigned i = 0; i < size; ++i) {
    append(terms[i]);
  }
}

// InputClause::InputClause(Array<InputTerm>terms) :
//  maxVariable(-1) {
//  unsigned size = terms.getSize();
//  for (unsigned i = 0; i < size; ++i) {
//    append(terms[i]);
//  }
//}

// jkunlin
InputClause::InputClause(const std::vector<unsigned> &symbols)
    : maxVariable(-1) {
  unsigned size = symbols.size();
  for (unsigned i = 0; i < size; ++i) {
    append(InputTerm(false, symbols[i]));
  }
}

// InputClause::InputClause(Array<unsigned>symbols) :
//  maxVariable(-1) {
//  unsigned size = symbols.getSize();
//  for (unsigned i = 0; i < size; ++i) {
//    append(InputTerm(false, symbols[i]));
//  }
//}

InputClause::~InputClause() {}

InputClause::operator vector<InputTerm> &() { return literals; }
InputClause::operator const vector<InputTerm> &() const { return literals; }

void InputClause::clear() {
  maxVariable = -1; // jkunlin
  literals.clear();
}

void InputClause::append(InputTerm term) {
  int variable = term.getVariable();
  if (variable > maxVariable) {
    maxVariable = variable;
  }
  // literals.push(term.isNegated() ? ~Minisat::mkLit(variable) : Minisat::mkLit(variable));
  literals.emplace_back(term);
}

void InputClause::undoAppend() { literals.pop_back(); }

// void SATSolver::reserve(int variables) {
//   // while (variables >= solver.nVars()) {
//   //   solver.newVar();
//   // }
//   cdcl_solver->reserve(variables);
// }

void SATSolver::addClause(InputClause &clause) {
  // reserve(clause.getMaxVariable());
  // solver.addClause(clause);
//   vector<vector<int> > cls;
//   cls.emplace_back(vector<int>());
//   for (InputTerm lit : clause.literals) {
//     if (lit.isNegated())
//       cls[0].emplace_back(-(lit.getVariable() + 1));
//     else
//       cls[0].emplace_back(lit.getVariable() + 1);
//   }
//   for (int item : cls[0])
//     printf("%d ", item);
//   puts("");
  // puts("qaq");
  // cdcl_solver->read_clauses(clause.getMaxVariable(), cls);

    vector<int> cls1;
    vector<bool> cls2;
    for (InputTerm lit : clause.literals) {
        cls1.emplace_back(lit.getVariable());
        cls2.emplace_back(lit.isNegated());
    }
    cdcl_solver->_addClause(cls1, cls2);
}

bool SATSolver::operator()(const InputKnown &known) {
  // if (disable) {
  //   return true;
  // }
  // reserve(known.getMaxVariable());
  // return solver.simplify() && solver.solve(known);

  // return 1;
  if (disable) {
    return true;
  }
  
  cdcl_solver->addVar(known.getMaxVariable());
  cdcl_solver->clear_assumptions();
  for (InputTerm lit : known.literals)
    cdcl_solver->add_assumption(lit.getVariable(), lit.isNegated() ^ 1);
  return cdcl_solver->solve();
}

vector<int> SATSolver::get_force_tuple_tc ( const vector<int> &tc, 
    const vector<unsigned> &target_tuple,
    const vector<pair<int, int> > &value_range ) {
    
    int nvar = (int) tc.size();
    cdcl_solver->addVar(nvar - 1);
    for (int i = 0; i < nvar; i ++)
      cdcl_solver->setPolarity(i, tc[i] == 0);
    cdcl_solver->clear_assumptions();

    int sz = (int) target_tuple.size();
    for (int i = 0; i < sz; i ++) {
      // std::cout << value_range[i].first << " " << value_range[i].second << " " << target_tuple[i] << "\n";
      for (int j = value_range[i].first; j <= value_range[i].second; j ++)
        if (j != target_tuple[i])
          cdcl_solver->add_assumption(j, 0);
      cdcl_solver->add_assumption(target_tuple[i], 1);
    }
    
    bool ret = cdcl_solver->solve();
    if(! ret) {
        std::cout << "c \033[1;31mError: SAT solve failing!\033[0m" << std::endl;
        exit(0);
    }

    vector<int> tc2(nvar, 0);
    cdcl_solver->get_solution(tc2);
    return tc2;
}

bool SATSolver::isUnique(int nvar, const unsigned encode, 
  const std::vector<unsigned> &tuple, const vector<pair<int, int> > &value_range) {

  cdcl_solver->addVar(nvar - 1);
  cdcl_solver->clear_assumptions();
  int sz = (int) tuple.size();
  for (int i = 0; i < sz; i ++) {
    for (int j = value_range[i].first; j <= value_range[i].second; j ++)
      if (j != tuple[i])
        cdcl_solver->add_assumption(j, 0);
    cdcl_solver->add_assumption(tuple[i], 1);
  }
  bool ret = cdcl_solver->solve();
  if(! ret) {
    std::cout << "c \033[1;31mError: SAT solve failing!\033[0m" << std::endl;
    exit(0);
  }
  vector<int> tc1(nvar, 0);
  cdcl_solver->get_solution(tc1);

  for (int i = 0; i < nvar; i ++)
    std::cout << tc1[i] << " ";
  puts("");

  InputClause clause;
  for (int i = 0; i < nvar; i ++) {
    InputTerm lit = InputTerm(tc1[i] != 0, i);
    std::cout << lit << " ";
    clause.append(lit);
  }
  puts("");
  addClause(clause);

  cdcl_solver->clear_assumptions();
  for (int i = 0; i < sz; i ++) {
    for (int j = value_range[i].first; j <= value_range[i].second; j ++)
      if (j != tuple[i])
        cdcl_solver->add_assumption(j, 0);
    cdcl_solver->add_assumption(tuple[i], 1);
  }

  ret = cdcl_solver->solve();
 
  vector<int> tc2(nvar, 0);
  cdcl_solver->get_solution(tc2);
  for (int i = 0; i < nvar; i ++)
    std::cout << tc2[i] << " ";
  puts("");
  
  int diff_num = 0;
  std::cout << "diff: ";
  for (int i = 0; i < nvar; i ++)
    if (tc1[i] != tc2[i])
      std::cout << i << " ", diff_num ++;
  puts("");

  if (! diff_num)
    puts("Error !"), exit(0);

  cdcl_solver->_removeClause_last();
  
  return (! ret);
}