#include "CoveringArray.h"

CoveringArray::CoveringArray(const SpecificationFile &specificationFile,
                             const ConstraintFile &constraintFile,
                             unsigned long long maxT, int seed,
                             int __forced_greedy_percent_, int use_weight_, int use_score2_,
                             int use_potential_score2_,
                             
                             const int use_Reset_, 
                             const int use_graded_tabu_, 
                             const int try_every_, 
                             const int _change_, 
                             const int score_every_cell_,
                             
                             const int use_withdraw_,
                             const int use_random_step_,
                             
                             const int thread_num_,
                             const int choices_num_,
                             const std::string outputFile_)

    : validater(specificationFile), satSolver(constraintFile.isEmpty()),
      specificationFile(specificationFile), coverage(specificationFile),
      entryTabu(4), lineTabu(10), maxTime(maxT), oneCoveredTuples(use_weight_, 1, 0), 
      twoCoveredTuples(use_weight_, 2, use_potential_score2_) {
  
  constr = constraintFile.getClauses();

  __forced_greedy_percent = __forced_greedy_percent_;
  use_weight = use_weight_;
  use_score2 = use_score2_;
  use_potential_score2 = use_potential_score2_;

  use_Reset = use_Reset_;
  use_graded_tabu = use_graded_tabu_;
  try_every = try_every_;
  _change = _change_;
  score_every_cell = score_every_cell_;

  use_withdraw = use_withdraw_;
  use_random_step = use_random_step_;

  thread_num = thread_num_;
  choices_num = choices_num_;
  outputFile = outputFile_;

  adaptive_adjustment();

  // std::cout << "__forced_greedy_percent = " << __forced_greedy_percent << "\n";
  // std::cout << "use_weight = " << use_weight << "\n";
  // std::cout << "use_score2 = " << use_score2 << "\n";
  // std::cout << "use_potential_score2 = " << use_potential_score2 << "\n";
  // std::cout << "use_Reset = " << use_Reset << ", reset_times = " << reset_times << "\n";
  // std::cout << "use_graded_tabu = " << use_graded_tabu << "\n";
  // std::cout << "try_every = " << try_every << "\n";
  // std::cout << "_change = " << _change << "\n";
  // std::cout << "score_every_cell = " << score_every_cell << "\n";

  // std::cout << "use_withdraw = " << use_withdraw << "\n";
  // std::cout << "use_random_step = " << use_random_step << "\n";

  // std::cout << "use_init_weight = " << use_init_weight << "\n";
  // std::cout << "only_use_init_weight = " << only_use_init_weight << "\n";

  // std::cout << "thread_num = " << thread_num << "\n";
  // std::cout << "choices_num = " << choices_num << "\n";

  std::cout << "the second pass of HSCA's optimization approach\n";

  clock_start = clock();
  const Options &options = specificationFile.getOptions();
  // add constraint into satSolver
  const std::vector<InputClause> &clauses = constraintFile.getClauses();
  for (unsigned i = 0; i < clauses.size(); ++i) {
    satSolver.addClause(const_cast<InputClause &>(clauses[i]));
  }
  const Valid::Formula &formula = constraintFile.getFormula();
#ifndef NVISIBLE
  formula.print();
#endif
  for (auto &c : formula) {
    validater.addClause(c);
  }

  option_constrained_indicator.clear();
  option_constrained_indicator.resize(options.size(), false);
  for (auto &c : formula) {
    for (auto &lit : c) {
      option_constrained_indicator[options.option(lit.variable())] = true;
    }
  }
  for (unsigned option = 0; option < options.size(); ++option) {
    if ((!option_constrained_indicator[option]) && (!__forced_greedy_percent)) {
      continue;
    }
    InputClause atLeast;
    for (unsigned j = options.firstSymbol(option),
                  limit = options.lastSymbol(option);
         j <= limit; ++j) {
      atLeast.append(InputTerm(false, j));
    }
    satSolver.addClause(atLeast);
    for (unsigned j = options.firstSymbol(option),
                  limit = options.lastSymbol(option);
         j <= limit; ++j) {
      for (unsigned k = j + 1; k <= limit; ++k) {
        InputClause atMost;
        atMost.append(InputTerm(true, j));
        atMost.append(InputTerm(true, k));
        satSolver.addClause(atMost);
      }
    }
  }

  // coverage.initialize(satSolver, option_constrained_indicator);
  // uncoveredTuples.initialize(specificationFile, coverage, true);
  coverage.unconstrained_initialize();
  uncoveredTuples.initialize(specificationFile, coverage);

  mersenne.seed(seed);
  step = 0;
}

void CoveringArray::potential_score2_init() {
  std::vector<unsigned> vec(coverage.tupleCount());
  std::iota(vec.begin(), vec.end(), 0);
  potential_set_size = std::max(potential_set_size, (std::size_t) coverage.tupleCount() / 100);
  unsigned swap_num = std::min(potential_set_size, vec.size());

  for (unsigned i = 0; i < swap_num; i ++) {
    unsigned offset = mersenne.next(vec.size() - i);
    std::swap(vec[i], vec[i + offset]);
  }
  for (unsigned i = 0; i < swap_num; i ++)
    twoCoveredTuples.add_potential_tuple(vec[i]);
}

void CoveringArray::greedyConstraintInitialize() {
  oneCoveredTuples.initialize(specificationFile, array.size());
  for (auto encode : uncoveredTuples) {
    const std::vector<unsigned> tuple = coverage.getTuple(encode);
    for (auto var : tuple) {
      varInUncovertuples.insert(var);
    }
  }
  const Options &options = specificationFile.getOptions();
  unsigned width = options.size();

  while (uncoveredTuples.size()) {
    oneCoveredTuples.addLine(options.allSymbolCount());
    array.push_back(std::vector<unsigned>(width));

    // reproduce it randomly, with at least one tuple covered
    unsigned encode =
        uncoveredTuples.encode(mersenne.next(uncoveredTuples.size()));
    mostGreedySatRow(array.size() - 1, encode);
  }
  entryTabu.initialize(Entry(array.size(), array.size()));
  lineTabu.initialize(Entry(array.size(), 0));
}

void CoveringArray::adaptive_adjustment() {
  const Options &options = specificationFile.getOptions();
  const unsigned &strenth = specificationFile.getStrenth();
  if (strenth == 5 && (options.size() == 30 || options.size() == 61 || options.size() == 27))
    use_weight = 0;
  if (strenth == 6 && (options.size() == 30 || options.size() == 29 || options.size() == 27))
    use_weight = 0;
}

void CoveringArray::arrayInitialize_thread(const std::vector<unsigned> &begin_columns, u_int64_t nums) {
  
  const Options &options = specificationFile.getOptions();
  const unsigned &strenth = specificationFile.getStrenth();

  std::vector<unsigned> tuple(strenth);
  for (unsigned lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
    auto &line = array[lineIndex];

    u_int64_t tmpnum = 0;
    for (std::vector<unsigned> columns = begin_columns;
        tmpnum < nums && columns[strenth - 1] < line.size(); combinadic.next(columns), tmpnum ++) {
      
      for (unsigned i = 0; i < strenth; ++i)
        tuple[i] = line[columns[i]];
      unsigned encode = coverage.encode(columns, tuple);
      if (coverByLineIndex[encode].first == -1)
        coverByLineIndex[encode].first = lineIndex;
      else
        coverByLineIndex[encode].second = lineIndex;
      coverage.cover(encode);
    }
  }
}

void CoveringArray::arrayInitialize(const std::string file_name) {
  const Options &options = specificationFile.getOptions();
  const unsigned &strenth = specificationFile.getStrenth();
  std::ifstream res_file(file_name);
  if (!res_file.is_open()) {
    std::cerr << "file open failed" << std::endl;
    exit(0);
  }
  std::string line;
  while (getline(res_file, line)) {
    array.push_back(std::vector<unsigned>(options.size()));
    std::vector<unsigned> &newRow = *array.rbegin();
    std::istringstream is(line);
    for (unsigned option = 0; option < options.size(); ++option) {
      unsigned value;
      is >> value;
      newRow[option] = value + options.firstSymbol(option);
    }
  }
  res_file.close();

  coverByLineIndex = std::vector<std::pair<int, int> > (coverage.tupleSize(), {-1, -1});

  // std::cout << "thread_num = " << thread_num << "\n";
  if (thread_num > 1) {

    std::vector<std::thread> Threads;

    u_int64_t columns_num = 1;
    for (int i = 0; i < strenth; i ++)
      columns_num *= (options.size() - i);
    for (int i = 1; i <= strenth; i ++)
      columns_num /= i;
    // std::cout << "columns_num = " << columns_num << "\n";

    u_int64_t every = columns_num / thread_num;
    std::vector<u_int64_t> every_num(thread_num, every);
    columns_num %= thread_num;
    for (int i = 0; i < columns_num; i ++)
      every_num[i] ++;
    
    // for (int i = 0; i < thread_num; i ++)
    //   std::cout << every_num[i] << ", ";
    // std::cout << "\n";

    u_int64_t tmpnum = 0;
    int thread_id = 0;

    std::vector<unsigned> begin_columns = combinadic.begin(strenth);
   
    for (std::vector<unsigned> columns = combinadic.begin(strenth);
         columns[strenth - 1] < options.size(); combinadic.next(columns)) {
      
      if (every_num[thread_id] <= 0)
        break ;
      
      tmpnum ++;
      if (tmpnum == every_num[thread_id]) {
        
        Threads.emplace_back(&CoveringArray:: arrayInitialize_thread, this, 
                begin_columns, every_num[thread_id]);
        
        // std::cout << "thread_id = "<< thread_id << "\n";
        // for (auto c : begin_columns)
        //   std::cout << c << ", ";
        // std::cout << "\n";
        // for (auto c : columns)
        //   std::cout << c << ", ";
        // std::cout << "\n";

        tmpnum = 0;
        thread_id ++;
        begin_columns = columns;
        combinadic.next(begin_columns);
      }
    }

    for (auto &thread : Threads)
      thread.join(); 
  }

  else {
    std::vector<unsigned> tuple(strenth);
    for (unsigned lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
      auto &line = array[lineIndex];
      for (std::vector<unsigned> columns = combinadic.begin(strenth);
          columns[strenth - 1] < line.size(); combinadic.next(columns)) {
        for (unsigned i = 0; i < strenth; ++i) {
          tuple[i] = line[columns[i]];
        }
        //			cover(coverage.encode(columns, tuple));
        unsigned encode = coverage.encode(columns, tuple);
        if (coverByLineIndex[encode].first == -1)
          coverByLineIndex[encode].first = lineIndex;
        else
          coverByLineIndex[encode].second = lineIndex;
        coverage.cover(encode);
      }
    }
  }
  
  coverage.set_zero_invalid();
  entryTabu.initialize(Entry(array.size(), array.size()));
  lineTabu.initialize(Entry(array.size(),0));
  
  oneCoveredTuples.initialize(specificationFile, array.size());
  if (use_weight && use_init_weight)
    oneCoveredTuples.set_init_weight(coverage, constr, specificationFile);
  oneCoveredTuples.pushOneCoveredTuple(coverage, coverByLineIndex);

  if (use_score2) {
    twoCoveredTuples.initialize(specificationFile, array.size());
    if (use_potential_score2)
      potential_score2_init();
    twoCoveredTuples.pushTwoCoveredTuple(coverage, coverByLineIndex);
  }
  
  check_score2(0);
}

void CoveringArray::simplifyInitialize() {
  
  // return ;

  // std::cout << "all tupleCount: " << coverage.tupleCount() << "\n";
  // std::cout << "oneCoveredTuples: " << oneCoveredTuples.size() << "\n";
  // std::cout << "twoCoveredTuples: " << twoCoveredTuples.size() << "\n";
  return ;

  const Options &options = specificationFile.getOptions();
  const unsigned &strenth = specificationFile.getStrenth();
  int nvar = options.lastSymbol(options.size() - 1) + 1;

  int unique_num = 0;
  for (auto tupleEncode : oneCoveredTuples) {
    const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
    std::vector<std::pair<int, int> > value_range;
    for (int i = 0; i < (int)tuple.size(); i ++) {
      int l = options.firstSymbol(options.option(tuple[i]));
      int r = options.lastSymbol(options.option(tuple[i]));
      value_range.emplace_back(std::make_pair(l, r));
    }
    
    // std::cout << tupleEncode << " (";
    // for (int i = 0; i < (int)tuple.size(); i ++)
    //   std::cout << tuple[i] << ", ";
    // std::cout << ")";

    bool ret = satSolver.isUnique(nvar, tupleEncode, tuple, value_range);

    // std::cout << "isUnique = " << ret << "\n";

    unique_num += ret;
  }
  // std::cout << unique_num << "/" << oneCoveredTuples.size() << "\n";
}

void CoveringArray::check_score2(int Flag) {

  return ;

  const Options &options = specificationFile.getOptions();
  const unsigned &strenth = specificationFile.getStrenth();
  std::vector<unsigned> tuple(strenth);

  for (unsigned lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
    auto &line = array[lineIndex];
    for (std::vector<unsigned> columns = combinadic.begin(strenth);
         columns[strenth - 1] < line.size(); combinadic.next(columns)) {
      for (unsigned i = 0; i < strenth; ++i) {
        tuple[i] = line[columns[i]];
      }
      unsigned encode = coverage.encode(columns, tuple);
      
      bool flag = false;
      if (coverage.coverCount(encode) == 1) {
        for (unsigned i = 0; i < strenth; ++i) {
          for (auto item : oneCoveredTuples.getECbyLineVar(lineIndex, tuple[i])) {
            if (item.column_index == i && item.encode == encode) {
              flag = true; break ;
            }
          }
        }
      }

      if (coverage.coverCount(encode) == 2) {
        for (unsigned i = 0; i < strenth; ++i) {
          for (auto item : twoCoveredTuples.getECbyLineVar(lineIndex, tuple[i])) {
            if (item.column_index == i && item.encode == encode) {
              flag = true; break ;
            }
          }
        }
      }

      else flag = true;

      if (! flag) {
        // std::cout << "check lineIndex = " << lineIndex << "\n";
        // for (unsigned i = 0; i < strenth; ++i)
        //   std::cout << tuple[i] << " ";
        // std::cout << "encode = " << encode << ", " << coverage.coverCount(encode) << "\n";
        
        // for (unsigned i = 0; i < strenth; ++i)
        //   for (auto item : twoCoveredTuples.getECbyLineVar(lineIndex, tuple[i]))
        //     std::cout << item.encode << " " << item.column_index << "\n";
        
        puts("\nError Meow !");
        std::cout << "Flag = " << Flag << "\n";
        exit(0);
      }
    }
  }
}

// void CoveringArray::actsInitialize(const std::string file_name) {
//   const Options &options = specificationFile.getOptions();
//   const unsigned &strenth = specificationFile.getStrenth();
//   std::ifstream res_file(file_name);
//   if (!res_file.is_open()) {
//     std::cerr << "file open failed" << std::endl;
//     exit(0);
//   }
//   std::string line;
//   while (getline(res_file, line)) {
//     if (line.find("Test Cases") != std::string::npos) {
//       break;
//     }
//   }
//   while (true) {
//     bool begin = false;
//     while (getline(res_file, line)) {
//       if (line[0] == '1') {
//         begin = true;
//         break;
//       }
//     }
//     if (!begin) {
//       break;
//     }
//     array.push_back(std::vector<unsigned>(options.size()));
//     std::vector<unsigned> &newRow = *array.rbegin();
//     for (unsigned option = 0; option < options.size(); ++option) {
//       unsigned value;
//       std::string value_str(
//           line.substr(line.find_last_of('=') + 1, line.size() - 1));
//       value = atoi(value_str.c_str());
//       newRow[option] = value + options.firstSymbol(option);
//       getline(res_file, line);
//     }
//   }
//   res_file.close();

//   clock_t start = clock();
//   std::vector<size_t> coverByLineIndex(coverage.tupleSize());
//   std::vector<unsigned> tuple(strenth);
//   for (unsigned lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
//     auto &line = array[lineIndex];
//     for (std::vector<unsigned> columns = combinadic.begin(strenth);
//          columns[strenth - 1] < line.size(); combinadic.next(columns)) {
//       for (unsigned i = 0; i < strenth; ++i) {
//         tuple[i] = line[columns[i]];
//       }
//       //			cover(coverage.encode(columns, tuple));
//       unsigned encode = coverage.encode(columns, tuple);
//       coverByLineIndex[encode] = lineIndex;
//       coverage.cover(encode);
//     }
//   }
//   coverage.set_zero_invalid();
//   std::cout << "actsInitialize: " << double(clock() - start) / CLOCKS_PER_SEC
//             << std::endl;
//   entryTabu.initialize(Entry(array.size(), array.size()));
//   tmpPrint();
//   oneCoveredTuples.initialize(specificationFile, array.size());
//   oneCoveredTuples.pushOneCoveredTuple(coverage, coverByLineIndex);
// #ifndef NDEBUG
//   int i = 0;
//   for (auto &line : array) {
//     std::cout << i++ << ": ";
//     for (auto v : line) {
//       std::cout << v << ' ';
//     }
//     std::cout << std::endl;
//   }
// #endif
// }

// void CoveringArray::call_acts(const ConstraintFile &constraintFile) {
//	std::string acts_inputfile_name = "input_acts";
//	std::ofstream acts_infile(acts_inputfile_name);
//	if (!acts_infile.is_open()) {
//		std::cerr << "open failed" << std::endl;
//		exit(0);
//	}
//
//	acts_infile << "[System]" << std::endl;
//	acts_infile << "Name: " << acts_inputfile_name << std::endl <<
// std::endl;
//	acts_infile << "[Parameter]" << std::endl;
//
//	const Options &options = specificationFile.getOptions();
//	for (unsigned option = 0; option < options.size(); ++option) {
//		acts_infile << 'p' << option << "(int): ";
//		acts_infile << 0;
//		for (unsigned var = 1; var < options.symbolCount(option); ++var)
//{
//			acts_infile << ',' << var;
//		}
//		acts_infile << std::endl;
//	}
//
//	acts_infile << std::endl << "[Constraint]" << std::endl;
//	const Valid::Formula &formula = constraintFile.getFormula();
//	for (auto &c : formula) {
//		const unsigned option = options.option(c[0].variable());
//		acts_infile << 'p' << option;
//		c[0].is_negative() ? acts_infile << "!=" : acts_infile << "=";
//		acts_infile << c[0].variable() - options.firstSymbol(option);
//		for (unsigned i = 1; i < c.size(); ++i) {
//			acts_infile <<  " || ";
//			const unsigned option = options.option(c[i].variable());
//			acts_infile << 'p' << option;
//			c[i].is_negative() ? acts_infile << "!=" : acts_infile
//<<
//"=";
//			acts_infile << c[i].variable() -
// options.firstSymbol(option);
//		}
//		acts_infile << std::endl;
//	}
//}

void CoveringArray::produceSatRow(std::vector<unsigned> &newLine,
                                  const unsigned encode) {
  const unsigned strength = specificationFile.getStrenth();
  const Options &options = specificationFile.getOptions();
  const unsigned width = options.size();
  assert(width == newLine.size());

  InputKnown known;
  const std::vector<unsigned> &ranTuple = coverage.getTuple(encode);
  const std::vector<unsigned> &ranTupleColumns = coverage.getColumns(encode);
  for (unsigned i = 0; i < strength; ++i) {
    newLine[ranTupleColumns[i]] = ranTuple[i];
    if (option_constrained_indicator[ranTupleColumns[i]]) {
      known.append(InputTerm(false, ranTuple[i]));
    }
  }
  std::vector<bool> columnStarted(width, false);
  std::vector<unsigned> columnBases(width);
  for (unsigned column = 0, passing = 0; column < width; ++column) {
    if (passing < strength && column == ranTupleColumns[passing]) {
      passing++;
      continue;
    }
    columnBases[column] = mersenne.next(options.symbolCount(column));
    newLine[column] = options.firstSymbol(column) + columnBases[column] - 1;
  }
  for (long column = 0, passing = 0; column < width; ++column) {
    if (passing < strength && column == ranTupleColumns[passing]) {
      passing++;
      continue;
    }
    const unsigned firstSymbol = options.firstSymbol(column);
    const unsigned lastSymbol = options.lastSymbol(column);
    while (true) {
      newLine[column]++;
      if (newLine[column] > lastSymbol) {
        newLine[column] = firstSymbol;
      }
      if (newLine[column] == firstSymbol + columnBases[column]) {
        // backtrack
        if (columnStarted[column]) {
          // std::cout << "backtracking" << std::endl;
          columnStarted[column] = false;
          // assign it the value before starting
          newLine[column]--;
          column--;
          while (passing > 0 && column == ranTupleColumns[passing - 1]) {
            column--;
            passing--;
          }
          if (option_constrained_indicator[column]) {
            known.undoAppend();
          }
          // undo column++ of the "for" loop
          column--;
          // the var of parent column is now unabled
          break;
        } else {
          columnStarted[column] = true;
        }
      }
      if (option_constrained_indicator[column]) {
        known.append(InputTerm(false, newLine[column]));
        if (satSolver(known)) {
          break;
        }
        known.undoAppend();
      } else {
        break;
      }
    }
  }
}

void CoveringArray::mostGreedySatRow(const unsigned lineIndex,
                                     const unsigned encode) {
  std::vector<unsigned> &newLine = array[lineIndex];
  const unsigned strength = specificationFile.getStrenth();
  const Options &options = specificationFile.getOptions();
  const unsigned width = options.size();
  assert(width == newLine.size());

  InputKnown known;
  const std::vector<unsigned> &ranTuple = coverage.getTuple(encode);
  const std::vector<unsigned> &ranTupleColumns = coverage.getColumns(encode);
  for (unsigned i = 0; i < strength; ++i) {
    known.append(InputTerm(false, ranTuple[i]));
  }
  cover(encode, lineIndex);

  std::vector<std::set<unsigned>> columnSymbols(width);
  for (unsigned column = 0, passing = 0; column < width; ++column) {
    if (passing < strength && column == ranTupleColumns[passing]) {
      passing++;
      continue;
    }
    for (unsigned symbol = options.firstSymbol(column);
         symbol <= options.lastSymbol(column); ++symbol) {
      known.append(InputTerm(false, symbol));
      if (satSolver(known)) {
        columnSymbols[column].insert(symbol);
      }
      known.undoAppend();
    }
    break;
  }
  std::vector<unsigned> assignment(width);
  for (unsigned i = 0; i < strength; ++i) {
    assignment[ranTupleColumns[i]] = ranTuple[i];
  }
  std::set<unsigned> assignedVar; // it is sorted
  for (auto var : ranTuple) {
    assignedVar.insert(var);
  }
  // note that covering should consider ranTuple
  for (unsigned column = 0, passing = 0; column < width; ++column) {
    if (passing < strength && column == ranTupleColumns[passing]) {
      passing++;
      continue;
    }
    while (true) {
      unsigned maxNewCoverCount = 0;
      std::vector<unsigned> bestVars;
      for (auto var : columnSymbols[column]) {
        // choose best one
        std::vector<unsigned> assignedVarTmp(assignedVar.begin(),
                                             assignedVar.end());
        std::vector<unsigned> tmpTuple(strength);
        std::vector<unsigned> tmpColumns(strength);
        unsigned newCoverCount = 0;
        for (std::vector<unsigned> columns = combinadic.begin(strength - 1);
             columns[strength - 2] < assignedVarTmp.size();
             combinadic.next(columns)) {
          for (unsigned i = 0; i < strength - 1; ++i) {
            tmpTuple[i] = assignedVarTmp[columns[i]];
          }
          tmpTuple[strength - 1] = var;
          std::sort(tmpTuple.begin(), tmpTuple.end());
          for (unsigned i = 0; i < strength; ++i) {
            tmpColumns[i] = options.option(tmpTuple[i]);
          }
          unsigned tmpEncode = coverage.encode(tmpColumns, tmpTuple);
          if (coverage.coverCount(tmpEncode) == 0) {
            newCoverCount++;
          }
        }
        if (newCoverCount > maxNewCoverCount) {
          maxNewCoverCount = newCoverCount;
          bestVars.clear();
          bestVars.push_back(var);
        } else if (newCoverCount == maxNewCoverCount) {
          bestVars.push_back(var);
        }
      }
      if (bestVars.size() == 0) {
        // backtrack, uncover tuples, undoAppend
        column--;
        while (passing > 0 && column == ranTupleColumns[passing - 1]) {
          column--;
          passing--;
        }
        unsigned backtrackVar = assignment[column];
        assignedVar.erase(backtrackVar);
        // uncover tuples
        std::vector<unsigned> tmpTuple(strength);
        std::vector<unsigned> tmpColumns(strength);
        std::vector<unsigned> assignedVarTmp(assignedVar.begin(),
                                             assignedVar.end());
        for (std::vector<unsigned> columns = combinadic.begin(strength - 1);
             columns[strength - 2] < assignedVarTmp.size();
             combinadic.next(columns)) {
          for (unsigned i = 0; i < strength - 1; ++i) {
            tmpTuple[i] = assignedVarTmp[columns[i]];
          }
          tmpTuple[strength - 1] = backtrackVar;
          std::sort(tmpTuple.begin(), tmpTuple.end());
          for (unsigned i = 0; i < strength; ++i) {
            tmpColumns[i] = options.option(tmpTuple[i]);
          }
          uncover(coverage.encode(tmpColumns, tmpTuple), lineIndex);
        }
        // undoAppend
        known.undoAppend();
      } else {
        // break tie randomly
        unsigned ranIndex = mersenne.next(bestVars.size());
        unsigned tmpVar = bestVars[ranIndex];
        known.append(InputTerm(false, tmpVar));
        // cover tuples
        std::vector<unsigned> tmpTuple(strength);
        std::vector<unsigned> tmpColumns(strength);
        std::vector<unsigned> assignedVarTmp(assignedVar.begin(),
                                             assignedVar.end());
        for (std::vector<unsigned> columns = combinadic.begin(strength - 1);
             columns[strength - 2] < assignedVarTmp.size();
             combinadic.next(columns)) {
          for (unsigned i = 0; i < strength - 1; ++i) {
            tmpTuple[i] = assignedVarTmp[columns[i]];
          }
          tmpTuple[strength - 1] = tmpVar;
          std::sort(tmpTuple.begin(), tmpTuple.end());
          for (unsigned i = 0; i < strength; ++i) {
            tmpColumns[i] = options.option(tmpTuple[i]);
          }
          cover(coverage.encode(tmpColumns, tmpTuple), lineIndex);
        }
        assignment[column] = tmpVar;
        assignedVar.insert(tmpVar);
        // close tmpvar
        columnSymbols[column].erase(tmpVar);
        // if it has, initial next column
        column++;
        while (passing < strength && column == ranTupleColumns[passing]) {
          passing++;
          column++;
        }
        column--;
        if (column < width - 1) {
          for (unsigned symbol = options.firstSymbol(column + 1);
               symbol <= options.lastSymbol(column + 1); ++symbol) {
            known.append(InputTerm(false, symbol));
            if (satSolver(known)) {
              columnSymbols[column + 1].insert(symbol);
            }
            known.undoAppend();
          }
        }
        break;
      }
    }
  }
  newLine.assign(assignedVar.begin(), assignedVar.end());
}

void CoveringArray::replaceRow(const unsigned lineIndex,
                               const unsigned encode) {
  std::vector<unsigned> &ranLine = array[lineIndex];
  const unsigned strength = specificationFile.getStrenth();
  std::vector<unsigned> tmpTuple(strength);
  // uncover the tuples
  for (std::vector<unsigned> columns = combinadic.begin(strength);
       columns[strength - 1] < ranLine.size(); combinadic.next(columns)) {
    for (unsigned i = 0; i < strength; ++i) {
      tmpTuple[i] = ranLine[columns[i]];
    }
    uncover(coverage.encode(columns, tmpTuple), lineIndex);
  }
  produceSatRow(ranLine, encode);
  // cover the tuples
  for (std::vector<unsigned> columns = combinadic.begin(strength);
       columns[strength - 1] < ranLine.size(); combinadic.next(columns)) {
    for (unsigned i = 0; i < strength; ++i) {
      tmpTuple[i] = ranLine[columns[i]];
    }
    cover(coverage.encode(columns, tmpTuple), lineIndex);
  }
  entryTabu.initialize(
      Entry(array.size(), specificationFile.getOptions().size()));
  lineTabu.initialize(Entry(array.size(), 0));
}

void CoveringArray::force_tuple(const unsigned lineIndex, const unsigned encode) {

  const std::vector<unsigned> &target_tuple = coverage.getTuple(encode);
  const std::vector<unsigned> &target_columns = coverage.getColumns(encode);
  const Options &options = specificationFile.getOptions();
  std::vector<unsigned> &nowLine = array[lineIndex];

  bool needChange = false;
  for (unsigned i = 0; i < target_tuple.size(); i ++)
    if (nowLine[target_columns[i]] != target_tuple[i]) {
      needChange = true; break ;
    }
  if (! needChange)
    return ;

  int nvar = options.lastSymbol(options.size() - 1) + 1;

  vector<int> tc(nvar, 0);
  for (int i = 0; i < options.size(); i ++)
    tc[nowLine[i]] = 1;
  vector<pair<int, int> > value_range;

  for (int i = 0; i < (int)target_tuple.size(); i ++) {
    int l = options.firstSymbol(options.option(target_tuple[i]));
    int r = options.lastSymbol(options.option(target_tuple[i]));
    value_range.emplace_back(std::make_pair(l, r));
  }
  vector<int> tc2 = satSolver.get_force_tuple_tc(tc, target_tuple, value_range);

  std::vector<unsigned> changedVars; 
  for (unsigned i = 0; i < nvar; i ++)
    if (tc2[i] && nowLine[options.option(i)] != i)
      changedVars.emplace_back(i);
  
  for (unsigned v : changedVars)
    replace(v, lineIndex);
}

void CoveringArray::WithdrawTuple(const unsigned encode) {

  const std::vector<unsigned> &tuple = coverage.getTuple(encode);
  const std::vector<unsigned> &columns = coverage.getColumns(encode);
  const Options &options = specificationFile.getOptions();

  std::vector<std::vector<unsigned> > canLines;
  for (auto &line : bestArray) {
    bool match = true;
    for (size_t i = 0; i < columns.size(); ++i) {
      if (line[columns[i]] != tuple[i]) {
        match = false;
        break;
      }
    }
    if (match) {
      canLines.emplace_back(line);
      break;
    }
  }

  int orgLines_num = std::min(1, choices_num / (int) canLines.size());
  orgLines_num = std::min(orgLines_num, (int) array.size());

  std::vector<int> line_id;
  for (int i = 0; i < (int) array.size(); i ++)
    line_id.emplace_back(i);
  for (int i = 0; i < orgLines_num; i ++) {
    int j = mersenne.next(line_id.size() - i);
    std::swap(line_id[i], line_id[i + j]);
  }

  std::vector<std::pair<unsigned, unsigned> > best_choice;
  long long bestScore = -1e18;
  for (unsigned i = 0; i < canLines.size(); i ++) {
    auto &line = canLines[i];
    for (unsigned j = 0; j < orgLines_num; j ++) {
      auto &orgLine = array[line_id[j]];

      std::vector<unsigned> changedVars; 
      for (unsigned i = 0; i < options.size(); i ++)
        if (line[i] != orgLine[i])
          changedVars.emplace_back(line[i]);
      
      long long tmpScore;
      if (changedVars.size() > 1)
        tmpScore = multiVarScoreOfRow2(changedVars, line_id[j]).first;
      else
        tmpScore = varScoreOfRow3(changedVars[0], line_id[j]).first;
      
      if (tmpScore > bestScore) {
        bestScore = tmpScore;
        best_choice.clear();
        best_choice.emplace_back(std::make_pair(i, line_id[j]));
      }
      else if (tmpScore == bestScore)
        best_choice.emplace_back(std::make_pair(i, line_id[j]));
    }
  }
  
  // std::cout << bestScore << "\n";

  int choice = mersenne.next(best_choice.size());
  
  const unsigned lineIndex = best_choice[choice].second;
  std::vector<unsigned> &ranLine = array[lineIndex];
  const unsigned strength = specificationFile.getStrenth();
  std::vector<unsigned> tmpTuple(strength);
  // uncover the tuples
  for (std::vector<unsigned> columns = combinadic.begin(strength);
       columns[strength - 1] < ranLine.size(); combinadic.next(columns)) {
    for (unsigned i = 0; i < strength; ++i) {
      tmpTuple[i] = ranLine[columns[i]];
    }
    uncover(coverage.encode(columns, tmpTuple), lineIndex);
  }
  
  ranLine = canLines[best_choice[choice].first];

  for (std::vector<unsigned> columns = combinadic.begin(strength);
       columns[strength - 1] < ranLine.size(); combinadic.next(columns)) {
    for (unsigned i = 0; i < strength; ++i) {
      tmpTuple[i] = ranLine[columns[i]];
    }
    cover(coverage.encode(columns, tmpTuple), lineIndex);
  }
  entryTabu.initialize(
      Entry(array.size(), specificationFile.getOptions().size()));
  lineTabu.initialize(Entry(array.size(), 0));
}

void CoveringArray::force_patching(const unsigned encode) {

  std::vector<int> line_id;
  for (int i = 0; i < (int) array.size(); i ++)
    line_id.emplace_back(i);

  int num = std::min((int) array.size(), sampling_num);
  for (int i = 0; i < num; i ++) {
    int j = mersenne.next(line_id.size() - i);
    std::swap(line_id[i], line_id[i + j]);
  }

  const std::vector<unsigned> &target_tuple = coverage.getTuple(encode);
  const std::vector<unsigned> &target_columns = coverage.getColumns(encode);
  const Options &options = specificationFile.getOptions();
  int nvar = options.lastSymbol(options.size() - 1) + 1;

  std::pair<long long, long long> bestScore = {-1e18, -1e18}, bestScore2 = {-1e18, -1e18};
  std::vector<unsigned> best_lineIndex, best_lineIndex2;
  std::vector<vector<unsigned> > best_changedVars, best_changedVars2;
  for (int i = 0; i < num; i ++) {
    unsigned lineIndex = line_id[i];
    std::vector<unsigned> &nowLine = array[lineIndex];

    vector<int> tc(nvar, 0);
    for (int i = 0; i < options.size(); i ++)
      tc[nowLine[i]] = 1;
    vector<pair<int, int> > value_range;

    bool isTabu = false;
    for (int i = 0; i < (int)target_tuple.size(); i ++) {
      int l = options.firstSymbol(options.option(target_tuple[i]));
      int r = options.lastSymbol(options.option(target_tuple[i]));
      value_range.emplace_back(std::make_pair(l, r));

      if (nowLine[options.option(target_tuple[i])] != target_tuple[i] 
        && entryTabu.isTabu(Entry(lineIndex, target_tuple[i]))) {

        isTabu = true;
        break ;
      }
    }  
    if (isTabu)
      continue ;
    vector<int> tc2 = satSolver.get_force_tuple_tc(tc, target_tuple, value_range);

    std::vector<unsigned> changedVars; 
    for (unsigned i = 0; i < nvar; i ++)
      if (tc2[i] && nowLine[options.option(i)] != i) {
        changedVars.emplace_back(i);
        if (entryTabu.isTabu(Entry(lineIndex, i))) {
          isTabu = true;
          break ;
        }
      }

    if (isTabu)
      continue ;
    
    std::pair<long long, long long> tmpScore = {-1, -1};
    if (changedVars.size() > 1)
      tmpScore = multiVarScoreOfRow2(changedVars, lineIndex);
    else
      tmpScore = varScoreOfRow3(changedVars[0], lineIndex);

    if (use_graded_tabu && (! lineTabu.isTabu(Entry(lineIndex, 0)))) {
      if (cmp(bestScore2, tmpScore) == -1) {
        bestScore2 = tmpScore;
        best_lineIndex2.clear(), best_lineIndex2.emplace_back(lineIndex);
        best_changedVars2.clear(), best_changedVars2.emplace_back(changedVars);
      }
      else if (cmp(bestScore2, tmpScore) == 0) {
        best_lineIndex2.emplace_back(lineIndex);
        best_changedVars2.emplace_back(changedVars);
      }
    }

    if (cmp(bestScore, tmpScore) == -1) {
      bestScore = tmpScore;
      best_lineIndex.clear(), best_lineIndex.emplace_back(lineIndex);
      best_changedVars.clear(), best_changedVars.emplace_back(changedVars);
    }
    else if (cmp(bestScore, tmpScore) == 0) {
      best_lineIndex.emplace_back(lineIndex);
      best_changedVars.emplace_back(changedVars);
    }
  }

  if (! best_lineIndex2.empty()) {
    int choice = mersenne.next(best_lineIndex2.size());
    unsigned lineIndex = best_lineIndex2[choice];
    const std::vector<unsigned> &changedVars = best_changedVars2[choice];
    for (unsigned v : changedVars)
      replace(v, lineIndex);
    return ;
  }

  if (! best_lineIndex.empty()) {
    int choice = mersenne.next(best_lineIndex.size());
    unsigned lineIndex = best_lineIndex[choice];
    const std::vector<unsigned> &changedVars = best_changedVars[choice];
    for (unsigned v : changedVars)
      replace(v, lineIndex);
    return ;
  }

  replaceRow2(mersenne.next(array.size()), encode);
}

void CoveringArray::replaceRow2(const unsigned lineIndex,
                                const unsigned encode) {
  
  if (use_withdraw)
    return WithdrawTuple(encode);

  std::vector<unsigned> &ranLine = array[lineIndex];
  const unsigned strength = specificationFile.getStrenth();
  std::vector<unsigned> tmpTuple(strength);
  // uncover the tuples
  for (std::vector<unsigned> columns = combinadic.begin(strength);
       columns[strength - 1] < ranLine.size(); combinadic.next(columns)) {
    for (unsigned i = 0; i < strength; ++i) {
      tmpTuple[i] = ranLine[columns[i]];
    }
    uncover(coverage.encode(columns, tmpTuple), lineIndex);
  }

  const std::vector<unsigned> &target_tuple = coverage.getTuple(encode);
  const std::vector<unsigned> &target_columns = coverage.getColumns(encode);
  
  // if (mersenne.next(100) < __forced_greedy_percent) {
  //   const Options &options = specificationFile.getOptions();
    
  //   int nvar = options.lastSymbol(options.size() - 1) + 1;
  //   vector<int> tc(nvar, 0);
  //   for (int i = 0; i < options.size(); i ++)
  //     tc[ranLine[i]] = 1;

  //   vector<pair<int, int> > value_range;
  //   for (int i = 0; i < (int)target_tuple.size(); i ++) {
  //     int l = options.firstSymbol(options.option(target_tuple[i]));
  //     int r = options.lastSymbol(options.option(target_tuple[i]));
  //     value_range.emplace_back(std::make_pair(l, r));
  //   }
    
  //   vector<int> tc2 = satSolver.get_force_tuple_tc(tc, target_tuple, value_range);

  //   for (int i = 0; i < options.size(); i ++)
  //     ranLine[i] = 0;
  //   for (int i = 0; i < nvar; i ++)
  //     if (tc2[i])
  //       ranLine[options.option(i)] = i;
  // }
  
  // else {
    for (auto &line : bestArray) {
      bool match = true;
      for (size_t i = 0; i < target_columns.size(); ++i) {
        if (line[target_columns[i]] != target_tuple[i]) {
          match = false;
          break;
        }
      }
      if (match) {
        ranLine = line;
        break;
      }
    }
  // }

  // cover the tuples
  for (std::vector<unsigned> columns = combinadic.begin(strength);
       columns[strength - 1] < ranLine.size(); combinadic.next(columns)) {
    for (unsigned i = 0; i < strength; ++i) {
      tmpTuple[i] = ranLine[columns[i]];
    }
    cover(coverage.encode(columns, tmpTuple), lineIndex);
  }
  entryTabu.initialize(
      Entry(array.size(), specificationFile.getOptions().size()));
  lineTabu.initialize(Entry(array.size(), 0));
  
  check_score2(1);
}

void CoveringArray::removeUselessRows2() {
  const Options &options = specificationFile.getOptions();
  const unsigned strength = specificationFile.getStrenth();
  std::vector<unsigned> tmpTuple(strength);

  for (size_t lineIndex = 0; lineIndex < array.size();) {
    if (oneCoveredTuples.oneCoveredCount(lineIndex) == 0) {
      const std::vector<unsigned> &line = array[lineIndex];
      for (std::vector<unsigned> columns = combinadic.begin(strength);
           columns[strength - 1] < options.size(); combinadic.next(columns)) {
        for (unsigned i = 0; i < strength; ++i) {
          tmpTuple[i] = line[columns[i]];
        }
        unsigned encode = coverage.encode(columns, tmpTuple);
        uncover(encode, lineIndex);
      }
      std::swap(array[lineIndex], array[array.size() - 1]);
      for (auto &entry : entryTabu) {
        if (entry.getRow() == array.size() - 1) {
          entry.setRow(lineIndex);
        }
        if (entry.getRow() == lineIndex) {
          entry.setRow(array.size() - 1);
        }
      }
      // meow QAQ
      for (auto &entry : lineTabu) {
        if (entry.getRow() == array.size() - 1) {
          entry.setRow(lineIndex);
        }
        if (entry.getRow() == lineIndex) {
          entry.setRow(array.size() - 1);
        }
      }

      oneCoveredTuples.exchange_row(lineIndex, array.size() - 1);
      oneCoveredTuples.pop_back_row();
      if (use_score2) {
        twoCoveredTuples.exchange_row(lineIndex, array.size() - 1);
        twoCoveredTuples.pop_back_row();
      }
      array.pop_back();
    } else {
      ++lineIndex;
    }
  }

  check_score2(2);
}

void CoveringArray::removeUselessRows() {
  const Options &options = specificationFile.getOptions();
  const unsigned strength = specificationFile.getStrenth();
  std::vector<unsigned> tmpTuple(strength);

  for (unsigned i = 0; i < array.size(); ++i) {
    bool useless = true;
    const std::vector<unsigned> &line = array[i];
    for (std::vector<unsigned> columns = combinadic.begin(strength);
         columns[strength - 1] < options.size(); combinadic.next(columns)) {
      for (unsigned j = 0; j < strength; ++j) {
        tmpTuple[j] = line[columns[j]];
      }
      unsigned encode = coverage.encode(columns, tmpTuple);
      if (coverage.coverCount(encode) == 1) {
        useless = false;
        break;
      }
    }
    if (useless) {
      for (std::vector<unsigned> columns = combinadic.begin(strength);
           columns[strength - 1] < options.size(); combinadic.next(columns)) {
        for (unsigned j = 0; j < strength; ++j) {
          tmpTuple[j] = line[columns[j]];
        }
        unsigned encode = coverage.encode(columns, tmpTuple);
        uncover(encode, i);
      }
      std::swap(array[i], array[array.size() - 1]);
      for (auto &entry : entryTabu) {
        if (entry.getRow() == array.size() - 1) {
          entry.setRow(i);
        }
        if (entry.getRow() == i) {
          entry.setRow(array.size() - 1);
        }
      }
      for (auto &entry : lineTabu) {
        if (entry.getRow() == array.size() - 1) {
          entry.setRow(i);
        }
        if (entry.getRow() == i) {
          entry.setRow(array.size() - 1);
        }
      }
      oneCoveredTuples.exchange_row(i, array.size() - 1);
      oneCoveredTuples.pop_back_row();
      if (use_score2) {
        twoCoveredTuples.exchange_row(i, array.size() - 1);
        twoCoveredTuples.pop_back_row();
      }
      array.pop_back();
      --i;
    }
  }
}

void CoveringArray::removeOneRow() {
  const Options &options = specificationFile.getOptions();
  const unsigned strength = specificationFile.getStrenth();
  unsigned minUncoverCount = std::numeric_limits<unsigned>::max();
  std::vector<unsigned> bestRowIndex;
  std::vector<unsigned> tmpTuple(strength);
  clock_t start = clock();
  for (unsigned i = 0; i < array.size(); ++i) {
    unsigned uncoverCount = 0;
    const std::vector<unsigned> &line = array[i];
    for (std::vector<unsigned> columns = combinadic.begin(strength);
         columns[strength - 1] < options.size(); combinadic.next(columns)) {
      for (unsigned j = 0; j < strength; ++j) {
        tmpTuple[j] = line[columns[j]];
      }
      unsigned encode = coverage.encode(columns, tmpTuple);
      if (coverage.coverCount(encode) == 1) {
        uncoverCount++;
      }
    }
    if (uncoverCount < minUncoverCount) {
      minUncoverCount = uncoverCount;
      bestRowIndex.clear();
      bestRowIndex.push_back(i);
    } else if (uncoverCount == minUncoverCount) {
      bestRowIndex.push_back(i);
    }
  }
  std::cout << "select time: " << double(clock() - start) / CLOCKS_PER_SEC
            << std::endl;

  unsigned rowToremoveIndex = bestRowIndex[mersenne.next(bestRowIndex.size())];
  for (std::vector<unsigned> columns = combinadic.begin(strength);
       columns[strength - 1] < options.size(); combinadic.next(columns)) {
    for (unsigned j = 0; j < strength; ++j) {
      tmpTuple[j] = array[rowToremoveIndex][columns[j]];
    }
    unsigned encode = coverage.encode(columns, tmpTuple);
    uncover(encode, rowToremoveIndex);
  }

  std::swap(array[array.size() - 1], array[rowToremoveIndex]);
  for (auto &entry : entryTabu) {
    if (entry.getRow() == array.size() - 1) {
      entry.setRow(rowToremoveIndex);
    }
    if (entry.getRow() == rowToremoveIndex) {
      entry.setRow(array.size() - 1);
    }
  }
  for (auto &entry : lineTabu) {
    if (entry.getRow() == array.size() - 1) {
      entry.setRow(rowToremoveIndex);
    }
    if (entry.getRow() == rowToremoveIndex) {
      entry.setRow(array.size() - 1);
    }
  }
  array.pop_back();
}

void CoveringArray::removeOneRowGreedy() {
  const Options &options = specificationFile.getOptions();
  const unsigned strength = specificationFile.getStrenth();

  unsigned rowToremoveIndex = 0;
  unsigned minOneCoveredCount = oneCoveredTuples.oneCoveredCount_Weight(0);
  std::vector<unsigned> bestLines;
  bestLines.emplace_back(0);

  for (size_t lineIndex = 1; lineIndex < array.size(); ++lineIndex) {
    unsigned oneCoveredCount = oneCoveredTuples.oneCoveredCount_Weight(lineIndex);
    if (minOneCoveredCount > oneCoveredCount) {
      minOneCoveredCount = oneCoveredCount;
      // rowToremoveIndex = lineIndex;
      bestLines.clear();
      bestLines.emplace_back(lineIndex);
    }
    else if(minOneCoveredCount == oneCoveredCount)
      bestLines.emplace_back(lineIndex);
  }
  rowToremoveIndex = bestLines[mersenne.next(bestLines.size())];

  std::vector<unsigned> tmpTuple(strength);
  for (std::vector<unsigned> columns = combinadic.begin(strength);
       columns[strength - 1] < options.size(); combinadic.next(columns)) {
    for (unsigned j = 0; j < strength; ++j) {
      tmpTuple[j] = array[rowToremoveIndex][columns[j]];
    }
    unsigned encode = coverage.encode(columns, tmpTuple);
    uncover(encode, rowToremoveIndex);
  }

  std::swap(array[array.size() - 1], array[rowToremoveIndex]);
  oneCoveredTuples.exchange_row(rowToremoveIndex, array.size() - 1);
  oneCoveredTuples.pop_back_row();
  if (use_score2) {
    twoCoveredTuples.exchange_row(rowToremoveIndex, array.size() - 1);
    twoCoveredTuples.pop_back_row();
  }
  for (auto &entry : entryTabu) {
    if (entry.getRow() == array.size() - 1) {
      entry.setRow(rowToremoveIndex);
    }
    if (entry.getRow() == rowToremoveIndex) {
      entry.setRow(array.size() - 1);
    }
  }
  for (auto &entry : lineTabu) {
    if (entry.getRow() == array.size() - 1) {
      entry.setRow(rowToremoveIndex);
    }
    if (entry.getRow() == rowToremoveIndex) {
      entry.setRow(array.size() - 1);
    }
  }
  array.pop_back();
}

void CoveringArray::removeOneRowRandom() {

  const Options &options = specificationFile.getOptions();
  const unsigned strength = specificationFile.getStrenth();

  unsigned rowToremoveIndex = mersenne.next(array.size());
  std::vector<unsigned> tmpTuple(strength);
  for (std::vector<unsigned> columns = combinadic.begin(strength);
       columns[strength - 1] < options.size(); combinadic.next(columns)) {
    for (unsigned j = 0; j < strength; ++j) {
      tmpTuple[j] = array[rowToremoveIndex][columns[j]];
    }
    unsigned encode = coverage.encode(columns, tmpTuple);
    uncover(encode, rowToremoveIndex);
    // puts("0.0");
  }

  // puts("> <");
  std::swap(array[array.size() - 1], array[rowToremoveIndex]);
  oneCoveredTuples.exchange_row(rowToremoveIndex, array.size() - 1);
  oneCoveredTuples.pop_back_row();
  if (use_score2) {
    twoCoveredTuples.exchange_row(rowToremoveIndex, array.size() - 1);
    twoCoveredTuples.pop_back_row();
  }
  for (auto &entry : entryTabu) {
    if (entry.getRow() == array.size() - 1) {
      entry.setRow(rowToremoveIndex);
    }
    if (entry.getRow() == rowToremoveIndex) {
      entry.setRow(array.size() - 1);
    }
  }
  for (auto &entry : lineTabu) {
    if (entry.getRow() == array.size() - 1) {
      entry.setRow(rowToremoveIndex);
    }
    if (entry.getRow() == rowToremoveIndex) {
      entry.setRow(array.size() - 1);
    }
  }
  array.pop_back();
  check_score2(3);
}

void CoveringArray::random_step () {

  const Options &options = specificationFile.getOptions();
  std::vector<std::pair<unsigned, unsigned> > operations;
  
  for (unsigned i = 0; i < array.size(); i ++) {
    auto &line = array[i];
    for (unsigned j = 0; j < options.size(); j ++)
      for (unsigned k = options.firstSymbol(j); k <= options.lastSymbol(j); k ++)
        if (line[j] != k)
          operations.emplace_back(std::make_pair(i, k));
  }
  
  std::random_shuffle(operations.begin(), operations.end());

  for (auto operation : operations) {
      
      unsigned lineIndex = operation.first;
      unsigned diffVar = operation.second;
      unsigned diffOption = options.option(diffVar);

      // std::cout << lineIndex << ", " << diffVar << ", " << diffOption << "\n";

      // Tabu
      if (entryTabu.isTabu(Entry(lineIndex, diffOption)))
        continue;
      // my check
      if (!validater.valida_change(array[lineIndex], diffVar))
        continue;
      
      replace(diffVar, lineIndex);
      break ;
  }

}

void CoveringArray::optimize() {
  int limit = 0;
  // std::cout << "optimize !\n";

  // while(array.size() > 0) {
  //   removeOneRowRandom();
  //   // removeOneRowGreedy();
  //   std::cout << uncoveredTuples.size() << ", ";
  // }
  // std::cout << "\n";
  // exit(0);

  while (true) {
    if ((double)(clock() - clock_start) / CLOCKS_PER_SEC > maxTime) {
      break;
    }
    if (uncoveredTuples.size() == 0) {
      removeUselessRows2();
      bestArray = array;
      tmpPrint();
      //			removeOneRow();
      removeOneRowRandom();
      // removeOneRowGreedy();

      add_weight_times = 0;
    }

    // if (use_random_step && mersenne.nextClosed() < 0.001) {
    //   puts("random_step");
    //   random_step();
    // }

    //		clock_t start = clock();
    // tabuStep();
    //		tabuStep2();
    //		tabuStep3();
    // tabuStep4();
    // tabuZero();
    tabugw();
    //		std::cout << "tabuStep time: " <<  double (clock() - start) /
    // CLOCKS_PER_SEC << std::endl;

    step++;
    continue;
  }

  if (uncoveredTuples.size() == 0) {
    removeUselessRows2();
    bestArray = array;
    tmpPrint();
  }

  // 顺便把结果打印
  // FILE *fout = fopen("a.out", "w");
  FILE *fout = fopen(outputFile.c_str(), "w");
  for (unsigned i = 0; i < bestArray.size(); ++i) {
    // std::cerr << i << "th  ";
    // for (auto x : bestArray[i]) {
    //   std::cerr << ' ' << x;
    // }
    // std::cerr << std::endl;
    int sz = (int) bestArray[i].size();
    for (int j = 0; j < sz; ++j)
      fprintf(fout, "%u%c", bestArray[i][j], " \n"[j == sz - 1]);
  }
  fclose(fout);
  
  // std::cout << "begin verify !\n";
  // if (!verify(bestArray)) {
    // std::cout << "wrong answer!!!!!" << std::endl;
    // return;
  // }
  // std::cout << "total steps: " << step << std::endl;

  //#ifndef NDEBUG
  std::cerr << "********Debuging CoveringArray::optimize*********" << std::endl;
  	// std::cerr << "printing bestArray..." << std::endl;
    // const Options &options = specificationFile.getOptions();
  	// for (unsigned i = 0; i < bestArray.size(); ++i) {
  	// 	std::cerr << i << "th  ";
    //     int j = 0;
  	// 	for (auto x : bestArray[i]) {
  	// 		std::cerr << ' ' << x - options.firstSymbol(j);
    //         j ++;
  	// 	}
  	// 	std::cerr << std::endl;
  	// }
  std::cerr << "total size : " << bestArray.size() << std::endl;
  std::cerr << "********End of Debuing CoveringArray::optimize********"
            << std::endl;
  //#endif
}

void CoveringArray::Reset() {

  // return ;
  // std::cout << "begin Reset: uncoveredTuples.size() = " << uncoveredTuples.size() << "\n";

  std::vector<unsigned> available_lineIndex;
  for (unsigned i = 0; i < array.size(); i ++)
    if (! lineTabu.isTabu(Entry(i, 0)))
      available_lineIndex.emplace_back(i);
  
  if (uncoveredTuples.size() > available_lineIndex.size()) {
    // std::cout << uncoveredTuples.size() << " " << available_lineIndex.size() << "\n";
    // puts("Reset Error !");
    // exit(0);

    puts("Reset Fail !");
    add_weight_times = 0;
    return ;
  }

  puts("Reset Success !");

  for (unsigned i = 0; i < uncoveredTuples.size(); i ++) {
    int offset = mersenne.next(available_lineIndex.size() - i);
    std::swap(available_lineIndex[i], available_lineIndex[i + offset]);
  }

  std::vector<unsigned> now_uncoveredTuples;
  for (auto tupleEncode : uncoveredTuples)
    now_uncoveredTuples.emplace_back(tupleEncode);

  unsigned i = 0;
  for (auto tupleEncode : now_uncoveredTuples) {
    unsigned ranLineIndex = available_lineIndex[i];
    unsigned ran = mersenne.next(100);
    if (ran < 95)
      force_tuple(ranLineIndex, tupleEncode);
    else
      replaceRow2(ranLineIndex, tupleEncode);
    i ++;
  }

  add_weight_times = 0;
  entryTabu.initialize(Entry(array.size(), specificationFile.getOptions().size()));
  lineTabu.initialize(Entry(array.size(), 0));
}

void CoveringArray::tabugw() {
  unsigned base = mersenne.next(uncoveredTuples.size());

  std::vector<unsigned> firstBestRows, firstBestRows2;
  std::vector<unsigned> firstBestVars, firstBestVars2;
  std::vector<unsigned> bestRows, bestRows2;
  std::vector<unsigned> bestVars, bestVars2;
  std::pair<long long, long long> bestScore = {-1e18, -1e18}, bestScore2 = {-1e18, -1e18};
  
  
  // long long sumWeight = 0;
  for (size_t i = 0; i < uncoveredTuples.size(); ++i) {
    const unsigned tupleEncode =
        uncoveredTuples.encode((base + i) % uncoveredTuples.size());
    
    // if (use_weight)
    //   sumWeight += oneCoveredTuples.getWeight(tupleEncode);
    
    const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
    const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);
    for (unsigned lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
      std::vector<unsigned> &line = array[lineIndex];
      unsigned diffCount = 0;
      unsigned diffVar;
      for (unsigned i = 0; i < tuple.size(); ++i) {
        if (line[columns[i]] != tuple[i]) {
          diffCount++;
          diffVar = tuple[i];
        }
      }
      if (diffCount > 1) {
        continue;
      }
      unsigned diffOption = specificationFile.getOptions().option(diffVar);
      // Tabu
      if (entryTabu.isTabu(Entry(lineIndex, diffOption))) {
        continue;
      }
      // my check
      if (!validater.valida_change(array[lineIndex], diffVar)) {
        continue;
      }
      std::pair<long long, long long> tmpScore = varScoreOfRow3(diffVar, lineIndex);

      if (use_graded_tabu && (! lineTabu.isTabu(Entry(lineIndex, 0)))) {
        if (cmp(bestScore2, tmpScore) == -1) {
          bestScore2 = tmpScore;
          if (i == 0) {
            firstBestRows2.clear();
            firstBestRows2.push_back(lineIndex);
            firstBestVars2.clear();
            firstBestVars2.push_back(diffVar);
          }
          bestRows2.clear();
          bestRows2.push_back(lineIndex);
          bestVars2.clear();
          bestVars2.push_back(diffVar);
        } 
        else if (cmp(bestScore2, tmpScore) == 0) {
          bestRows2.push_back(lineIndex);
          bestVars2.push_back(diffVar);
          if (i == 0) {
            firstBestRows2.push_back(lineIndex);
            firstBestVars2.push_back(diffVar);
          }
        }
      }

      if (cmp(bestScore, tmpScore) == -1) {
        bestScore = tmpScore;
        if (i == 0) {
          firstBestRows.clear();
          firstBestRows.push_back(lineIndex);
          firstBestVars.clear();
          firstBestVars.push_back(diffVar);
        }
        bestRows.clear();
        bestRows.push_back(lineIndex);
        bestVars.clear();
        bestVars.push_back(diffVar);
      } else if (cmp(bestScore, tmpScore) == 0) {
        bestRows.push_back(lineIndex);
        bestVars.push_back(diffVar);
        if (i == 0) {
          firstBestRows.push_back(lineIndex);
          firstBestVars.push_back(diffVar);
        }
      }
    }
  }

  // std::cout << bestScore << " " << uncoveredTuples.size() << " " << sumWeight << "\n";
  // std::cout << bestScore.first << " " << bestScore.second << "\n";
  // std::cout << bestScore.first << " " << bestScore.second << "\n";

  if (use_graded_tabu && bestScore2.first > 0) {
    unsigned ran = mersenne.next(bestRows2.size());
    replace(bestVars2[ran], bestRows2[ran]);
    return ;
  }

  if (bestScore.first > 0) {
    unsigned ran = mersenne.next(bestRows.size());
    replace(bestVars[ran], bestRows[ran]);
    return ;
  }
  
  // weight
  if (use_weight && (! only_use_init_weight)) {

    // std::cout << uncoveredTuples.size() << "\n";
    add_weight_times ++;
    // if (add_weight_times % 500 == 0)
    //   std::cout << "add_weight_times = " << add_weight_times << "\n";
    if (use_Reset && add_weight_times >= reset_times)
      Reset();

    for (auto tupleEncode : uncoveredTuples) {
      oneCoveredTuples.addWeight(tupleEncode);
      if (use_score2) {
        if ((! use_potential_score2) || twoCoveredTuples.is_potential(tupleEncode))
          twoCoveredTuples.addWeight(tupleEncode);
      }
    }
  }

  if (use_random_step && mersenne.nextClosed() < 0.001) {
    // puts("random_step");
    random_step();
  }

  if (mersenne.next(100) < __forced_greedy_percent) {
    force_patching(uncoveredTuples.encode(base));
    return;
  }

  else if ((! __forced_greedy_percent) && mersenne.nextClosed() < 0.001) {
    replaceRow2(mersenne.next(array.size()), uncoveredTuples.encode(base));
    return;
  }

  if ((! _change || mersenne.next(100) < 90) && use_graded_tabu && firstBestRows2.size() != 0) {
    unsigned ran = mersenne.next(firstBestRows2.size());
    replace(firstBestVars2[ran], firstBestRows2[ran]);
    return;
  }

  if ((! _change || mersenne.next(100) < 90) && firstBestRows.size() != 0) {
    unsigned ran = mersenne.next(firstBestRows.size());
    replace(firstBestVars[ran], firstBestRows[ran]);
    return;
  }

  // const unsigned tupleEncode =
  //     uncoveredTuples.encode(base % uncoveredTuples.size());
  // const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
  // const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);

  bestRows.clear();
  bestRows2.clear();
  std::vector<unsigned> bestTuple, bestTuple2;

  std::vector<unsigned> changedVars;
  std::vector<unsigned> org_vars;

  for (size_t i = 0; i < uncoveredTuples.size(); ++i) {
    const unsigned tupleEncode = uncoveredTuples.encode((base + i) % uncoveredTuples.size());
    // std::cout << i << " " << tupleEncode << "\n";
    const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
    const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);

    for (unsigned lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
      org_vars.clear();
      changedVars.clear();
      std::vector<unsigned> &line = array[lineIndex];
      for (unsigned i = 0; i < tuple.size(); ++i) {
        if (line[columns[i]] != tuple[i]) {
          org_vars.push_back(line[columns[i]]);
          changedVars.push_back(tuple[i]);
        }
      }
      if (changedVars.size() == 0) {
        continue;
      }
      // Tabu
      bool isTabu = false;
      for (auto v : changedVars) {
        unsigned diffOption = specificationFile.getOptions().option(v);
        if (entryTabu.isTabu(Entry(lineIndex, diffOption))) {
          isTabu = true;
          break;
        }
      }
      if (isTabu) {
        continue;
      }
      // check constraint, before tmpScore or after it?
      bool need_to_check = false;
      for (auto v : changedVars) {
        const Options &options = specificationFile.getOptions();
        if (option_constrained_indicator[options.option(v)]) {
          need_to_check = true;
          break;
        }
      }

      // my check
      if (need_to_check &&
          !validater.valida_changes(array[lineIndex], changedVars)) {
        continue;
      }
      // greedy
      std::pair<long long, long long> tmpScore;
      if (changedVars.size() > 1) {
        tmpScore = multiVarScoreOfRow2(changedVars, lineIndex);
      } else {
        tmpScore = varScoreOfRow3(changedVars[0], lineIndex);
      }

      if (use_graded_tabu && (! lineTabu.isTabu(Entry(lineIndex, 0)))) {
        if (cmp(bestScore2, tmpScore) == -1) {
          bestScore2 = tmpScore;
          bestRows2.clear();
          bestRows2.push_back(lineIndex);
          bestTuple2.clear();
          bestTuple2.emplace_back(tupleEncode);
        }
        else if (cmp(bestScore2, tmpScore) == 0) {
          bestRows2.push_back(lineIndex);
          bestTuple2.emplace_back(tupleEncode);
        }
      }

      if (cmp(bestScore, tmpScore) == -1) {
        bestScore = tmpScore;
        bestRows.clear();
        bestRows.push_back(lineIndex);
        bestTuple.clear();
        bestTuple.emplace_back(tupleEncode);
      }
      else if (cmp(bestScore, tmpScore) == 0) {
        bestRows.push_back(lineIndex);
        bestTuple.emplace_back(tupleEncode);
      }
    }

    if (! try_every)
      break ;
  }

  if (use_graded_tabu && bestRows2.size() != 0) {
    unsigned ran = mersenne.next(bestRows2.size());
    unsigned lineIndex = bestRows2[ran];
    unsigned tupleEncode = bestTuple2[ran];
    const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
    const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);

    changedVars.clear();
    for (unsigned i = 0; i < tuple.size(); ++i) {
      if (array[lineIndex][columns[i]] != tuple[i]) {
        changedVars.push_back(tuple[i]);
      }
    }
    if (changedVars.size() > 1) {
      // multiVarReplace(changedVars, lineIndex);
      for (auto v : changedVars) {
        replace(v, lineIndex);
      }
    } else {
      replace(changedVars[0], lineIndex);
    }
    return;
  }

  // need to handle when bestRows.size() == 0
  if (bestRows.size() != 0) {
    // std::cout << "change_mutivar" << std::endl;
    unsigned ran = mersenne.next(bestRows.size());
    unsigned lineIndex = bestRows[ran];
    unsigned tupleEncode = bestTuple[ran];
    const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
    const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);

    changedVars.clear();
    for (unsigned i = 0; i < tuple.size(); ++i) {
      if (array[lineIndex][columns[i]] != tuple[i]) {
        changedVars.push_back(tuple[i]);
      }
    }
    if (changedVars.size() > 1) {
      // multiVarReplace(changedVars, lineIndex);
      for (auto v : changedVars) {
        replace(v, lineIndex);
      }
    } else {
      replace(changedVars[0], lineIndex);
    }
    return;
  }

  if (__forced_greedy_percent)
    force_patching(uncoveredTuples.encode(base));
  else
    replaceRow2(mersenne.next(array.size()), uncoveredTuples.encode(base));
    // WithdrawTuple(uncoveredTuples.encode(base));
}

std::pair<long long, long long>
CoveringArray::multiVarScoreOfRow2(const std::vector<unsigned> &sortedMultiVars,
                                   const unsigned lineIndex) {
  const Options &options = specificationFile.getOptions();
  std::vector<unsigned> &line = array[lineIndex];
  std::vector<unsigned> changedColumns;
  for (auto var : sortedMultiVars) {
    int c = options.option(var);
    changedColumns.push_back(c);
  }

  long long coverChangeCount = 0, coverChangeCount2 = 0;
  for (auto tupleEncode : uncoveredTuples) {
    const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
    const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);

    bool needChange = true;
    size_t i = 0, j = 0;
    while (i != columns.size() && j != changedColumns.size()) {
      if (columns[i] == changedColumns[j]) {
        if (tuple[i] != sortedMultiVars[j]) {
          needChange = false;
          break;
        }
        i++;
        j++;
      } else {
        if (columns[i] < changedColumns[j]) {
          if (tuple[i] != line[columns[i]]) {
            needChange = false;
            break;
          }
          i++;
        } else {
          j++;
        }
      }
    }
    for (; i != columns.size(); ++i) {
      if (tuple[i] != line[columns[i]]) {
        needChange = false;
        break;
      }
    }
    if (needChange) {
      coverChangeCount += use_weight ? oneCoveredTuples.getWeight(tupleEncode) : 1;
    }
  }

  if (use_score2) {
    for (auto tupleEncode : oneCoveredTuples) {
      if (use_potential_score2 && (! twoCoveredTuples.is_potential(tupleEncode)))
        continue ;

      const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
      const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);

      bool needChange = true;
      size_t i = 0, j = 0;
      while (i != columns.size() && j != changedColumns.size()) {
        if (columns[i] == changedColumns[j]) {
          if (tuple[i] != sortedMultiVars[j]) {
            needChange = false;
            break;
          }
          i++;
          j++;
        } else {
          if (columns[i] < changedColumns[j]) {
            if (tuple[i] != line[columns[i]]) {
              needChange = false;
              break;
            }
            i++;
          } else {
            j++;
          }
        }
      }
      for (; i != columns.size(); ++i) {
        if (tuple[i] != line[columns[i]]) {
          needChange = false;
          break;
        }
      }
      if (needChange) {
        coverChangeCount2 += use_weight ? twoCoveredTuples.getWeight(tupleEncode) : 1;
      }
    }
  }

  std::vector<long long> overlapCount(changedColumns.size() + 1, 0);
  std::vector<long long> overlapCount2(changedColumns.size() + 1, 0);
  for (auto cc : changedColumns) {
    for (auto &ecEntry : oneCoveredTuples.getECbyLineVar(lineIndex, line[cc])) {
      unsigned tupleEncode = ecEntry.encode;
      const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
      const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);
      bool needChange = true;
      for (size_t i = 0; i < columns.size(); ++i) {
        if (line[columns[i]] != tuple[i]) {
          needChange = false;
          break;
        }
      }
      if (needChange) {
        long long olc = 0;
        size_t j = 0, k = 0;
        while (j != columns.size() && k != changedColumns.size()) {
          if (columns[j] == changedColumns[k]) {
            olc++;
            j++;
            k++;
          } else {
            columns[j] < changedColumns[k] ? j++ : k++;
          }
        }
        overlapCount[olc] += use_weight ? oneCoveredTuples.getWeight(tupleEncode) : 1;
      }
    }

    if (! use_score2)
      continue ;

    for (auto &ecEntry : twoCoveredTuples.getECbyLineVar(lineIndex, line[cc])) {
      unsigned tupleEncode = ecEntry.encode;
      const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
      const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);
      bool needChange = true;
      for (size_t i = 0; i < columns.size(); ++i) {
        if (line[columns[i]] != tuple[i]) {
          needChange = false;
          break;
        }
      }
      if (needChange) {
        long long olc = 0;
        size_t j = 0, k = 0;
        while (j != columns.size() && k != changedColumns.size()) {
          if (columns[j] == changedColumns[k]) {
            olc++;
            j++;
            k++;
          } else {
            columns[j] < changedColumns[k] ? j++ : k++;
          }
        }
        overlapCount2[olc] += use_weight ? twoCoveredTuples.getWeight(tupleEncode) : 1;
      }
    }

  }
  for (size_t olc = 1; olc < overlapCount.size(); ++olc) {
    if (! score_every_cell)
      coverChangeCount -= overlapCount[olc] / olc;
    else
      coverChangeCount -= overlapCount[olc] / olc / olc;
    if (use_score2) {
      if (! score_every_cell)
        coverChangeCount2 -= overlapCount2[olc] / olc;
      else
        coverChangeCount2 -= overlapCount2[olc] / olc / olc;
    }
  }
  return {coverChangeCount, coverChangeCount2};
}

long long
CoveringArray::multiVarRow(const std::vector<unsigned> &sortedMultiVars,
                           const unsigned lineIndex, const bool change) {
  const Options &options = specificationFile.getOptions();
  const unsigned strength = specificationFile.getStrenth();
  long long score = 0;

  std::vector<unsigned> varColumns;
  for (auto var : sortedMultiVars) {
    varColumns.push_back(options.option(var));
  }
  if (change) {
    lineTabu.insert(Entry(lineIndex, 0));
    for (auto column : varColumns) {
      entryTabu.insert(Entry(lineIndex, column));
    }
    lineTabu.insert(Entry(lineIndex, 0));
  }
  std::vector<unsigned> &line = array[lineIndex];

  // must from the end to the begining
  for (unsigned i = sortedMultiVars.size(); i--;) {
    std::swap(line[line.size() - sortedMultiVars.size() + i],
              line[varColumns[i]]);
  }

  std::vector<unsigned> tmpSortedColumns(strength);
  std::vector<unsigned> tmpSortedTupleToCover(strength);
  std::vector<unsigned> tmpSortedTupleToUncover(strength);
  unsigned tmpToCoverEncode;
  unsigned tmpToUncoverEncode;

  if (sortedMultiVars.size() >= strength) {
    for (std::vector<unsigned> changedColums = combinadic.begin(strength);
         changedColums[strength - 1] < sortedMultiVars.size();
         combinadic.next(changedColums)) {
      for (unsigned i = 0; i < strength; ++i) {
        tmpSortedTupleToCover[i] = sortedMultiVars[changedColums[i]];
        tmpSortedTupleToUncover[i] =
            line[line.size() - sortedMultiVars.size() + changedColums[i]];
      }
      std::sort(tmpSortedTupleToCover.begin(), tmpSortedTupleToCover.end());
      std::sort(tmpSortedTupleToUncover.begin(), tmpSortedTupleToUncover.end());
      for (unsigned i = 0; i < strength; ++i) {
        tmpSortedColumns[i] = options.option(tmpSortedTupleToCover[i]);
      }
      tmpToCoverEncode =
          coverage.encode(tmpSortedColumns, tmpSortedTupleToCover);
      tmpToUncoverEncode =
          coverage.encode(tmpSortedColumns, tmpSortedTupleToUncover);
      if (change) {
        cover(tmpToCoverEncode, lineIndex);
        uncover(tmpToUncoverEncode, lineIndex);
      } else {
        if (coverage.coverCount(tmpToCoverEncode) == 0) {
          ++score;
        }
        if (coverage.coverCount(tmpToUncoverEncode) == 1) {
          --score;
        }
      }
    }
  }

  for (unsigned curRelevantCount = 1,
                maxRelevantCount = std::min(
                    strength - 1,
                    (const unsigned)(line.size() - sortedMultiVars.size()));
       curRelevantCount <= maxRelevantCount; ++curRelevantCount) {
    for (std::vector<unsigned> relevantColumns =
             combinadic.begin(curRelevantCount);
         relevantColumns[curRelevantCount - 1] <
         line.size() - sortedMultiVars.size();
         combinadic.next(relevantColumns)) {
      for (std::vector<unsigned> changedColums =
               combinadic.begin(strength - curRelevantCount);
           changedColums[strength - curRelevantCount - 1] <
           sortedMultiVars.size();
           combinadic.next(changedColums)) {

        for (unsigned i = 0; i < curRelevantCount; ++i) {
          tmpSortedTupleToCover[i] = tmpSortedTupleToUncover[i] =
              line[relevantColumns[i]];
        }
        for (unsigned i = 0; i < strength - curRelevantCount; ++i) {
          tmpSortedTupleToCover[curRelevantCount + i] =
              sortedMultiVars[changedColums[i]];
          tmpSortedTupleToUncover[curRelevantCount + i] =
              line[line.size() - sortedMultiVars.size() + changedColums[i]];
        }
        std::sort(tmpSortedTupleToCover.begin(), tmpSortedTupleToCover.end());
        std::sort(tmpSortedTupleToUncover.begin(),
                  tmpSortedTupleToUncover.end());
        for (unsigned i = 0; i < strength; ++i) {
          tmpSortedColumns[i] = options.option(tmpSortedTupleToCover[i]);
        }
        tmpToCoverEncode =
            coverage.encode(tmpSortedColumns, tmpSortedTupleToCover);
        tmpToUncoverEncode =
            coverage.encode(tmpSortedColumns, tmpSortedTupleToUncover);
        if (change) {
          cover(tmpToCoverEncode, lineIndex);
          uncover(tmpToUncoverEncode, lineIndex);
        } else {
          if (coverage.coverCount(tmpToCoverEncode) == 0) {
            ++score;
          }
          if (coverage.coverCount(tmpToUncoverEncode) == 1) {
            --score;
          }
        }
      }
    }
  }
  // must from the begining to the end
  for (unsigned i = 0; i < sortedMultiVars.size(); ++i) {
    std::swap(line[line.size() - sortedMultiVars.size() + i],
              line[varColumns[i]]);
  }

  if (change) {
    for (unsigned i = 0; i < sortedMultiVars.size(); ++i) {
      line[varColumns[i]] = sortedMultiVars[i];
    }
  }
  return score;
}

long long
CoveringArray::multiVarScoreOfRow(const std::vector<unsigned> &sortedMultiVars,
                                  const unsigned lineIndex) {
  return multiVarRow(sortedMultiVars, lineIndex, false);
}
void CoveringArray::multiVarReplace(
    const std::vector<unsigned> &sortedMultiVars, const unsigned lineIndex) {
  const Options &options = specificationFile.getOptions();
  std::vector<unsigned> org_vars;
  for (auto var : sortedMultiVars) {
    auto option = options.option(var);
    org_vars.push_back(array[lineIndex][option]);
  }
  multiVarRow(sortedMultiVars, lineIndex, true);
}

long long CoveringArray::varScoreOfRow(const unsigned var,
                                       const unsigned lineIndex) {
  const Options &options = specificationFile.getOptions();
  const unsigned strength = specificationFile.getStrenth();
  std::vector<unsigned> &line = array[lineIndex];
  const unsigned varOption = options.option(var);
  if (line[varOption] == var) {
    return 0;
  }
  std::swap(line[line.size() - 1], line[varOption]);

  long long coverChangeCount = 0;
  std::vector<unsigned> tmpSortedColumns(strength);
  std::vector<unsigned> tmpSortedTupleToCover(strength);
  std::vector<unsigned> tmpSortedTupleToUncover(strength);
  for (std::vector<unsigned> columns = combinadic.begin(strength - 1);
       columns[strength - 2] < line.size() - 1; combinadic.next(columns)) {
    for (unsigned i = 0; i < columns.size(); ++i) {
      tmpSortedTupleToUncover[i] = tmpSortedTupleToCover[i] = line[columns[i]];
    }
    tmpSortedTupleToCover[strength - 1] = var;
    tmpSortedTupleToUncover[strength - 1] = line[line.size() - 1];
    std::sort(tmpSortedTupleToCover.begin(), tmpSortedTupleToCover.end());
    std::sort(tmpSortedTupleToUncover.begin(), tmpSortedTupleToUncover.end());
    for (unsigned i = 0; i < tmpSortedTupleToCover.size(); ++i) {
      tmpSortedColumns[i] = options.option(tmpSortedTupleToCover[i]);
    }
    unsigned tmpTupleToCoverEncode =
        coverage.encode(tmpSortedColumns, tmpSortedTupleToCover);
    unsigned tmpTupleToUncoverEncode =
        coverage.encode(tmpSortedColumns, tmpSortedTupleToUncover);
    if (coverage.coverCount(tmpTupleToCoverEncode) == 0) {
      coverChangeCount++;
    }
    if (coverage.coverCount(tmpTupleToUncoverEncode) == 1) {
      coverChangeCount--;
    }
  }

  std::swap(line[line.size() - 1], line[varOption]);

  return coverChangeCount;
}

long long CoveringArray::varScoreOfRow2(const unsigned var,
                                        const unsigned lineIndex) {
  const Options &options = specificationFile.getOptions();
  std::vector<unsigned> &line = array[lineIndex];
  const unsigned varOption = options.option(var);
  if (line[varOption] == var) {
    return 0;
  }
  long long coverChangeCount = 0;
  for (auto tupleEncode : uncoveredTuples) {
    const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
    const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);
    bool match = false;
    bool needChange = true;
    for (size_t i = 0; i < columns.size(); ++i) {
      if (columns[i] == varOption) {
        match = true;
        if (tuple[i] != var) {
          needChange = false;
          break;
        }
      } else if (line[columns[i]] != tuple[i]) {
        needChange = false;
        break;
      }
    }
    if (match && needChange) {
      coverChangeCount++;
    }
  }

  for (auto &ecEntry :
       oneCoveredTuples.getECbyLineVar(lineIndex, line[varOption])) {
    unsigned tupleEncode = ecEntry.encode;
    const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
    const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);
    bool needChange = true;
    for (size_t i = 0; i < columns.size(); ++i) {
      if (line[columns[i]] != tuple[i]) {
        needChange = false;
        break;
      }
    }
    if (needChange) {
      coverChangeCount--;
    }
  }
  return coverChangeCount;
}

std::pair<long long, long long> CoveringArray::varScoreOfRow3(const unsigned var,
                                        const unsigned lineIndex) {
  const Options &options = specificationFile.getOptions();
  std::vector<unsigned> &line = array[lineIndex];
  const unsigned varOption = options.option(var);
  if (line[varOption] == var) {
    return {0, 0};
  }
  long long coverChangeCount = 0, coverChangeCount2 = 0;
  for (auto tupleEncode : uncoveredTuples) {
    const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
    const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);
    bool match = false;
    bool needChange = true;
    for (size_t i = 0; i < columns.size(); ++i) {
      if (columns[i] == varOption) {
        match = true;
        if (tuple[i] != var) {
          needChange = false;
          break;
        }
      } else if (line[columns[i]] != tuple[i]) {
        needChange = false;
        break;
      }
    }
    if (match && needChange) {
      coverChangeCount += use_weight ? oneCoveredTuples.getWeight(tupleEncode) : 1;
    }
  }

  if (use_score2) {
    for (auto tupleEncode : oneCoveredTuples) {
      if (use_potential_score2 && (! twoCoveredTuples.is_potential(tupleEncode)))
        continue ;
      const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
      const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);
      bool match = false;
      bool needChange = true;
      for (size_t i = 0; i < columns.size(); ++i) {
        if (columns[i] == varOption) {
          match = true;
          if (tuple[i] != var) {
            needChange = false;
            break;
          }
        } else if (line[columns[i]] != tuple[i]) {
          needChange = false;
          break;
        }
      }
      if (match && needChange) {
        coverChangeCount2 += use_weight ? twoCoveredTuples.getWeight(tupleEncode) : 1;
      }
    }
  }

  if (! use_weight) {
    long long score1 = coverChangeCount -
          oneCoveredTuples.getECbyLineVar(lineIndex, line[varOption]).size();
    
    long long score2 = 0;
    if (use_score2)
      score2 = coverChangeCount2 - twoCoveredTuples.getECbyLineVar(lineIndex, line[varOption]).size();
    
    return {score1, score2};
    // return score1 + 0.1 * score2;
  }
  else {
    
    // long long sumWeight = 0;
    // const std::vector<ECEntry> vec = oneCoveredTuples.getECbyLineVar(lineIndex, line[varOption]);
    // for (auto item  : vec)
    //     sumWeight += oneCoveredTuples.getWeight(item.encode);

    // if (sumWeight != oneCoveredTuples.getECWeightbyLineVar(lineIndex, line[varOption])) {
    //     std::cout << sumWeight << " " << oneCoveredTuples.getECWeightbyLineVar(lineIndex, line[varOption]) << "\n";
    //     puts("Error Meow !");
    //     exit(0);
    // }
    long long score1 = coverChangeCount -
          oneCoveredTuples.getECWeightbyLineVar(lineIndex, line[varOption]);
    long long score2 = 0;
    if (use_score2)
      score2 = coverChangeCount2 - twoCoveredTuples.getECWeightbyLineVar(lineIndex, line[varOption]);

    return {score1, score2};
    // return score1 + 0.1 * score2;
  }
}

bool CoveringArray::checkCovered(unsigned encode) {
  for (auto ec : uncoveredTuples) {
    if (ec == encode) {
      return true;
    }
  }
  return false;
}

// quite similar to varScoreOfRow function
void CoveringArray::replace(const unsigned var, const unsigned lineIndex) {

  // puts("replace");

  const Options &options = specificationFile.getOptions();
  const unsigned strength = specificationFile.getStrenth();
  std::vector<unsigned> &line = array[lineIndex];
  const unsigned varOption = options.option(var);

  entryTabu.insert(Entry(lineIndex, varOption));
  lineTabu.insert(Entry(lineIndex, 0));

  if (line[varOption] == var) {
    return;
  }
  std::swap(line[line.size() - 1], line[varOption]);

  std::vector<unsigned> tmpSortedColumns(strength);
  std::vector<unsigned> tmpSortedTupleToCover(strength);
  std::vector<unsigned> tmpSortedTupleToUncover(strength);
  for (std::vector<unsigned> columns = combinadic.begin(strength - 1);
       columns[strength - 2] < line.size() - 1; combinadic.next(columns)) {
    for (unsigned i = 0; i < columns.size(); ++i) {
      tmpSortedTupleToUncover[i] = tmpSortedTupleToCover[i] = line[columns[i]];
    }
    tmpSortedTupleToCover[strength - 1] = var;
    tmpSortedTupleToUncover[strength - 1] = line[line.size() - 1];
    std::sort(tmpSortedTupleToCover.begin(), tmpSortedTupleToCover.end());
    std::sort(tmpSortedTupleToUncover.begin(), tmpSortedTupleToUncover.end());
    for (unsigned i = 0; i < tmpSortedTupleToCover.size(); ++i) {
      tmpSortedColumns[i] = options.option(tmpSortedTupleToCover[i]);
    }
    unsigned tmpTupleToCoverEncode =
        coverage.encode(tmpSortedColumns, tmpSortedTupleToCover);
    unsigned tmpTupleToUncoverEncode =
        coverage.encode(tmpSortedColumns, tmpSortedTupleToUncover);
    // need not check coverCount, cover(encode) will do this
    // puts("> <");
    cover(tmpTupleToCoverEncode, lineIndex);
    // puts("0.0");
    uncover(tmpTupleToUncoverEncode, lineIndex);
    // puts("0v0");
  }
  std::swap(line[line.size() - 1], line[varOption]);
  line[varOption] = var;

  check_score2(4);

  // puts("replace End");
}

void CoveringArray::cover(const unsigned encode, const unsigned oldLineIndex) {
  coverage.cover(encode);
  unsigned coverCount = coverage.coverCount(encode);
  if (coverCount == 1) {
    uncoveredTuples.pop(encode);
    oneCoveredTuples.push(encode, oldLineIndex, coverage.getTuple(encode));
  }
  if (coverCount == 2) {
    const std::vector<unsigned> &tuple = coverage.getTuple(encode);
    const std::vector<unsigned> &columns = coverage.getColumns(encode);
    for (size_t lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
      if (lineIndex == oldLineIndex) {
        continue;
      }
      auto &line = array[lineIndex];
      bool match = true;
      for (size_t i = 0; i < columns.size(); ++i) {
        if (tuple[i] != line[columns[i]]) {
          match = false;
          break;
        }
      }
      if (match) {
        oneCoveredTuples.pop(encode, lineIndex, tuple);
        if (use_score2)
          twoCoveredTuples.push(encode, lineIndex, tuple);
        break;
      }
    }
    if (use_score2)
      twoCoveredTuples.push(encode, oldLineIndex, tuple);
  }

  if (coverCount == 3 && use_score2) { // 2 -> 3
    // std::cout << oldLineIndex << "\n";
    const std::vector<unsigned> &tuple = coverage.getTuple(encode);
    const std::vector<unsigned> &columns = coverage.getColumns(encode);
    for (size_t lineIndex = 0, num = 0; lineIndex < array.size(); ++lineIndex) {
      if (lineIndex == oldLineIndex) {
        continue;
      }
      auto &line = array[lineIndex];
      bool match = true;
      for (size_t i = 0; i < columns.size(); ++i) {
        if (tuple[i] != line[columns[i]]) {
          match = false;
          break;
        }
      }
      if (match) {
        twoCoveredTuples.pop(encode, lineIndex, tuple);
        num ++;
        if (num == 2)
          break ;
      }
    }
  }
}

void CoveringArray::uncover(const unsigned encode,
                            const unsigned oldLineIndex) {

  coverage.uncover(encode);
  unsigned coverCount = coverage.coverCount(encode);
  if (coverCount == 0) { // 1 -> 0
    uncoveredTuples.push(encode);
    oneCoveredTuples.pop(encode, oldLineIndex, coverage.getTuple(encode));
  }
  if (coverCount == 1) { // 2 -> 1
    const std::vector<unsigned> &tuple = coverage.getTuple(encode);
    const std::vector<unsigned> &columns = coverage.getColumns(encode);
    for (size_t lineIndex = 0; lineIndex < array.size(); ++lineIndex) {
      if (lineIndex == oldLineIndex) {
        continue;
      }
      auto &line = array[lineIndex];
      bool match = true;
      for (size_t i = 0; i < columns.size(); ++i) {
        if (tuple[i] != line[columns[i]]) {
          match = false;
          break;
        }
      }
      if (match) {
        oneCoveredTuples.push(encode, lineIndex, tuple);
        if (use_score2)
          twoCoveredTuples.pop(encode, lineIndex, tuple);
        break ;
      }
    }
    if (use_score2)
      twoCoveredTuples.pop(encode, oldLineIndex, tuple);
  }
  if (coverCount == 2 && use_score2) { // 3 -> 2
    const std::vector<unsigned> &tuple = coverage.getTuple(encode);
    const std::vector<unsigned> &columns = coverage.getColumns(encode);
    for (size_t lineIndex = 0, num = 0; lineIndex < array.size(); ++lineIndex) {
      if (lineIndex == oldLineIndex) {
        continue;
      }
      auto &line = array[lineIndex];
      bool match = true;
      for (size_t i = 0; i < columns.size(); ++i) {
        if (tuple[i] != line[columns[i]]) {
          match = false;
          break;
        }
      }
      if (match) {
        twoCoveredTuples.push(encode, lineIndex, tuple);
        num ++;
        // std::cout << num << " " << lineIndex << " " << encode << "\n";
        if (num == 2)
          break ;
      }
    }
  }
}

void CoveringArray::tmpPrint() {
  // std::cout << (double)(clock() - clock_start) / CLOCKS_PER_SEC << '\t'
  //           << array.size() << '\t' << step << std::endl;

  const unsigned strength = specificationFile.getStrenth();

  std::cout << "\033[;32mc current " << strength << "-wise CA size: " 
                << array.size()
                << ", step #" << step
                << ", time #" << (double)(clock() - clock_start) / CLOCKS_PER_SEC
                << " \033[0m" << std::endl;
}

bool CoveringArray::verify(
    const std::vector<std::vector<unsigned>> &resultArray) {
  const unsigned strength = specificationFile.getStrenth();
  const Options &options = specificationFile.getOptions();
  Coverage tmpCoverage(specificationFile);
  tmpCoverage.initialize(satSolver, option_constrained_indicator);

  std::vector<unsigned> tuple(strength);
  unsigned lineIndex = 0;
  for (auto &line : resultArray) {
    for (unsigned column = 0; column < line.size(); ++column) {
      if (line[column] < options.firstSymbol(column) ||
          line[column] > options.lastSymbol(column)) {
        std::cerr << "error line: " << lineIndex;
        std::cerr << " option: " << column << std::endl;
        std::cerr << "should be " << options.firstSymbol(column)
                  << " <= var <= " << options.lastSymbol(column) << std::endl;
        for (auto var : line) {
          std::cerr << var << ' ';
        }
        std::cerr << std::endl;
        return false;
      }
    }
    for (std::vector<unsigned> columns = combinadic.begin(strength);
         columns[strength - 1] < line.size(); combinadic.next(columns)) {
      for (unsigned i = 0; i < strength; ++i) {
        tuple[i] = line[columns[i]];
      }
      unsigned encode = tmpCoverage.encode(columns, tuple);
      if (tmpCoverage.coverCount(encode) < 0) {
        for (unsigned i = 0; i < strength; ++i)
          std::cerr << tuple[i] << " ";
        std::cerr << "\n";
        std::cerr << "violate constraints" << std::endl;
        return false;
      }
      tmpCoverage.cover(encode);
    }
    ++lineIndex;
  }
  return tmpCoverage.allIsCovered();
}
