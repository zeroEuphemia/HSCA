#include "LineVarTupleSet.h"

void LineVarTupleSet::initialize(const SpecificationFile &specificationFile,
                                 const unsigned array_size) {

  const Options &options = specificationFile.getOptions();
  unsigned strength = specificationFile.getStrenth();

  unsigned MaxSize = 0;
  for (std::vector<unsigned> tuple = combinadic.begin(strength);
       tuple[strength - 1] < options.size(); combinadic.next(tuple)) {
    unsigned blockSize = 1;
    for (unsigned i = 0; i < strength; ++i) {
      blockSize *= options.symbolCount(tuple[i]);
    }
    MaxSize += blockSize;
  }
  mapping.resize(MaxSize);
  if (Num == 1)
    varMapping.resize(MaxSize, std::vector<size_t>(strength));
  if (Num == 2) {
    varMapping2 = std::vector<std::vector<pair4> > (MaxSize, std::vector<pair4>(strength, {-1, -1, -1, -1}));
    for (int i = 0; i < array_size; i ++)
      line_idx_map.emplace_back(i);
    
    if (use_potential_score2) {
      is_potential_tuple.resize(MaxSize);
    }
  }

  // std::cout << "Num = " << Num << ", use_weight = " << use_weight << "\n";
  if (use_weight) {
    weight = std::vector<u_int64_t>(MaxSize, 1);
    sum_weight = std::vector<std::vector<u_int64_t> >(array_size, vector<u_int64_t>(options.allSymbolCount(), 0));
    sum_weight_line = std::vector<u_int64_t>(array_size, 0);
  }

  lineVarTupleSet.resize(array_size);
  for (unsigned i = 0; i < array_size; ++i) {
    lineVarTupleSet[i].resize(options.allSymbolCount());
  }

  lineOneCoveredCount.resize(array_size, 0);
}

void LineVarTupleSet::pushOneCoveredTuple(
    const Coverage &coverage, const std::vector<std::pair<int, int> > &coverByLineindex) {
  for (unsigned encode = 0; encode < coverage.tupleCount(); ++encode) {
    if (coverage.coverCount(encode) == 1) {
      push(encode, coverByLineindex[encode].first, coverage.getTuple(encode));
    }
  }
}

void LineVarTupleSet::pushTwoCoveredTuple(
    const Coverage &coverage, const std::vector<std::pair<int, int> > &coverByLineindex) {
  for (unsigned encode = 0; encode < coverage.tupleCount(); ++encode) {
    if (use_potential_score2 && (! is_potential_tuple[encode]))
      continue ;
    if (coverage.coverCount(encode) == 2) {
      push(encode, coverByLineindex[encode].first, coverage.getTuple(encode));
      push(encode, coverByLineindex[encode].second, coverage.getTuple(encode));
    }
  }
}

void LineVarTupleSet::push(const unsigned encode, const unsigned lineIndex,
                           const std::vector<unsigned> &tuple) {
  
  if (Num == 2 && use_potential_score2 && (! is_potential(encode)))
    return ;

  mapping[encode] = tupleSet.size();
  tupleSet.push_back(encode);
  lineOneCoveredCount[lineIndex]++;
  if (use_weight)
    sum_weight_line[lineIndex] += weight[encode];
  
  for (size_t i = 0; i < tuple.size(); ++i) {
    unsigned var = tuple[i];
    if (Num == 1)
      varMapping[encode][i] = lineVarTupleSet[lineIndex][var].size();

    else if (Num == 2) {
      if (varMapping2[encode][i].line1 == -1) {
        varMapping2[encode][i].line1 = line_idx_map[lineIndex];
        varMapping2[encode][i].pos1 = lineVarTupleSet[lineIndex][var].size();
      }
      else if (varMapping2[encode][i].line2 == -1) {
        varMapping2[encode][i].line2 = line_idx_map[lineIndex];
        varMapping2[encode][i].pos2 = lineVarTupleSet[lineIndex][var].size();
      }
      else
        puts("1 Error > < !"), exit(0);
    }

    lineVarTupleSet[lineIndex][var].push_back({encode, i});
    if (use_weight)
      sum_weight[lineIndex][var] += weight[encode];
  }
}

void LineVarTupleSet::pop(const unsigned encode, const unsigned lineIndex,
                          const std::vector<unsigned> &tuple) {
  
  if (Num == 2 && use_potential_score2 && (! is_potential(encode)))
    return ;

  tupleSet[mapping[encode]] = tupleSet[tupleSet.size() - 1];
  mapping[tupleSet[tupleSet.size() - 1]] = mapping[encode];
  tupleSet.pop_back();

  lineOneCoveredCount[lineIndex]--;
  if (use_weight)
    sum_weight_line[lineIndex] -= weight[encode];
  for (size_t i = 0; i < tuple.size(); ++i) {
    unsigned var = tuple[i];
    if (use_weight)
      sum_weight[lineIndex][var] -= weight[encode];
    
    std::vector<ECEntry> &varTS = lineVarTupleSet[lineIndex][var];

    if (Num == 1) {
      varTS[varMapping[encode][i]] = *varTS.rbegin();
      varMapping[varTS.rbegin()->encode][varTS.rbegin()->column_index] =
        varMapping[encode][i];
      varTS.pop_back();
    }    
    
    else if (Num == 2) {

      // std::cout << lineIndex << " " << varMapping2[encode][i].line1 << " " << varMapping2[encode][i].line2 << " " << var << "\n";

      if (line_idx_map[lineIndex] == varMapping2[encode][i].line1) {

        varTS[varMapping2[encode][i].pos1] = *varTS.rbegin();

        if (varMapping2[varTS.rbegin()->encode][varTS.rbegin()->column_index].line1 == line_idx_map[lineIndex])
          varMapping2[varTS.rbegin()->encode][varTS.rbegin()->column_index].pos1 = varMapping2[encode][i].pos1;
        else if (varMapping2[varTS.rbegin()->encode][varTS.rbegin()->column_index].line2 == line_idx_map[lineIndex])
          varMapping2[varTS.rbegin()->encode][varTS.rbegin()->column_index].pos2 = varMapping2[encode][i].pos1;
        else
          puts("2 Error > < !"), exit(0);

        varTS.pop_back();
        varMapping2[encode][i].line1 = varMapping2[encode][i].pos1 = -1;
      }

      else if (line_idx_map[lineIndex] == varMapping2[encode][i].line2) {

        varTS[varMapping2[encode][i].pos2] = *varTS.rbegin();

        if (varMapping2[varTS.rbegin()->encode][varTS.rbegin()->column_index].line1 == line_idx_map[lineIndex])
          varMapping2[varTS.rbegin()->encode][varTS.rbegin()->column_index].pos1 = varMapping2[encode][i].pos2;
        else if (varMapping2[varTS.rbegin()->encode][varTS.rbegin()->column_index].line2 == line_idx_map[lineIndex])
          varMapping2[varTS.rbegin()->encode][varTS.rbegin()->column_index].pos2 = varMapping2[encode][i].pos2;
        else
          puts("3 Error > < !"), exit(0);

        varTS.pop_back();
        varMapping2[encode][i].line2 = varMapping2[encode][i].pos2 = -1;
      }

      else
        puts("4 Error > < !"), exit(0);
      
      // std::cout << lineIndex << " " << varMapping2[encode][i].line1 << " " << varMapping2[encode][i].line2 << " " << var << "\n";
      // puts("-----------------");
    }
    
  }
}


void LineVarTupleSet::set_init_weight(const Coverage &coverage, const std::vector<InputClause> &constrs,
                                      const SpecificationFile &specificationFile) {

  if (! use_weight)
    return ;

  const Options &options = specificationFile.getOptions();
  int nvar = options.allSymbolCount();
  std::vector<unsigned> var_weight(nvar, 0);

  for (InputClause cls : constrs)
    for (InputTerm lit : cls.literals)
      if (lit.isNegated())
        var_weight[lit.getVariable()] ++;

  // std::cout << "nvar = " << nvar << "\n";
  // for (int i = 0; i < nvar; i ++)
  //   std::cout << var_weight[i] << " ";
  // std::cout << "\n";

  for (unsigned encode = 0; encode < coverage.tupleCount(); ++encode) {
    const std::vector<unsigned> &tuple = coverage.getTuple(encode);
    u_int64_t w = 0;
    for (int i = 0; i < tuple.size(); i ++) {
      // w = std::max(w, (u_int64_t)var_weight[tuple[i]] * 10 / constrs.size()); // output
      // w += (u_int64_t)var_weight[tuple[i]] * 10 / constrs.size(); // output2
      // w = std::max(w, (u_int64_t)var_weight[tuple[i]]); // output4
      w += (u_int64_t)var_weight[tuple[i]];
    }
    addWeight(encode, w);
    // std::cout << weight[encode] << " ";
  }
}
