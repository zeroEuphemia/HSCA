#ifndef LINEVARTUPLESET_H_MBOTX5KJ
#define LINEVARTUPLESET_H_MBOTX5KJ

#include <vector>

#include "Combinadic.h"
#include "Coverage.h"
#include "SpecificationFile.h"

struct ECEntry {
public:
  unsigned encode;
  size_t column_index;
};

struct pair4 {
  int line1, pos1, line2, pos2;
};

class LineVarTupleSet {
public:
  LineVarTupleSet(int use_wieght_, int num, int use_potential_score2_){
    use_weight = use_wieght_;
    Num = num;
    use_potential_score2 = use_potential_score2_;
  };
  void initialize(const SpecificationFile &specificationFile,
                  const unsigned array_size);
  void pop(const unsigned encode, const unsigned lineIndex,
           const std::vector<unsigned> &tuple);
  void push(const unsigned encode, const unsigned lineIndex,
            const std::vector<unsigned> &tuple);
  unsigned encode(const unsigned index) { return tupleSet[index]; }
  unsigned size() const { return tupleSet.size(); }
  std::vector<unsigned>::const_iterator begin() const {
    return tupleSet.begin();
  }
  std::vector<unsigned>::const_iterator end() const { return tupleSet.end(); }

  void pushOneCoveredTuple(const Coverage &coverage,
                           const std::vector<std::pair<int, int> > &coverByLineindex);
  void pushTwoCoveredTuple(const Coverage &coverage,
                           const std::vector<std::pair<int, int> > &coverByLineindex);

  void exchange_row(unsigned lineIndex1, unsigned lineIndex2) {
    lineVarTupleSet[lineIndex1].swap(lineVarTupleSet[lineIndex2]);
    std::swap(lineOneCoveredCount[lineIndex1], lineOneCoveredCount[lineIndex2]);

    if (use_weight) {
      sum_weight[lineIndex1].swap(sum_weight[lineIndex2]);
      std::swap(sum_weight_line[lineIndex1], sum_weight_line[lineIndex2]);
    }

    if (Num == 2)
      std::swap(line_idx_map[lineIndex1], line_idx_map[lineIndex2]);
  }
  void pop_back_row() {
    lineVarTupleSet.pop_back();
    lineOneCoveredCount.pop_back();

    if (use_weight) {
      sum_weight.pop_back();
      sum_weight_line.pop_back();
    }

    if (Num == 2)
      line_idx_map.pop_back();
  }

  void addLine(unsigned allSymbolCount) {
    lineVarTupleSet.resize(lineVarTupleSet.size() + 1);
    lineVarTupleSet.rbegin()->resize(allSymbolCount);
    lineOneCoveredCount.push_back(0);

    if (use_weight) {
      sum_weight.resize(sum_weight.size() + 1);
      sum_weight[sum_weight.size() - 1] = vector<u_int64_t>(allSymbolCount, 0);
      sum_weight_line.push_back(0);
    }
  }

  const std::vector<ECEntry> &getECbyLineVar(unsigned lineIndex, unsigned var) {
    return lineVarTupleSet[lineIndex][var];
  }
  u_int64_t getECWeightbyLineVar(unsigned lineIndex, unsigned var) {
    return sum_weight[lineIndex][var];
  }

  unsigned oneCoveredCount(unsigned lineIndex) {
    return lineOneCoveredCount[lineIndex];
  }
  u_int64_t oneCoveredCount_Weight(unsigned lineIndex) {
    return sum_weight_line[lineIndex];
  }

  void addWeight(unsigned encode) {
    weight[encode] ++;
  }
  void addWeight(unsigned encode, u_int64_t w) {
    weight[encode] += w;
  }
  u_int64_t getWeight(unsigned encode) {
    return weight[encode];
  }

  void add_potential_tuple(unsigned encode) {
    if (! use_potential_score2)
      puts("Error in add_potential_tuple !"), exit(0);
    
    if (! is_potential_tuple[encode]) {
      is_potential_tuple[encode] = 1;
    }
  }
  bool is_potential(unsigned encode) {
    if (! use_potential_score2)
      puts("Error in add_potential_tuple !"), exit(0);
    return is_potential_tuple[encode];
  }

  void set_init_weight(const Coverage &coverage, const std::vector<InputClause> &constrs,
                                      const SpecificationFile &specificationFile);

private:
  int Num;
  int use_potential_score2 = 0;

  std::vector<unsigned> tupleSet; // contents of tuple encode
  std::vector<std::vector<unsigned>::size_type>
      mapping; // encode -> index in tupleSet
  
  // std::vector<unsigned> potential_tuple_set;
  std::vector<bool> is_potential_tuple;
  

  std::vector<size_t> lineOneCoveredCount;

  std::vector<std::vector<std::vector<ECEntry>>> lineVarTupleSet;
  std::vector<std::vector<size_t>> varMapping;
  std::vector<std::vector<pair4> > varMapping2;

  int use_weight;
  std::vector<u_int64_t> weight;
  std::vector<u_int64_t> sum_weight_line;
  std::vector<std::vector<u_int64_t> > sum_weight;

  vector<int> line_idx_map;
};

#endif /* end of include guard: LINEVARTUPLESET_H_MBOTX5KJ */
