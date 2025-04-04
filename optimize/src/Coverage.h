#ifndef COVERAGE_H
#define COVERAGE_H

#include <algorithm>
#include <vector>

#include "Combinadic.h"
#include "PascalTriangle.h"
#include "SAT.h"
#include "SpecificationFile.h"

class Coverage {
public:
  Coverage(const SpecificationFile &specificationFile);
  void unconstrained_initialize();
  void set_invalid(const unsigned encode);
  void set_zero_invalid();
  void initialize(SATSolver &satSovler,
                  const std::vector<bool> &option_constrained_indicator);
  void initialize(const std::vector<std::vector<unsigned>> &array);
  int coverCount(const unsigned encode) const { return contents[encode]; }
  void cover(const unsigned encode) { ++contents[encode]; }
  void uncover(const unsigned encode) { --contents[encode]; }
  unsigned encode(const std::vector<unsigned> &sortedColumns,
                  const std::vector<unsigned> &sortedSubset);
  const std::vector<unsigned> &getColumns(const unsigned encode) const;
  const std::vector<unsigned> &getTuple(unsigned encode) const {
    return tuples[encode];
  }
  size_t tupleSize() { return contents.size(); }
  bool allIsCovered();
  unsigned tupleCount() const { return contents.size(); }
  void error() {
    for (unsigned i = 0; i < contents.size(); ++i) {
      if (contents[i] < 0) {
        std::cerr << "encode: " << i << "count: " << contents[i] << std::endl;
        abort();
      }
    }
  }
  void print();

private:
  const SpecificationFile &specificationFile;
  std::vector<unsigned> offsets;
  std::vector<std::vector<unsigned>> columns;
  std::vector<int> contents;
  std::vector<std::vector<unsigned>> tuples;
};

#endif /* end of include guard: COVERAGE_H */
