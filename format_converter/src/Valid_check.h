#ifndef VALID_CHECK_H_ONAP2ASG
#define VALID_CHECK_H_ONAP2ASG

#include "SpecificationFile.h"
#include <algorithm>
#include <iostream>
#include <vector>

namespace Valid {

#define NVISIBLE

class Literal {
private:
    bool neg;
    int var;

public:
    Literal() {}
    Literal(bool is_negative, int v) : neg(is_negative), var(v) {}

    void assign(bool is_negative, int v) {
        neg = is_negative;
        var = v;
    }
    bool is_negative() { return neg; }
    bool is_negative() const { return neg; }
    unsigned int variable() { return var; }
    unsigned int variable() const { return var; }
    int operator<(const Literal &other) const { return var < other.var; }
    bool operator==(const Literal &other) const {
        return neg == other.neg && var == other.var;
    }

#ifndef NVISIBLE
    void print() const {
        if (neg)
            std::cout << '-';
        std::cout << var;
    }
#endif
};

typedef std::vector<Literal> Clause;
typedef std::vector<Clause> Formula;

} // namespace Valid
#endif /* end of include guard: VALID_CHECK_H_ONAP2ASG */
