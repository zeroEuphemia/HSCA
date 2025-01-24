#ifndef LOCALSEARCH_OPTIMIZER_INCLUDE_H
#define LOCALSEARCH_OPTIMIZER_INCLUDE_H
#endif

#include "Optimizer.h"

#include "../minisat_ext/BlackBoxSolver.h"
#include "../minisat_ext/Ext.h"
#include "MyList_uint64.h"

#include <vector>
#include <numeric>
#include <random>
#include <utility>
#include <algorithm>
#include <set>
#include <map>
#include <iostream>
#include <cstring>
#include <unistd.h>
#include <chrono>

using std::vector;
using std::pair;
using std::set;
using std::map;
using std::string;

struct t_combnum {
    int idx[mx_strength + 1];
    void print (int sz) {
        for (int i = 0; i < sz; i ++)
            std::cout << idx[i] << " ";
        std::cout << "\n";
    }
};

class GlobalOptimizer {

public:
    GlobalOptimizer (const Expandor &expandor, const vector<vector<int> >& init_tc, const Argument &params);
    ~GlobalOptimizer ();

    void search ();
    void SaveTestcaseSet ();

private:

    int strength, seed = 1;
    std::mt19937_64 gen;
    int stop_length = 50000;

    bool use_group, use_weight;
    int *weight;

    string output_testcase_path;
    int nvar, mvar, nclauses, group_num;
    vector<vector<int> > clauses;
    vector<int> group_id, group_id2, group_to_index;
    vector<vector<int> > member;

    vector<vector<int> > testcases;
    int testcase_size;

    int *covered_times;

    // vector<u_int64_t> uncovered_tuples;
    vector<t_tuple> uncovered_tuples;

    vector<vector<int> > pos_in_cls;
    vector<vector<int> > neg_in_cls;
    vector<vector<int> > clauses_cov;

    int greedy_limit;
    int testcase_taboo = 4;
    int __forced_greedy_percent = 90;
    vector<int> last_greedy_time;
    vector<vector<int> > last_greedy_time_cell;
    int taboo_method = 1; // 1: cell 2: tc

    CDCLSolver:: Solver *cdcl_solver;

    vector<vector<int> > bestArray;


    vector<int> testcase_pos_to_idx, testcase_idx_to_pos;
    int testcase_idx;

    // vector<vector<MyList_64> > unique_covered_tuples;
    vector<vector<vector<u_int64_t> > > unique_covered_tuples;
    int *unique_idx;
    vector<u_int64_t> all_unique_covered_tuples;
    vector<int> covered_tcid;
    vector<t_combnum> unique_pos;
    
    int *tmp_covered_tcid;


    u_int64_t **combnum, *sum_all_possible_;
    vector<vector<vector<u_int64_t> > > base_combnum;

    // t_tuple *tuples_U;
    vector<t_tuple> tuples_U;
    u_int64_t tuples_num;

    u_int64_t GetBase (int t, int n, t_combnum comb) {
        if (t == 2)
            return (2ll * n - comb.idx[0] - 1) * comb.idx[0] / 2 + comb.idx[1] - comb.idx[0] - 1;
        u_int64_t res = combnum[t][n] - combnum[t][n - comb.idx[0]];

        t_combnum c;
        for (int i = 0; i < t - 1; i ++)
            c.idx[i] = comb.idx[i + 1] - comb.idx[0] - 1;
        return res + GetBase(t - 1, n - comb.idx[0] - 1, c);
    }

    t_combnum BaseToComb (int t, int n, u_int64_t base) {
        
        if (t == 1) {
            t_combnum c;
            c.idx[0] = (int) base;
            return c;
        }
        /*
            C(n, t) - C(n, t)
            idx[0] = 0   C(n - 1, t - 1)

            C(n, t) - C(n - 1, t)
            idx[0] = 1   C(n - 2, t - 1)
            
            C(n, t) - C(n - 2, t)
            idx[0] = 2
        */
        int c0 = upper_bound(base_combnum[t][n].begin(), base_combnum[t][n].end(), base) - base_combnum[t][n].begin() - 1;
        u_int64_t base0 = base - base_combnum[t][n][c0];
        t_combnum c = BaseToComb (t - 1, n - c0 - 1, base0);
        for (int i = t - 1; i; i --)
            c.idx[i] = c.idx[i - 1] + c0 + 1;
        c.idx[0] = c0;
        return c;
    }

    // u_int64_t TupleToIndex (t_tuple t) {

    //     t_combnum comb;
    //     for (int i = 0; i < strength; i ++)
    //         comb.idx[i] = abs(t.v[i]) - 1;

    //     u_int64_t base1 = GetBase(strength, nvar, comb), base2 = 0;
    //     for (int i = 0; i < strength; i ++) {
    //         int v = t.v[i] > 0;
    //         base2 |= v << (strength - i - 1);
    //     }
    //     return base2 * combnum[strength][nvar] + base1;
    // }

    // t_tuple IndexToTuple (u_int64_t index) {
    //     u_int64_t value = index / combnum[strength][nvar];
    //     u_int64_t base = index % combnum[strength][nvar];
    //     t_combnum comb = BaseToComb(strength, nvar, base);
    //     t_tuple t;
    //     for (int i = strength - 1; i >= 0; i ++, value >>= 1) {
    //         int pi = comb.idx[i], vi = (value & 1);
    //         t.v[i] = vi ? (pi + 1) : (- pi - 1);
    //     }
    //     return t;
    // }

    u_int64_t TupleToIndex (t_tuple t) {
        
        t_combnum comb;
        vector<pair<int, int> > res;
        for (int i = 0; i < strength; i ++) {
            int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
            res.emplace_back(make_pair(group_id2[pi], t.v[i]));
        }
        sort(res.begin(), res.end());

        for (int i = 0; i < strength; i ++)
            comb.idx[i] = res[i].first;
        u_int64_t base = GetBase(strength, mvar, comb);
        
        u_int64_t offset = 0;
        for (int i = 0; i < strength; i ++) {
            int g = res[i].first;
            int pi = abs(res[i].second) - 1, vi = res[i].second > 0;
            int num = g < group_num ? member[g].size() : 2;
            if (~group_id[pi]) {
                for (vi = 0; vi < (int)member[g].size(); vi ++)
                    if (member[g][vi] == pi)
                        break ;
            }
            offset = offset * num + vi;
        }
        return sum_all_possible_[base] + offset;
    }

    t_tuple IndexToTuple (u_int64_t index) {
        return tuples_U[index];
    }

    t_combnum begin_combnum (int strength) {
        t_combnum comb;
        for (int i = 0; i < strength; i ++)
            comb.idx[i] = i;
        return comb;
    }

    t_combnum end_combnum (int strength, int nvar) {
        t_combnum comb;
        for (int i = 0; i < strength; i ++)
            comb.idx[i] = nvar - strength + i; // nvar - strength, ..., nvar - strength + strength - 1
        return comb;
    }

    t_tuple end_tuple (int strength, int nvar) {
        t_tuple t;
        for (int i = 0; i < strength; i ++)
            t.v[i] = nvar - strength + i + 1;
        return t;
    }

    bool next_value (int strength, vector<int>& values, const vector<int>& mx_values) {
        for (int i = strength - 1; i >= 0; i --) {
            if (values[i] < mx_values[i] - 1) {
                values[i] ++; return true;
            }
            values[i] = 0;
        }
        return false;
    }

    t_tuple get_tuple (int strength, const t_combnum& comb, const vector<int>& value) {
        
        vector<pair<int, int> > tp;
        for (int i = 0; i < strength; i ++) {
            int pi = comb.idx[i], vi = value[i];
            if (pi < group_num)
                pi = member[pi][vi], vi = 1;
            else
                pi = group_to_index[pi];
            
            tp.emplace_back(make_pair(pi, vi));
        }
        sort(tp.begin(), tp.end());

        t_tuple t;
        for (int i = 0; i < strength; i ++) {
            int pi = tp[i].first, vi = tp[i].second;
            t.v[i] = vi ? (pi + 1) : (- pi - 1);
        }
        return t;
    }

    void next_comb (t_combnum &comb, int strength) {
        int limit = strength - 1, ceiling = comb.idx[0];
        for (int i = 0; i < limit; i ++) {
            int entry = ceiling + 1;
            ceiling = comb.idx[i + 1];
            if (entry < ceiling) {
                comb.idx[i] = entry;
                return ;
            }
            comb.idx[i] = i;
        }
        comb.idx[limit] ++;
    }

    bool is_covered (const vector<int> &tc, const t_tuple &t) {
        for (int i = 0; i < strength; i ++) {
            int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
            if (tc[pi] != vi)
                return false;
        }
        return true;
    }

    bool is_covered_changeone (const vector<int> &tc, const t_tuple &t, int opt, int value) {
        for (int i = 0; i < strength; i ++) {
            int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
            if ((pi != opt) && (tc[pi] != vi))
                return false;
        }
        return tc[opt] == value;
    }

    bool is_covered_changetwo (const vector<int> &tc, const t_tuple &t, int opt1, int value1, int opt2, int value2) {
        for (int i = 0; i < strength; i ++) {
            int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
            if ((pi != opt1) && (pi != opt2) && (tc[pi] != vi))
                return false;
        }
        return (tc[opt1] == value1) && (tc[opt2] == value2);
    }

    pair<int, int> get_different_var (const vector<int> &tc, const t_tuple &t) {
        int num = 0, first_var = -1;
        for (int i = 0; i < strength; i ++) {
            int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
            if (tc[pi] != vi) {
                num ++;
                if (first_var == -1)
                    first_var = pi;
            }
        }
        return {num, first_var};
    }

    vector<int> get_different_vars (const vector<int> &tc, const t_tuple &t) {
        vector<int> vars;
        for (int i = 0; i < strength; i ++) {
            int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
            if (tc[pi] != vi)
                vars.emplace_back(pi);
        }
        return vars;
    }

    void updateinfo_remove_testcase (int tcid_idx);

    int get_which_remove ();
    void remove_testcase_greedily ();
    void remove_testcase_randomly ();

    void change_bit (int v, int ad, const vector<int>& tc, vector<int>& cur_clauses_cov);
    bool check_force_tuple (const vector<int>& tc, t_tuple tp, vector<int>& cur_clauses_cov);

    pair<u_int64_t, u_int64_t> get_gain_for_flip_one (int tcid, int bit);
    pair<u_int64_t, u_int64_t> get_gain_for_flip_more (int tcid, const vector<int>& bits);
    
    void forcetestcase (int tcid, const vector<int>& tc2);
    void force_flip_one (int tcid, int bit);
    void force_flip_more (int tcid, const vector<int>& bits);

    void backtrack_row (int tcid, t_tuple tp);

    void remove_unnecessary_tc ();
    void random_greedy_step (t_tuple tp, long long maxi);
    void random_greedy_step (long long maxi);
    void gradient_descent ();

    bool is_taboo (int tcid, const t_tuple &t) {
        
        if (taboo_method == 2)
            return greedy_limit - last_greedy_time[tcid] <= testcase_taboo;
        
        const vector<int>& tc = testcases[tcid];
        for (int i = 0; i < strength; i ++) {
            int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
            if (tc[pi] == vi)
                continue ;
            
            if (greedy_limit - last_greedy_time_cell[tcid][pi] <= testcase_taboo)
                return true;
            if (use_group && (~group_id[pi]))
                for (int x : member[group_id[pi]])
                    if (greedy_limit - last_greedy_time_cell[tcid][x] <= testcase_taboo)
                        return true;
        }
        return false;
    }

    void set_taboo (int tcid, const t_tuple &t) {

        ++ greedy_limit;
        last_greedy_time[tcid] = greedy_limit;
        const vector<int>& tc = testcases[tcid];
        for (int i = 0; i < strength; i ++) {
            int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
            if (tc[pi] != vi) {
                last_greedy_time_cell[tcid][pi] = greedy_limit;

                if (use_group && (~group_id[pi]))
                    for (int x : member[group_id[pi]])
                        last_greedy_time_cell[tcid][x] = greedy_limit;
            }
        }
    }

    void set_taboo (int tcid, int var) {
        ++ greedy_limit;
        last_greedy_time[tcid] = greedy_limit;

        last_greedy_time_cell[tcid][var] = greedy_limit;
        if (~group_id[var])
            for (int x : member[group_id[var]])
                last_greedy_time_cell[tcid][x] = greedy_limit;
    }

    void set_taboo (int tcid, const vector<int>& tc2) {
        ++ greedy_limit;
        last_greedy_time[tcid] = greedy_limit;
        const vector<int>& tc = testcases[tcid];
        for (int i = 0; i < nvar; i ++)
            if (tc[i] != tc2[i])
                last_greedy_time_cell[tcid][i] = greedy_limit;
    }

    void set_covered_times ();
    vector<thread> Threads;
    int thread_num = 32;

    void set_sum_all_possible_ ();

    void SaveTestcaseSet(string result_path);

    vector<u_int64_t> sum_weight;
    vector<vector<u_int64_t> > sum_weight_var;
    bool greedy_step_forced (t_tuple tp);
    pair<u_int64_t, u_int64_t> get_gain_for_forcetestcase (int tcid, const vector<int>& tc2);


    void check_unique (int state);
};
