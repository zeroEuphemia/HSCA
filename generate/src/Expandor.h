#ifndef LOCALSEARCH_OPTIMIZER_INCLUDE_H
#define LOCALSEARCH_OPTIMIZER_INCLUDE_H
#endif
#include "./core/Solver.h"
#include "../minisat_ext/BlackBoxSolver.h"
#include "../minisat_ext/Ext.h"
#include "Argument.h"
#include "cnfinfo.h"
#include "MyBitSet.h"

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
#include <thread>

using std::vector;
using std::pair;
using std::set;
using std::map;
using std::string;
using std::thread;

const int mx_strength = 6;
struct t_tuple {
    int v[mx_strength + 1];
    t_tuple () { }
    t_tuple (int strength) { }
    t_tuple (const vector<int> &vec) {
        for (int i = 0; i < (int) vec.size(); i ++)
            v[i] = vec[i];
    }
    void print(int strength) {
        std::cout << "(";
        for (int i = 0; i < strength; i ++)
            std::cout << v[i] << ", ";
        std::cout << ")\n";
    }
};

class Expandor {
public:
    Expandor (const Argument &params);
    ~Expandor ();
    void Expand ();
    void SaveTestcaseSet ();
    void Free ();
    vector<vector<int> > get_final_tc ();
    
    int strength;

    string output_testcase_path;
    
    int seed, nvar, nclauses, group_num;
    vector<int> group_id;
    vector<vector<int> > member;
    
    vector<vector<int> > clauses;
    vector<vector<int> > pos_in_cls;
    vector<vector<int> > neg_in_cls;

    vector<t_tuple> tuples_U;
    vector<vector<int> > old_testcase;
    vector<vector<int> > new_testcase;

    bool need_all_valid_tuples = false;
    vector<t_tuple> vaild_tuples;

private:
    int turnCNF(int x) {
        int v = (x >> 1) + 1;
        return (x & 1) ? v : -v;
    }
    int turnIdx(int x, bool v) {
        return x << 1 | v;
    }

    string input_cnf_path;
	string reduced_cnf_path;
	string init_CA_file_path;

    int candidate_set_size;
    int use_cnf_reduction;

    bool use_invalid_expand = true;
    bool use_cache = true;

    vector<t_tuple> uncovered_tuples;
    MyBitSet *covered_last_strength_bitmap;
    MyBitSet *covered_now_strength_bitmap;

    vector<t_tuple> t_clauses;

    vector<int> count_each_var_uncovered[2];

    vector<vector<int> > candidate_testcase_set_;

    int uncovered_nums;

    std::mt19937_64 gen;

    vector<MyBitSet> CA;

    vector<bool> group_flag;
    vector<vector<bool> > group_flag_array;

    vector<pair<int, int> > gain;

    void set_covered_now_strength_bitmap (const MyBitSet &tc, int idx, int last, t_tuple tuple);
    void get_remaining_valid_tuples (int thread_id, int value, int idx, int last, int _ed, t_tuple tuple);
    void thread_work(int thread_id, int st_value, int ed_value, int st_dfs, int ed_dfs);
    bool check_part_invalid(t_tuple tp);
    bool check_clauses_invalid(t_tuple tp);

    void GenerateCandidateTestcaseSet();
    void get_gain_thread(int st, int ed);
    int get_gain(const vector<int> &testcase);
    int SelectTestcaseFromCandidateSetByTupleNum();
    int GenerateTestcase();
    void Update_t_TupleInfo(const vector<int> &testcase, bool setmap);
    void Update_t_TupleInfo(int st, const vector<int> &testcase, const vector<int> &sidx);
    void ReplenishTestCase();
    void GenerateCoveringArray();
    void SaveTestcaseSet(string result_path);

    u_int64_t **combnum, num_combination_all_possible_;    
    u_int64_t GetBase(int t, int n, const vector<int> &vec) {
        if (t == 2)
            return (2ll * n - vec[0] - 1) * vec[0] / 2 + vec[1] - vec[0] - 1;
        long long res = combnum[t][n] - combnum[t][n - vec[0]];
        vector<int> v(t - 1);
        for (int i = 0; i < t - 1; i ++)
            v[i] = vec[i + 1] - vec[0] - 1;
        return res + GetBase(t - 1, n - vec[0] - 1, v);
    }
    u_int64_t TupleToIndex(int strength, t_tuple t) {
        
        if(strength == 1) {
            int idx = abs(t.v[0]) - 1;
            return t.v[0] < 0 ? idx : nvar + idx;
        }

        vector<int> vec(strength);
        for (int i = 0; i < strength; i ++)
            vec[i] = abs(t.v[i]) - 1;
        u_int64_t base1 = GetBase(strength, nvar, vec);
        u_int64_t base2 = 0;
        for (int i = 0; i < strength; i ++) {
            int v = t.v[i] > 0;
            base2 |= v << (strength - i - 1);
        }
        return base2 * combnum[strength][nvar] + base1;
    }

    bool is_covered (const vector<int> &tc, const t_tuple &t) {
        for (int i = 0; i < strength; i ++) {
            int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
            if (tc[pi] != vi)
                return false;
        }
        return true;
    }

    int thread_num = 32;
    u_int64_t *all_tuples, *invalid_nums, *covered_nums;
    // CDCLSolver:: Solver *cdcl_solver;
    vector<CDCLSolver:: Solver *> cdcl_solver_array;
    vector<thread> Threads;
    vector<vector<t_tuple> > uncovered_tuples_array;
    vector<vector<t_tuple> > vaild_tuples_array;

    // ExtMinisat:: SamplingSolver *cdcl_sampler;
    vector<ExtMinisat:: SamplingSolver *> cdcl_sampler_array;
    void cdcl_sampler_thread_work(int thread_id, const vector<pair<int, int> > &prob, int st, int ed);
};