#ifndef LOCALSEARCH_OPTIMIZER_INCLUDE_H
#define LOCALSEARCH_OPTIMIZER_INCLUDE_H
#endif

#include "Expandor.h"

#include "../minisat_ext/BlackBoxSolver.h"
#include "../minisat_ext/Ext.h"
#include "MyList.h"

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

class Optimizer {

public:
    Optimizer (const Expandor &expandor, const Argument &params, bool whole);
    ~Optimizer ();

    void search ();
    void SaveTestcaseSet ();
    void Free ();
    vector<vector<int> > get_final_tc ();

private:

    int strength, seed = 1;
    std::mt19937_64 gen;
    int stop_length = 10000;

    bool use_sampling = false;
    int sampling_num = 100;
    bool use_group, use_weight;
    int *weight;

    string output_testcase_path;
    int nvar, nclauses, group_num;
    vector<vector<int> > clauses;
    vector<int> group_id;
    vector<vector<int> > member;
    
    vector<vector<int> > last_strengh_testcases;
    vector<vector<int> > testcases;
    int testcase_size;

    vector<t_tuple> tuples_U;
    int tuples_nums;
    vector<int> covered_times;

    vector<int> uncovered_tuples, covered_tuples;
    int covered_tuples_nums;

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

    vector<MyList> unique_covered_tuples;
    vector<MyList> tc_covered_tuples;
    vector<MyList> covered_testcases;

    vector<int> testcase_pos_to_idx;
    vector<int> testcase_idx_to_pos;
    int testcase_idx;

    vector<int> tuples_idx_to_pos;

    vector<LIST_node*> unique_node;
    vector<u_int64_t> sum_weight;

    void update_covered_testcases (int tpid);
    void UpdateInfo_remove_testcase (int tcid_idx);
    // u_int64_t update_unique_covered (int tcid);
    int get_which_remove ();
    void remove_testcase_greedily ();
    void remove_testcase_randomly ();

    void change_bit (int v, int ad, const vector<int>& tc, vector<int>& cur_clauses_cov);
    bool check_force_tuple (const vector<int>& tc, t_tuple tp, vector<int>& cur_clauses_cov);
    pair<u_int64_t, u_int64_t> get_gain_for_forcetestcase (int tcid, const vector<int>& tc2);
    pair<bool, pair<u_int64_t, u_int64_t> > get_gain_for_forcetuple (int tcid, t_tuple tp);
    pair<u_int64_t, u_int64_t> get_gain_for_forcetestcase_2 (int tcid, const vector<int>& tc2);
    pair<bool, pair<u_int64_t, u_int64_t> > get_gain_for_forcetuple_2 (int tcid, t_tuple tp);
    void forcetestcase (int tcid, const vector<int>& tc2);
    void forcetuple (int tcid, t_tuple tp);
    void backtrack_row (int tcid, t_tuple tp);
    bool gradient_descent ();

    void adaptive_adjustment (int strength, int nvar, int nclauses, int group_num);
    bool greedy_step_forced (t_tuple tp);
    void random_greedy_step (t_tuple tp, long long maxi);
    void random_greedy_step ();
    void random_greedy_step (long long maxi);

    bool check_for_flip(int tcid, int vid);
    void flip_bit(int tid, int vid);
    void random_step ();

    void remove_unnecessary_tc ();

    void SaveTestcaseSet(string result_path);
    
    bool is_covered (const vector<int> &tc, const t_tuple &t) {
        for (int i = 0; i < strength; i ++) {
            int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
            if (tc[pi] != vi)
                return false;
        }
        return true;
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

    vector<int> get_new_tc (const vector<int>& tc, const t_tuple &t) {
        vector<int> tc2 = tc;
        for (int i = 0; i < strength; i ++) {
            int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
            tc2[pi] = vi;
            if (use_group && (~group_id[pi]))
                for (int x : member[group_id[pi]])
                    if (x != pi)
                        tc2[x] = 0;
        }
        return tc2;
    }

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

    void set_taboo (int tcid, const vector<int>& tc2) {
        ++ greedy_limit;
        last_greedy_time[tcid] = greedy_limit;
        const vector<int>& tc = testcases[tcid];
        for (int i = 0; i < nvar; i ++)
            if (tc[i] != tc2[i])
                last_greedy_time_cell[tcid][i] = greedy_limit;
    }
};
