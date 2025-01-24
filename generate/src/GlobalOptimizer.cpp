#include "GlobalOptimizer.h"

#include <iostream>
#include <unistd.h>
#include <chrono>
#include <string>
#include <vector>
#include <queue>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <numeric>
#include <iterator>

using std::vector;
using std::string;
using std::pair;
using std::mt19937_64;

GlobalOptimizer:: GlobalOptimizer (const Expandor &expandor, const vector<vector<int> >& init_tc, const Argument &params) {

    /*
        基础参数
    */
    strength = expandor.strength;
    gen.seed(seed = expandor.seed);
    stop_length = params.stop_length;

    __forced_greedy_percent = params.forced_greedy_percent;
    taboo_method = params.taboo_method;
    use_group = params.use_group;
    use_weight = params.use_weight;
    std::cout << "forced_greedy_percent = " << __forced_greedy_percent << "\n";
    std::cout << "taboo_method = " << taboo_method << "\n";
    std::cout << "use_group = " << use_group << "\n";
    std::cout << "use_weight = " << use_weight << "\n";

    output_testcase_path = expandor.output_testcase_path;

    nvar = expandor.nvar;
    nclauses = expandor.nclauses;
    clauses = expandor.clauses;
    group_num = expandor.group_num;
    group_id = expandor.group_id;
    member = expandor.member;

    group_id2.resize(nvar);
    group_to_index.resize(nvar);
    mvar = group_num;
    for (int i = 0; i < nvar; i ++) {
        group_id2[i] = (group_id[i] == -1) ? (mvar ++) : group_id[i];
        if (group_id[i] == -1)
            group_to_index[group_id2[i]] = i;
        else 
            group_to_index[group_id2[i]] = group_id2[i];
    }

    testcases = init_tc;
    testcase_size = testcases.size();
    testcase_idx = testcase_size - 1;
    std::cout << "testcase_size = " << testcase_size << "\n";
    testcase_pos_to_idx = vector<int>(testcase_size, 0);
    testcase_idx_to_pos = vector<int>(testcase_size, 0);
    for (int i = 0; i < testcase_size; i ++) {
        testcase_pos_to_idx[i] = i;
        testcase_idx_to_pos[i] = i;
    }

    combnum = new u_int64_t *[mx_strength + 1];
    for (int i = 2; i <= mx_strength; i ++)
        combnum[i] = new u_int64_t[nvar + 1]{0};
    for (int i = 2; i <= nvar; i ++) {
        combnum[2][i] = 1ll * i * (i - 1) / 2;
        for (int j = 3; j <= mx_strength; j ++)
            combnum[j][i] = combnum[j][i - 1] + combnum[j - 1][i - 1];
    }

    base_combnum = vector<vector<vector<u_int64_t> > > (mx_strength + 1, 
        vector<vector<u_int64_t> >(nvar + 1, vector<u_int64_t>()));
    for (int i = 2; i <= mx_strength; i ++) {
        for (int j = 1; j <= nvar; j ++) {
            base_combnum[i][j].resize(j);
            for (int k = 0; k < j; k ++)
                base_combnum[i][j][k] = combnum[i][j] - combnum[i][j - k];
        }
    }

    std::cout << "mvar = " << mvar << "\n";
    std::cout << "group_num = " << group_num << ", " << member.size() << "\n";
    for (int i = 0; i < group_num; i ++)
        std::cout << member[i].size() << ", ";
    std::cout << "\n";

    u_int64_t combination_all_possible_ = combnum[strength][mvar];
    std::cout << "combination_all_possible_ = " << combination_all_possible_ << "\n";

    sum_all_possible_ = new u_int64_t[combination_all_possible_ + 2] {0};
    set_sum_all_possible_();
    sum_all_possible_[0] = 0;


    for (u_int64_t i = 1; i <= combination_all_possible_; i ++)
        sum_all_possible_[i] += sum_all_possible_[i - 1];
    std::cout << "set_sum_all_possible_ success\n";

    tuples_num = sum_all_possible_[combination_all_possible_];
    // tuples_U = new t_tuple [tuples_num + 1];
    tuples_U.resize(tuples_num + 2);
    std::cout << "tuples_num = " << tuples_num << "\n";

    cout << "set_covered_times unqiue_covered_tuples\n";
    set_covered_times();
    cout << "set_covered_times unqiue_covered_tuples end\n";

    pos_in_cls = expandor.pos_in_cls;
    neg_in_cls = expandor.neg_in_cls;

    clauses_cov = vector<vector<int> >(testcase_size, vector<int>(nclauses, 0));
    for (int i = 0; i < testcase_size; i ++) {
        const vector<int>& tc = testcases[i];
        vector<int>& cur_clauses_cov = clauses_cov[i];
        for (int j = 0; j < nvar; j ++) {
            const vector<int>& vec = (tc[j] ? pos_in_cls[j + 1]: neg_in_cls[j + 1]);
            for (int x: vec) ++ cur_clauses_cov[x];
        }
    }

    greedy_limit = 0;
    last_greedy_time = vector<int>(testcase_size, -testcase_taboo - 1);
    // last_greedy_time_cell = vector<vector<int> > (testcase_size, vector<int>(nvar, -testcase_taboo - 1));
    last_greedy_time_cell.resize(testcase_size + 1);
    for (int i = 0; i < testcase_size; i ++) {
        last_greedy_time_cell[i].resize(nvar + 1);
        for (int j = 0; j < nvar; j ++)
            last_greedy_time_cell[i][j] = -testcase_taboo - 1;
    }

    cdcl_solver = new CDCLSolver::Solver;
    cdcl_solver->read_clauses(nvar, clauses);

    if (use_weight) {
        weight = new int[tuples_num + 2];
        for (int i = 0; i < tuples_num; i ++)
            weight[i] = 1;
    }

    for (int i = 0; i < tuples_num; i ++) {
        t_tuple t = tuples_U[i];
        if (TupleToIndex(t) != i)
            std::cout << "! ERROR ! \n", exit(0);
    }

    cout << "Optimizer init success" << endl;
}

GlobalOptimizer:: ~GlobalOptimizer () {
}

void GlobalOptimizer:: set_sum_all_possible_ () {

    u_int64_t combination_all_possible_ = combnum[strength][mvar];

    for (t_combnum comb = begin_combnum(strength); 
        comb.idx[strength - 1] < mvar; next_comb(comb, strength)) {    
        u_int64_t base = GetBase(strength, mvar, comb) + 1;
        
        sum_all_possible_[base] = 1;
        for (int i = 0; i < strength; i ++) {
            int values_num = comb.idx[i] < group_num ? member[comb.idx[i]].size() : 2;
            sum_all_possible_[base] *= values_num;
        }
    }
}

void GlobalOptimizer:: set_covered_times () {
    
    unique_covered_tuples = vector<vector<vector<u_int64_t> > >(testcase_size,
        vector<vector<u_int64_t> >(nvar, vector<u_int64_t>()));
    
    covered_times = new int [tuples_num + 1] {0};

    sum_weight = vector<u_int64_t> (testcase_size, 0);
    sum_weight_var = vector<vector<u_int64_t> > (testcase_size, vector<u_int64_t>(nvar, 0));

    // set tuples_U
    vector<int> mx_values = vector<int>(strength, 0);
    for (t_combnum comb = begin_combnum(strength); 
        comb.idx[strength - 1] < mvar; next_comb(comb, strength)) {
        
        vector<int> values = vector<int>(strength, 0);
        for (int i = 0; i < strength; i ++)
            mx_values[i] = comb.idx[i] < group_num ? member[comb.idx[i]].size() : 2;
        
        do {
            t_tuple t = get_tuple(strength, comb, values);
            u_int64_t tpid = TupleToIndex(t);
            tuples_U[tpid] = t;
        } while(next_value(strength, values, mx_values));
    }

    // set covered_times
    tmp_covered_tcid = new int[tuples_num + 1];
    for (t_combnum comb = begin_combnum(strength); 
        comb.idx[strength - 1] < mvar; next_comb(comb, strength)) {
        
        vector<int> values = vector<int>(strength, 0);
        for (int i = 0; i < testcase_size; i ++) {
            
            const vector<int>& tc = testcases[i];
            for (int j = 0; j < strength; j ++) {
                int pj, vj;
                if (comb.idx[j] < group_num)
                    for (pj = comb.idx[j], vj = 0; tc[member[pj][vj]] == 0; vj ++);
                else
                    pj = group_to_index[comb.idx[j]], vj = tc[pj];
                values[j] = vj;
            }
            
            t_tuple t = get_tuple(strength, comb, values);
            u_int64_t tpid = TupleToIndex(t);
            if (tpid >= tuples_num)
                std::cout << "ERROR\n", exit(0);
            covered_times[tpid] ++;
            tmp_covered_tcid[tpid] = i;
        }
    }

    unique_idx = new int[tuples_num + 1] {0};
    for (int tpid = 0; tpid < tuples_num; tpid ++)
        if (covered_times[tpid] == 1) {
            t_tuple t = IndexToTuple(tpid);
            int tcid = tmp_covered_tcid[tpid];

            unique_idx[tpid] = all_unique_covered_tuples.size();
            all_unique_covered_tuples.emplace_back(tpid);
            covered_tcid.emplace_back(tcid);
            sum_weight[tcid] ++;
            t_combnum comb;
            for (int i = 0; i < strength; i ++) {
                int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
                
                comb.idx[i] = unique_covered_tuples[tcid][pi].size();
                unique_covered_tuples[tcid][pi].emplace_back(tpid);
                sum_weight_var[tcid][pi] ++;
            }
            unique_pos.emplace_back(comb);
        }
    
    delete tmp_covered_tcid;
}

void GlobalOptimizer:: check_unique (int state) {
    for (int i = 0; i < testcase_size; i ++)
        for (int j = 0; j < nvar; j ++) {
            int sz = unique_covered_tuples[i][j].size();
            for (int k = 0; k < sz; k ++) {
                u_int64_t tpid = unique_covered_tuples[i][j][k];
                t_tuple t = IndexToTuple(tpid);
                int upos = unique_idx[tpid];
                for (int r = 0; r < strength; r ++)
                    if (abs(t.v[r]) - 1 == j) {
                        if (unique_pos[upos].idx[r] != k)
                            std::cout << state << " ERROR > < !" << std::endl, exit(0);
                        break ;
                    }
            }
        }
}

void GlobalOptimizer:: updateinfo_remove_testcase (int tcid) {

    vector<u_int64_t> unique_id;
    vector<t_tuple> unique_tp;

    const vector<int> &tc = testcases[tcid];
    
    vector<int> columns;
    for (int i = 0; i < nvar; i ++)
        if (group_id[i] == -1 || tc[i] == 1)
            columns.emplace_back(i);
    int columns_num = columns.size();

    for (t_combnum comb = begin_combnum(strength); 
        comb.idx[strength - 1] < columns_num; next_comb(comb, strength)) {

        t_tuple t;
        for (int i = 0; i < strength; i ++) {
            int pi = columns[comb.idx[i]];
            int vi = tc[pi];
            t.v[i] = vi ? (pi + 1) : (- pi - 1);
        }
        u_int64_t tpid = TupleToIndex(t);
        covered_times[tpid] --;
        
        if (covered_times[tpid] == 0) {
            uncovered_tuples.emplace_back(t);

            int upos = unique_idx[tpid];
            int sz = all_unique_covered_tuples.size();
            if (upos != sz) {
                swap(all_unique_covered_tuples[upos], all_unique_covered_tuples[sz - 1]);
                swap(covered_tcid[upos], covered_tcid[sz - 1]);
                swap(unique_pos[upos], unique_pos[sz - 1]);
                unique_idx[all_unique_covered_tuples[upos]] = upos;
            }
            all_unique_covered_tuples.pop_back();
            covered_tcid.pop_back();
            unique_pos.pop_back();
        }

        else if (covered_times[tpid] == 1) {
            unique_id.emplace_back(tpid);
            unique_tp.emplace_back(t);
        }
    }

    int new_unique_num = unique_id.size();
    for (int i = 0; i < new_unique_num; i ++) {
        u_int64_t tpid = unique_id[i];
        t_tuple t = unique_tp[i];

        unique_idx[tpid] = all_unique_covered_tuples.size();
        all_unique_covered_tuples.emplace_back(tpid);

        for (int j = 0; j < testcase_size; j ++)
            if (is_covered(testcases[j], t) && j != tcid) {

                covered_tcid.emplace_back(testcase_pos_to_idx[j]);
                sum_weight[j] += use_weight ? weight[tpid] : 1;
                t_combnum comb;
                for (int k = 0; k < strength; k ++) {
                    int p = abs(t.v[k]) - 1;
                    // unique_covered_tuples[j][p].insert(tpid);
                    comb.idx[k] = unique_covered_tuples[j][p].size();
                    unique_covered_tuples[j][p].emplace_back(tpid);
                    sum_weight_var[j][p] += use_weight ? weight[tpid] : 1;
                }
                unique_pos.emplace_back(comb);
                break ;
            }
    }

    if (tcid != testcase_size - 1) {
                
        testcases[tcid] = testcases[testcase_size - 1];
        last_greedy_time[tcid] = last_greedy_time[testcase_size - 1];
        last_greedy_time_cell[tcid] = last_greedy_time_cell[testcase_size - 1];
        clauses_cov[tcid] = clauses_cov[testcase_size - 1];
        unique_covered_tuples[tcid] = unique_covered_tuples[testcase_size - 1];
        
        sum_weight[tcid] = sum_weight[testcase_size - 1];
        sum_weight_var[tcid] = sum_weight_var[testcase_size - 1];

        int idx = testcase_pos_to_idx[testcase_size - 1];
        testcase_idx_to_pos[idx] = tcid;
        testcase_pos_to_idx[tcid] = idx;
    }

    testcases.pop_back();
    last_greedy_time.pop_back();
    last_greedy_time_cell.pop_back();
    clauses_cov.pop_back();
    sum_weight.pop_back();
    sum_weight_var.pop_back();

    vector<vector<u_int64_t> >().swap(unique_covered_tuples[testcase_size - 1]);
    unique_covered_tuples.pop_back();

    testcase_size --;

    // check_unique (1);
}

int GlobalOptimizer:: get_which_remove () {

    u_int64_t mini = 0;
    vector<int> besttcs;

    for (int i = 0; i < testcase_size; i ++) {
        u_int64_t res = sum_weight[i];
        if (besttcs.empty() || res == mini)
            mini = res, besttcs.emplace_back(i);
        else if (res < mini) {
            mini = res;
            besttcs.clear(), besttcs.emplace_back(i);
        }
    }
    return besttcs[gen() % besttcs.size()];
}

void GlobalOptimizer:: remove_testcase_greedily () {
    int tcid = get_which_remove();
    updateinfo_remove_testcase(tcid);
}

void GlobalOptimizer:: remove_testcase_randomly () {
    updateinfo_remove_testcase(gen() % testcase_size);
}

void GlobalOptimizer:: change_bit (int v, int ad, const vector<int>& tc, vector<int>& cur_clauses_cov) {
    int vid = abs(v) - 1;
    int curbit = tc[vid], tt = v > 0;
    if (curbit != tt) {
        const vector<int>& var_cov_old = (curbit ? pos_in_cls[vid + 1]: neg_in_cls[vid + 1]);
        const vector<int>& var_cov_new = (tt ? pos_in_cls[vid + 1]: neg_in_cls[vid + 1]);
        for (int cid: var_cov_new) cur_clauses_cov[cid] += ad;
        for (int cid: var_cov_old) cur_clauses_cov[cid] -= ad;
    }
}

bool GlobalOptimizer:: check_force_tuple (const vector<int>& tc, t_tuple tp, vector<int>& cur_clauses_cov) {

    for (int i = 0; i < strength; i ++) {
        change_bit(tp.v[i], 1, tc, cur_clauses_cov);
        if (use_group) {
            int pi = abs(tp.v[i]) - 1, vi = tp.v[i] > 0;
            if (~group_id[pi])
                for (int x : member[group_id[pi]])
                    if (x != pi)
                        change_bit(-(x + 1), 1, tc, cur_clauses_cov);
        }
    }
    
    bool has0 = false;
    for (int i = 0; i < nclauses; ++i){
        if (cur_clauses_cov[i] == 0){
            has0 = true; break;
        }
    }
    for (int i = 0; i < strength; i ++) {
        change_bit(tp.v[i], -1, tc, cur_clauses_cov);
        if (use_group) {
            int pi = abs(tp.v[i]) - 1, vi = tp.v[i] > 0;
            if (~group_id[pi])
                for (int x : member[group_id[pi]])
                    if (x != pi)
                        change_bit(-(x + 1), -1, tc, cur_clauses_cov);
        }
    }
    return has0 ? false : true;
}

pair<u_int64_t, u_int64_t> GlobalOptimizer:: get_gain_for_flip_one (int tcid, int bit) {
    
    const vector<int>& tc = testcases[tcid];

    u_int64_t break_cnt = sum_weight_var[tcid][bit];
    int other = -1;
    if (~group_id[bit]) // value = 1
        for (int x : member[group_id[bit]])
            if (tc[x] == 1)
                other = x, break_cnt += sum_weight_var[tcid][x];

    u_int64_t gain_cnt = 0;
    for (t_tuple t : uncovered_tuples) {
        if (~group_id[bit]) {
            if (is_covered_changetwo(tc, t, bit, tc[bit] ^ 1, other, 0)) {
                u_int64_t tpid = TupleToIndex(t);
                gain_cnt += use_weight ? weight[tpid] : 1;
            }
        }
        else {
            if (is_covered_changeone(tc, t, bit, tc[bit] ^ 1)) {
                u_int64_t tpid = TupleToIndex(t);
                gain_cnt += use_weight ? weight[tpid] : 1;
            }
        }
    }
    return {break_cnt, gain_cnt};
}

pair<u_int64_t, u_int64_t> GlobalOptimizer:: get_gain_for_flip_more (int tcid, const vector<int>& bits) {

    // std::cout << "get_gain_for_flip_one begin" << std::endl;
    const vector<int>& tc = testcases[tcid];

    vector<int> tc2 = tc;
    for (int x : bits) {
        tc2[x] ^= 1;
        if (~group_id[x])
            for (int y : member[group_id[x]])
                if (x != y)
                    tc2[y] = 0;
    }

    u_int64_t break_cnt = 0, gain_cnt = 0;
    for (t_tuple t : uncovered_tuples)
        if (is_covered(tc2, t)) {
            u_int64_t tpid = TupleToIndex(t);
            gain_cnt += use_weight ? weight[tpid] : 1;
        }

    int different_num = 0;
    for (int i = 0; i < nvar; i ++)
        different_num += (tc[i] != tc2[i]);
    
    vector<u_int64_t> break_cnt_num(different_num + 2, 0);
    for (int bit : bits) {
        int opt = bit;
        if (~group_id[bit]) {// 0 -> 1
            for (int x : member[group_id[bit]])
                if (tc[x] == 1) {
                    opt = bit; break ;
                }
        }
        for (u_int64_t tpid : unique_covered_tuples[tcid][opt]) {
            t_tuple t = IndexToTuple(tpid);

            int num = get_different_var (tc, t).first;
            break_cnt_num[num] += use_weight ? weight[tpid] : 1;
        }
    }
    for (int i = 1; i <= different_num; i ++)
        break_cnt += break_cnt_num[i] / i;
    
    // std::cout << "get_gain_for_flip_one end" << std::endl;
    return {break_cnt, gain_cnt};
}

void GlobalOptimizer:: forcetestcase (int tcid, const vector<int>& tc2) {
    
    const vector<int>& tc = testcases[tcid];

    vector<u_int64_t> unique_id;
    vector<t_tuple> unique_tp;

    vector<int> columns;
    for (int i = 0; i < nvar; i ++)
        if (group_id[i] == -1 || tc[i] == 1)
            columns.emplace_back(i);
    int columns_num = columns.size();

    for (t_combnum comb = begin_combnum(strength); 
        comb.idx[strength - 1] < columns_num; next_comb(comb, strength)) {

        t_tuple t;
        for (int i = 0; i < strength; i ++) {
            int pi = columns[comb.idx[i]], vi = tc[pi];
            t.v[i] = vi ? (pi + 1) : (- pi - 1);
        }
        u_int64_t tpid = TupleToIndex(t);
        covered_times[tpid] --;

        if (covered_times[tpid] == 0) {
            uncovered_tuples.emplace_back(t);

            int upos = unique_idx[tpid];
            int sz = all_unique_covered_tuples.size();
            if (upos != sz) {
                swap(all_unique_covered_tuples[upos], all_unique_covered_tuples[sz - 1]);
                swap(covered_tcid[upos], covered_tcid[sz - 1]);
                swap(unique_pos[upos], unique_pos[sz - 1]);
                unique_idx[all_unique_covered_tuples[upos]] = upos;
            }
            all_unique_covered_tuples.pop_back();
            covered_tcid.pop_back();
            unique_pos.pop_back();
        }
        if (covered_times[tpid] == 1) {
            unique_id.emplace_back(tpid);
            unique_tp.emplace_back(t);
        }
    }

    int new_unique_num = unique_id.size();
    for (int i = 0; i < new_unique_num; i ++) {
        u_int64_t tpid = unique_id[i];
        t_tuple t = unique_tp[i];

        unique_idx[tpid] = all_unique_covered_tuples.size();
        all_unique_covered_tuples.emplace_back(tpid);

        for (int j = 0; j < testcase_size; j ++)
            if (is_covered(testcases[j], t) && j != tcid) {

                covered_tcid.emplace_back(testcase_pos_to_idx[j]);
                sum_weight[j] += use_weight ? weight[tpid] : 1;
                t_combnum comb;
                for (int k = 0; k < strength; k ++) {
                    int p = abs(t.v[k]) - 1;
                    // unique_covered_tuples[j][p].insert(tpid);
                    comb.idx[k] = unique_covered_tuples[j][p].size();
                    unique_covered_tuples[j][p].emplace_back(tpid);
                    sum_weight_var[j][p] += use_weight ? weight[tpid] : 1;
                }
                unique_pos.emplace_back(comb);
                break ;
            }
    }


    for (int i = 0; i < nvar; i ++) {
        unique_covered_tuples[tcid][i].clear();
        sum_weight_var[tcid][i] = 0;
    }
    sum_weight[tcid] = 0;
    testcases[tcid] = tc2;

    columns.clear();
    for (int i = 0; i < nvar; i ++)
        if (group_id[i] == -1 || tc2[i] == 1)
            columns.emplace_back(i);
    columns_num = columns.size();
    
    for (t_combnum comb = begin_combnum(strength); 
        comb.idx[strength - 1] < columns_num; next_comb(comb, strength)) {

        t_tuple t;
        for (int i = 0; i < strength; i ++) {
            int pi = columns[comb.idx[i]], vi = tc2[pi];
            t.v[i] = vi ? (pi + 1) : (- pi - 1);
        }
        u_int64_t tpid = TupleToIndex(t);

        covered_times[tpid] ++;
        if (covered_times[tpid] == 1) {

            unique_idx[tpid] = all_unique_covered_tuples.size();
            all_unique_covered_tuples.emplace_back(tpid);
            covered_tcid.emplace_back(testcase_pos_to_idx[tcid]);
            sum_weight[tcid] += use_weight ? weight[tpid] : 1;
            t_combnum comb;
            for (int i = 0; i < strength; i ++) {
                int p = abs(t.v[i]) - 1;
                // unique_covered_tuples[tcid][p].insert(tpid);
                comb.idx[i] = unique_covered_tuples[tcid][p].size();
                unique_covered_tuples[tcid][p].emplace_back(tpid);
                sum_weight_var[tcid][p] += use_weight ? weight[tpid] : 1;
            }
            unique_pos.emplace_back(comb);
        }
        else if (covered_times[tpid] == 2) {
            int unique_p = unique_idx[tpid];

            int tc_idx = covered_tcid[unique_p];
            int tc_pos = testcase_idx_to_pos[tc_idx];
            sum_weight[tc_pos] -= use_weight ? weight[tpid] : 1;
            t_combnum comb = unique_pos[unique_p];
            for (int i = 0; i < strength; i ++) {
                int ci = abs(t.v[i]) - 1;
                int pi = comb.idx[i];
                int sz = unique_covered_tuples[tc_pos][ci].size();

                if (pi >= sz) {
                    std::cout << "tpid =  " << tpid << ", ", t.print(strength);
                    comb.print(strength);
                    std::cout << "ci = " << ci << std::endl;
                    std::cout << pi << " " << sz << std::endl;
                    std::cout << "meow meow ! ERROR !" << std::endl, exit(0);
                }

                if (pi != sz - 1) {
                    swap(unique_covered_tuples[tc_pos][ci][pi], unique_covered_tuples[tc_pos][ci][sz - 1]);
                    
                    u_int64_t tt_pid = unique_covered_tuples[tc_pos][ci][pi];
                    t_tuple tt = IndexToTuple(tt_pid);
                    int tt_upos = unique_idx[tt_pid];
                    int k = 0;
                    for ( ; k < strength; k ++)
                        if (abs(tt.v[k]) - 1 == ci)
                            break ;
                    unique_pos[tt_upos].idx[k] = pi;
                }
                unique_covered_tuples[tc_pos][ci].pop_back();
                sum_weight_var[tc_pos][ci] -= use_weight ? weight[tpid] : 1;
            }

            int sz = all_unique_covered_tuples.size();
            if (unique_p != sz - 1) {
                swap(all_unique_covered_tuples[unique_p], all_unique_covered_tuples[sz - 1]);
                swap(covered_tcid[unique_p], covered_tcid[sz - 1]);
                swap(unique_pos[unique_p], unique_pos[sz - 1]);

                unique_idx[all_unique_covered_tuples[unique_p]] = unique_p;
            }
            all_unique_covered_tuples.pop_back();
            covered_tcid.pop_back();
            unique_pos.pop_back();
        }
    }
    vector<int> break_pos;
    int uncovered_cnt = uncovered_tuples.size();
    for (int i = 0; i < uncovered_cnt; i ++) {
        t_tuple t = uncovered_tuples[i];
        u_int64_t tpid = TupleToIndex(t);
        if (is_covered(tc2, t))
            break_pos.emplace_back(i);
    }
    int break_num = break_pos.size();
    for(int i = break_num - 1; i >= 0; i --) {
        int p = break_pos[i];
        if(p != uncovered_cnt)
            uncovered_tuples[p] = uncovered_tuples[uncovered_cnt - 1];
        uncovered_tuples.pop_back();
        uncovered_cnt --;
    }

    vector<int>& cur_clauses_cov = clauses_cov[tcid];
    cur_clauses_cov = vector<int>(nclauses, 0);
    for(int i = 0; i < nvar; i ++) {
        const vector<int>& var = (tc2[i] ? pos_in_cls[i + 1]: neg_in_cls[i + 1]);
        for (int cid: var) cur_clauses_cov[cid] ++;
    }

    // check_unique (2);
}

void GlobalOptimizer:: force_flip_one (int tcid, int bit) {

    vector<int>& tc = testcases[tcid];

    vector<u_int64_t> unique_id;
    vector<t_tuple> unique_tp;

    vector<int> columns;
    int g = group_id[bit] == -1 ? -2 : group_id[bit];
    for (int i = 0; i < nvar; i ++) {
        if (i == bit || group_id[i] == g)
            continue ;
        if (group_id[i] == -1 || tc[i] == 1)
            columns.emplace_back(i);
    }
    int columns_num = columns.size();

    int b = bit;
    if (~group_id[bit])
        for (int x : member[group_id[bit]])
            if (tc[x] == 1) {
                b = x; break ;
            }

    for (t_combnum comb = begin_combnum(strength - 1); 
        comb.idx[strength - 2] < columns_num; next_comb(comb, strength - 1)) {
        
        t_tuple t;
        int i = 0;
        for (; i < strength - 1 && columns[comb.idx[i]] < b; i ++) {
            int pi = columns[comb.idx[i]], vi = tc[pi];
            t.v[i] = vi ? (pi + 1) : (- pi - 1);
        }
        t.v[i] = tc[b] ? (b + 1) : (- b - 1);
        for (; i < strength - 1; i ++) {
            int pi = columns[comb.idx[i]], vi = tc[pi];
            t.v[i + 1] = vi ? (pi + 1) : (- pi - 1);
        }

        u_int64_t tpid = TupleToIndex(t);
        covered_times[tpid] --;

        if (covered_times[tpid] == 0) {
            uncovered_tuples.emplace_back(t);
            sum_weight[tcid] -= use_weight ? weight[tpid] : 1;
            
            int upos = unique_idx[tpid];
            t_combnum comb = unique_pos[upos];
            for (int i = 0; i < strength; i ++) {
                int p = abs(t.v[i]) - 1;
                sum_weight_var[tcid][p] -= use_weight ? weight[tpid] : 1;

                int pos = comb.idx[i];
                int sz = unique_covered_tuples[tcid][p].size();
                if (pos < sz) {
                    swap(unique_covered_tuples[tcid][p][pos], unique_covered_tuples[tcid][p][sz - 1]);

                    u_int64_t tt_pid = unique_covered_tuples[tcid][p][pos];
                    t_tuple tt = IndexToTuple(tt_pid);
                    int tt_upos = unique_idx[tt_pid];
                    int k = 0;
                    for ( ; k < strength; k ++)
                        if (abs(tt.v[k]) - 1 == p)
                            break ;
                    unique_pos[tt_upos].idx[k] = pos;
                }
                unique_covered_tuples[tcid][p].pop_back();
            }
        
            
            int sz = all_unique_covered_tuples.size();
            if (upos >= sz)
                std::cout << "ERROR !" << std::endl, exit(0);
            if (upos != sz) {
                swap(all_unique_covered_tuples[upos], all_unique_covered_tuples[sz - 1]);
                swap(covered_tcid[upos], covered_tcid[sz - 1]);
                swap(unique_pos[upos], unique_pos[sz - 1]);
                unique_idx[all_unique_covered_tuples[upos]] = upos;
            }
            all_unique_covered_tuples.pop_back();
            covered_tcid.pop_back();
            unique_pos.pop_back();
        }
        if (covered_times[tpid] == 1) {
            unique_id.emplace_back(tpid);
            unique_tp.emplace_back(t);
        }
    }

    int new_unique_num = unique_id.size();
    for (int i = 0; i < new_unique_num; i ++) {
        u_int64_t tpid = unique_id[i];
        t_tuple t = unique_tp[i];

        unique_idx[tpid] = all_unique_covered_tuples.size();
        all_unique_covered_tuples.emplace_back(tpid);

        for (int j = 0; j < testcase_size; j ++)
            if (is_covered(testcases[j], t) && j != tcid) {

                covered_tcid.emplace_back(testcase_pos_to_idx[j]);
                sum_weight[j] += use_weight ? weight[tpid] : 1;
                t_combnum comb;
                for (int k = 0; k < strength; k ++) {
                    int p = abs(t.v[k]) - 1;
                    // unique_covered_tuples[j][p].insert(tpid);
                    comb.idx[k] = unique_covered_tuples[j][p].size();
                    unique_covered_tuples[j][p].emplace_back(tpid);
                    sum_weight_var[j][p] += use_weight ? weight[tpid] : 1;
                }
                unique_pos.emplace_back(comb);
                break ;
            }
    }

    unique_covered_tuples[tcid][b].clear();

    vector<int>& cur_clauses_cov = clauses_cov[tcid];

    change_bit(tc[bit] ? (- bit - 1) : (bit + 1), 1, tc, cur_clauses_cov);
    tc[bit] ^= 1;
    if (b != bit) {
        change_bit(- b - 1, 1, tc, cur_clauses_cov);
        tc[b] ^= 1;
    }
    
    for (t_combnum comb = begin_combnum(strength - 1); 
        comb.idx[strength - 2] < columns_num; next_comb(comb, strength - 1)) {
        
        t_tuple t;
        int i = 0;
        for (; i < strength - 1 && columns[comb.idx[i]] < bit; i ++) {
            int pi = columns[comb.idx[i]], vi = tc[pi];
            t.v[i] = vi ? (pi + 1) : (- pi - 1);
        }
        t.v[i] = tc[bit] ? (bit + 1) : (- bit - 1);
        for (; i < strength - 1; i ++) {
            int pi = columns[comb.idx[i]], vi = tc[pi];
            t.v[i + 1] = vi ? (pi + 1) : (- pi - 1);
        }
        u_int64_t tpid = TupleToIndex(t);
        // std::cout << bit << ", " << tpid << ", ", t.print(strength);
        covered_times[tpid] ++;

        if (covered_times[tpid] == 1) {
            unique_idx[tpid] = all_unique_covered_tuples.size();
            all_unique_covered_tuples.emplace_back(tpid);
            covered_tcid.emplace_back(testcase_pos_to_idx[tcid]);
            sum_weight[tcid] += use_weight ? weight[tpid] : 1;
            t_combnum comb;
            for (int i = 0; i < strength; i ++) {
                int p = abs(t.v[i]) - 1;
                // unique_covered_tuples[tcid][p].insert(tpid);
                comb.idx[i] = unique_covered_tuples[tcid][p].size();
                unique_covered_tuples[tcid][p].emplace_back(tpid);
                sum_weight_var[tcid][p] += use_weight ? weight[tpid] : 1;
            }
            unique_pos.emplace_back(comb);
        }
        
        else if (covered_times[tpid] == 2) {
            int unique_p = unique_idx[tpid];

            int tc_idx = covered_tcid[unique_p];
            int tc_pos = testcase_idx_to_pos[tc_idx];
            sum_weight[tc_pos] -= use_weight ? weight[tpid] : 1;
            t_combnum comb = unique_pos[unique_p];
            for (int i = 0; i < strength; i ++) {
                int ci = abs(t.v[i]) - 1;
                int pi = comb.idx[i];
                int sz = unique_covered_tuples[tc_pos][ci].size();
                
                if (pi >= sz) {
                    std::cout << "tpid =  " << tpid << ", ", t.print(strength);
                    comb.print(strength);
                    std::cout << "bit, b = " << bit << ", " << b << std::endl;
                    std::cout << "ci = " << ci << std::endl;
                    std::cout << pi << " " << sz << std::endl;
                    std::cout << "MEOW ! ERROR !" << std::endl, exit(0);
                }
                if (pi != sz - 1) {
                    swap(unique_covered_tuples[tc_pos][ci][pi], unique_covered_tuples[tc_pos][ci][sz - 1]);

                    u_int64_t tt_pid = unique_covered_tuples[tc_pos][ci][pi];
                    t_tuple tt = IndexToTuple(tt_pid);
                    int tt_upos = unique_idx[tt_pid];
                    int k = 0;
                    for ( ; k < strength; k ++)
                        if (abs(tt.v[k]) - 1 == ci)
                            break ;
                    unique_pos[tt_upos].idx[k] = pi;
                }
                unique_covered_tuples[tc_pos][ci].pop_back();
                sum_weight_var[tc_pos][ci] -= use_weight ? weight[tpid] : 1;
            }

            int sz = all_unique_covered_tuples.size();
            if (unique_p != sz - 1) {
                swap(all_unique_covered_tuples[unique_p], all_unique_covered_tuples[sz - 1]);
                swap(covered_tcid[unique_p], covered_tcid[sz - 1]);
                swap(unique_pos[unique_p], unique_pos[sz - 1]);
                unique_idx[all_unique_covered_tuples[unique_p]] = unique_p;
            }
            all_unique_covered_tuples.pop_back();
            covered_tcid.pop_back();
            unique_pos.pop_back();
        }
    }

    vector<int> break_pos;
    int uncovered_cnt = uncovered_tuples.size();
    for (int i = 0; i < uncovered_cnt; i ++) {
        t_tuple t = uncovered_tuples[i];
        if (is_covered(tc, t))
            break_pos.emplace_back(i);
    }
    int break_num = break_pos.size();
    for(int i = break_num - 1; i >= 0; i --) {
        int p = break_pos[i];
        if(p != uncovered_cnt)
            uncovered_tuples[p] = uncovered_tuples[uncovered_cnt - 1];
        uncovered_tuples.pop_back();
        uncovered_cnt --;
    }

    // check_unique (3);
}

void GlobalOptimizer:: force_flip_more (int tcid, const vector<int>& bits) { // group : 0 -> 1
    for (int x : bits)
        force_flip_one(tcid, x);
}

void GlobalOptimizer:: backtrack_row (int tcid, t_tuple tp) {

    for (const vector<int>& tc : bestArray)
        if (is_covered (tc, tp)) {
            forcetestcase (tcid, tc);
            break ;
        }
}

void GlobalOptimizer:: remove_unnecessary_tc () {
    
    int mid = gen() % testcase_size;
    for (int i = mid; i < testcase_size; ) {
        if (! sum_weight[i]) {
            updateinfo_remove_testcase(i);
            std::cout << "\033[;32mc current " << "remove " << i << ", ";
            std::cout << "\033[;32mc current " << strength << "-wise CA size: " 
                << testcase_size << "\033[0m" << std::endl;
            return ;
        }
        else i ++;
    }

    for (int i = mid - 1; i >= 0; i --)
        if (! sum_weight[i]) {
            updateinfo_remove_testcase(i);
            std::cout << "\033[;32mc current " << "remove " << i << ", ";
            std::cout << "\033[;32mc current " << strength << "-wise CA size: " 
                << testcase_size << "\033[0m" << std::endl;
            return ;
        }
}

void GlobalOptimizer:: random_greedy_step (t_tuple tp, long long maxi) {

    vector<int> besttcids;

    for (int i = 0; i < testcase_size; i ++) {    
        if (is_taboo(i, tp))
            continue ;
        
        const vector<int>& tc = testcases[i];
        vector<int>& cur_clauses_cov = clauses_cov[i];
        if (! check_force_tuple(tc, tp, cur_clauses_cov))
            continue ;
        
        vector<int> different_vars = get_different_vars(tc, tp);
        if (different_vars.empty())
            continue ;
        
        auto res = get_gain_for_flip_more(i, different_vars);
        long long net_gain = res.second - res.first;
        if (net_gain == maxi)
            besttcids.emplace_back(i);
        else if (net_gain > maxi) {
            maxi = net_gain;
            besttcids.clear();
            besttcids.emplace_back(i);
        }
    }

    // std::cout << besttcids.size() << "\n";
    if (! besttcids.empty()) {
        int besttcid = besttcids[gen() % besttcids.size()];
        set_taboo (besttcid, tp);
        vector<int> different_vars = get_different_vars(testcases[besttcid], tp);
        force_flip_more (besttcid, different_vars);
        return ;
    }
    
    backtrack_row (gen() % testcase_size, tp);
}

void GlobalOptimizer:: random_greedy_step (long long maxi) {

    int uncovered_cnt = uncovered_tuples.size();
    int picked_tuple = gen() % uncovered_cnt;
    t_tuple tp = uncovered_tuples[picked_tuple];

    vector<int> besttcids;

    for (int i = 0; i < testcase_size; i ++) {    
        if (is_taboo(i, tp))
            continue ;
        
        const vector<int>& tc = testcases[i];
        vector<int>& cur_clauses_cov = clauses_cov[i];
        if (! check_force_tuple(tc, tp, cur_clauses_cov))
            continue ;
        
        vector<int> different_vars = get_different_vars(tc, tp);
        if (different_vars.empty())
            continue ;
        
        auto res = get_gain_for_flip_more(i, different_vars);
        long long net_gain = res.second - res.first;
        if (net_gain == maxi)
            besttcids.emplace_back(i);
        else if (net_gain > maxi) {
            maxi = net_gain;
            besttcids.clear();
            besttcids.emplace_back(i);
        }
    }

    if (gen() % 100 < 10 && greedy_step_forced (tp))
        return ;

    // std::cout << besttcids.size() << "\n";
    if (! besttcids.empty()) {
        int besttcid = besttcids[gen() % besttcids.size()];
        set_taboo (besttcid, tp);
        vector<int> different_vars = get_different_vars(testcases[besttcid], tp);
        force_flip_more (besttcid, different_vars);
        return ;
    }

    if (gen() % 100 < __forced_greedy_percent && greedy_step_forced (tp))
        return ;
    
    backtrack_row (gen() % testcase_size, tp);
}

void GlobalOptimizer:: gradient_descent () {
    int uncovered_cnt = uncovered_tuples.size();
    int offset = gen() % uncovered_cnt;
    if (offset != 0)
        swap(uncovered_tuples[0], uncovered_tuples[offset]);

    long long maxi = 0;
    vector<pair<int, int> > best_choices;
    vector<pair<int,int> > first_best;

    for (int i = 0; i < uncovered_cnt; i ++) {

        t_tuple tp = uncovered_tuples[i];
        u_int64_t tpid = TupleToIndex(tp);

        for (int j = 0; j < testcase_size; j ++) {
            if (is_taboo (j, tp))
                continue ;

            const vector<int> &tc = testcases[j];

            pair<int, int> difference = get_different_var(tc, tp);
            if (difference.first != 1)
                continue ;
            
            int different_var = difference.second;
            vector<int>& cur_clauses_cov = clauses_cov[j];

            bool legal = check_force_tuple(tc, tp, cur_clauses_cov);
            if (! legal)
                continue ;

            auto res = get_gain_for_flip_one(j, different_var);
            long long net_gain = res.second - res.first;

            if (best_choices.empty() || net_gain == maxi) {
                maxi = net_gain;
                best_choices.emplace_back(make_pair(different_var, j));
                if (i == 0)
                    first_best.emplace_back(make_pair(different_var, j));
            }
            else if (net_gain > maxi) {
                best_choices.clear();
                best_choices.emplace_back(make_pair(different_var, j));
                maxi = net_gain;
                if (i == 0) {
                    first_best.clear();
                    first_best.emplace_back(make_pair(different_var, j));
                }
            }
        }
    }
    
    if (maxi > 0 && (! best_choices.empty())) {

        int choice = gen() % best_choices.size();
        int var = best_choices[choice].first;
        int besttcid = best_choices[choice].second;

        set_taboo (besttcid, var);
        force_flip_one (besttcid, var);
        return ;
    }
    
    if (use_weight)
        for (t_tuple t : uncovered_tuples) {
            u_int64_t tpid = TupleToIndex(t);
            weight[tpid] ++;
        }

    if (gen () % 1000 < 1) {
        t_tuple tp = uncovered_tuples[0];
        backtrack_row (gen() % testcase_size, tp);
        return ;
    }

    if (gen () % 1000 < 1) {
        t_tuple tp = uncovered_tuples[0];
        greedy_step_forced(tp);
        return ;
    }

    if (! first_best.empty()) {
        int choice = gen() % first_best.size();
        int var = first_best[choice].first;
        int besttcid = first_best[choice].second;
        set_taboo (besttcid, var);
        force_flip_one (besttcid, var);
        return ;
    }
    t_tuple tp = uncovered_tuples[0];
    return random_greedy_step (maxi);
}

pair<u_int64_t, u_int64_t> GlobalOptimizer:: get_gain_for_forcetestcase (int tcid, const vector<int>& tc2) {

    const vector<int>& tc = testcases[tcid];

    u_int64_t break_cnt = sum_weight[tcid];
    u_int64_t gain_cnt = 0;

    map<u_int64_t, bool> visited;
    for (int var = 0; var < nvar; var ++)
        for (u_int64_t tpid : unique_covered_tuples[tcid][var]) {
            if (visited.count(tpid))
                continue ;
            visited[tpid] = true;
            t_tuple t = IndexToTuple(tpid);
            if (is_covered(tc2, t))
                gain_cnt += use_weight ? weight[tpid] : 1;
        }
    visited.clear();

    for (t_tuple t : uncovered_tuples)
        if (is_covered(tc2, t))
            gain_cnt += use_weight ? weight[TupleToIndex(t)] : 1;
    return {break_cnt, gain_cnt};
}

bool GlobalOptimizer:: greedy_step_forced (t_tuple tp) {
    
    cdcl_solver->clear_assumptions();
    for (int i = 0; i < strength; i ++) {
        int pi = abs(tp.v[i]) - 1, vi = tp.v[i] > 0;
        cdcl_solver->add_assumption(pi, vi);
    }

    long long maxi = 0;
    vector<int> besttcids;
    vector<vector<int> > besttc2s;

    for (int i = 0; i < testcase_size; i ++) {

        if (is_taboo(i, tp))
            continue ;

        for(int j = 0; j < nvar; j ++)
            cdcl_solver->setPolarity(j, testcases[i][j] == 0);
        vector<int> tc2 = vector<int>(nvar, 0);
        bool ret = cdcl_solver->solve();
        if(! ret) {
            cout << "c \033[1;31mError: SAT solve failing!\033[0m" << endl;
            exit(0);
        }
        cdcl_solver->get_solution(tc2);

        auto res = get_gain_for_forcetestcase(i, tc2);
        long long net_gain = res.second - res.first;

        if (besttcids.empty() || net_gain == maxi) {
            maxi = net_gain;
            besttcids.emplace_back(i);
            besttc2s.emplace_back(tc2);
        }
        else if(net_gain > maxi) {
            maxi = net_gain;
            besttcids.clear(), besttcids.emplace_back(i);
            besttc2s.clear(), besttc2s.emplace_back(tc2);
        }
    }

    if (besttcids.empty())
        return false;
    
    int idx = gen() % besttcids.size();
    int besttcid = besttcids[idx];
    vector<int> besttc2 = besttc2s[idx];

    set_taboo(besttcid, besttc2);
    forcetestcase(besttcid, besttc2);

    return true;
}

void GlobalOptimizer:: search () {

    int cur_step = 0, last_success_step = 0;

    for ( ; last_success_step < stop_length; ) {
        
        // remove_unnecessary_tc();

        if(uncovered_tuples.empty()) {

            remove_unnecessary_tc();
            
            bestArray = testcases;
            cout << "\033[;32mc current " << strength << "-wise CA size: " 
                << testcase_size << ", step #" << cur_step << " \033[0m" << endl;
            remove_testcase_greedily ();
            // remove_testcase_randomly ();
            last_success_step = 0;
            continue ;
        }

        cur_step ++;
        last_success_step ++;
        gradient_descent ();
    }

    if(! uncovered_tuples.empty()) {
        testcases = bestArray;
        testcase_size = testcases.size();
    }
}

void GlobalOptimizer:: SaveTestcaseSet() {
    SaveTestcaseSet(output_testcase_path);
}

void GlobalOptimizer:: SaveTestcaseSet(string result_path) {
    
    ofstream res_file(result_path);
    cout << testcases.size() << "\n";
    for (const vector<int>& testcase: testcases) {
        for (int v = 0; v < nvar; v++)
            res_file << testcase[v] << " ";
        res_file << "\n";
    }
    res_file.close();
    cout << "c Testcase set saved in " << result_path << endl;
}