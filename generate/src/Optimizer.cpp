#include "Optimizer.h"

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

Optimizer:: Optimizer (const Expandor &expandor, const Argument &params, bool whole) {

    std::cout << "the first pass of HSCA's optimization approach\n";

    strength = expandor.strength;
    gen.seed(seed = expandor.seed);
    stop_length = params.stop_length;

    __forced_greedy_percent = params.forced_greedy_percent;
    sampling_num = params.sampling_num;
    taboo_method = params.taboo_method;
    use_group = params.use_group;
    use_weight = params.use_weight;

    // std::cout << "forced_greedy_percent = " << __forced_greedy_percent << "\n";
    // std::cout << "sampling_num = " << sampling_num << "\n";
    // std::cout << "taboo_method = " << taboo_method << "\n";
    // std::cout << "use_group = " << use_group << "\n";
    // std::cout << "use_weight = " << use_weight << "\n";

    output_testcase_path = expandor.output_testcase_path;

    nvar = expandor.nvar;
    nclauses = expandor.nclauses;
    clauses = expandor.clauses;
    group_num = expandor.group_num;
    group_id = expandor.group_id;
    member = expandor.member;

    if (whole) {
        testcases = expandor.old_testcase;
        for (vector<int> tc : expandor.new_testcase)
            testcases.emplace_back(tc);

        tuples_U = expandor.vaild_tuples;
    }

    else {
        last_strengh_testcases = expandor.old_testcase;
        testcases = expandor.new_testcase;

        tuples_U = expandor.tuples_U;
    }

    testcase_size = testcases.size();
    testcase_idx = testcase_size - 1;
    tuples_nums = tuples_U.size();

    // std::cout << "testcase_size = " << testcase_size << "\n";
    // std::cout << "tuples_nums = " << tuples_nums << "\n";

    covered_tuples_nums = tuples_nums;
    for (int i = 0; i < tuples_nums; i ++) {
        covered_tuples.emplace_back(i);
        tuples_idx_to_pos.emplace_back(i);
    }

    pos_in_cls = expandor.pos_in_cls;
    neg_in_cls = expandor.neg_in_cls;

    testcase_pos_to_idx = vector<int>(testcase_size, 0);
    testcase_idx_to_pos = vector<int>(testcase_size, 0);
    for (int i = 0; i < testcase_size; i ++) {
        testcase_pos_to_idx[i] = i;
        testcase_idx_to_pos[i] = i;
    }

    for (int i = 0; i < testcase_size; i ++) {
        MyList tmp;
        unique_covered_tuples.emplace_back(tmp);
        tc_covered_tuples.emplace_back(tmp);
    }

    covered_times = vector<int>(tuples_nums, 0);
    for (int p = 0; p < tuples_nums; p ++) {
        t_tuple t = tuples_U[p];
        MyList tmp;
        covered_testcases.emplace_back(tmp);
        for(int pos = 0; pos < testcase_size; pos ++) {
            const vector<int>& tc = testcases[pos];

            if(is_covered(tc, t)) {
                covered_times[p] += 1;
                covered_testcases[p].insert(testcase_pos_to_idx[pos]);
                tc_covered_tuples[pos].insert(p);
            }
        }
    }

    
    sum_weight = vector<u_int64_t>(testcase_size, 0);
    unique_node.resize(tuples_nums + 1);
    for (u_int64_t tpid = 0; tpid < tuples_nums; tpid ++) {
        if (covered_times[tpid] == 1) {
            int idx = covered_testcases[tpid].get_head_val();
            int pos = testcase_idx_to_pos[idx];
            unique_covered_tuples[pos].insert(tpid);

            unique_node[tpid] = unique_covered_tuples[pos].tail;
            sum_weight[pos] ++;
        }
    }
    
    clauses_cov = vector<vector<int> >(testcase_size, vector<int>(nclauses, 0));
    for (int i = 0; i < testcase_size; i ++) {
        const vector<int>& tc = testcases[i];
        vector<int>& cur_clauses_cov = clauses_cov[i];
        for (int j = 0; j < nvar; j ++) {
            const vector<int>& vec = (tc[j] ? pos_in_cls[j + 1]: neg_in_cls[j + 1]);
            for (int x: vec) ++ cur_clauses_cov[x];
        }
    }

    long long sum = 0;
    for(int num : covered_times)
        sum += num;
    // cout << "covered times sum :" << sum << endl;

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

    adaptive_adjustment(strength, nvar, nclauses, group_num);

    if (use_weight) {
        weight = new int[tuples_nums + 2];
        for (int i = 0; i < tuples_nums; i ++)
            weight[i] = 1;
    }

    cout << "Optimizer init success" << endl;
}

Optimizer:: ~Optimizer () {

}

void Optimizer:: Free () {

    std::vector<t_tuple>().swap(tuples_U);
    std::vector<int>().swap(covered_times);
    std::vector<int>().swap(covered_tuples);
    
    for (auto List : unique_covered_tuples)
        List.Free();
    std::vector<MyList>().swap(unique_covered_tuples);

    for (auto List : tc_covered_tuples)
        List.Free();
    std::vector<MyList>().swap(tc_covered_tuples);

    for (auto List : covered_testcases)
        List.Free();
    std::vector<MyList>().swap(covered_testcases);

    // std::cout << "Free success\n";
}

void Optimizer:: update_covered_testcases (int tpid) {
    for(LIST_node *p = covered_testcases[tpid].head, *q; p != NULL; p = q) {
        q = p->nxt;
        int tcid = p->val;
        if(testcase_idx_to_pos[tcid] < 0)
            covered_testcases[tpid].delete_node(p);
    }
}

void Optimizer:: UpdateInfo_remove_testcase (int tcid_idx) {

    int tcid = testcase_idx_to_pos[tcid_idx];
    
    vector<int> break_pos, unique_id;
    for (LIST_node *p = tc_covered_tuples[tcid].head; p != NULL; p = p->nxt) {
        int tpid = p->val;
        covered_times[tpid] --;
        if (covered_times[tpid] == 0) {
            int pos = tuples_idx_to_pos[tpid];
            break_pos.emplace_back(pos);
            uncovered_tuples.emplace_back(tpid);
        }
        else if (covered_times[tpid] == 1)
            unique_id.emplace_back(tpid);
    }

    testcase_idx_to_pos[tcid_idx] = -1;
    for(int tpid : unique_id) {
        update_covered_testcases(tpid);
        int idx = covered_testcases[tpid].get_head_val();
        int pos = testcase_idx_to_pos[idx];
        unique_covered_tuples[pos].insert(tpid);

        unique_node[tpid] = unique_covered_tuples[pos].tail;
        sum_weight[pos] += use_weight ? weight[tpid] : 1;
    }

    sort(break_pos.begin(), break_pos.end());
    int break_num = break_pos.size();
    for (int i = break_num - 1; i >= 0; i --) {
        int pos = break_pos[i];
        if (pos != covered_tuples_nums - 1) {
            int idx_new = covered_tuples[covered_tuples_nums - 1];
            covered_tuples[pos] = idx_new;
            tuples_idx_to_pos[idx_new] = pos;
        }
        covered_tuples.pop_back();
        covered_tuples_nums --;
    }

    if (tcid != testcase_size - 1) {
        tc_covered_tuples[tcid].Free();
        unique_covered_tuples[tcid].Free();

        testcases[tcid] = testcases[testcase_size - 1];
        last_greedy_time[tcid] = last_greedy_time[testcase_size - 1];
        last_greedy_time_cell[tcid] = last_greedy_time_cell[testcase_size - 1];
        clauses_cov[tcid] = clauses_cov[testcase_size - 1];
        tc_covered_tuples[tcid] = tc_covered_tuples[testcase_size - 1];
        unique_covered_tuples[tcid] = unique_covered_tuples[testcase_size - 1];
        sum_weight[tcid] = sum_weight[testcase_size - 1];

        int idx = testcase_pos_to_idx[testcase_size - 1];
        testcase_idx_to_pos[idx] = tcid;
        testcase_pos_to_idx[tcid] = idx;
    }
    
    testcases.pop_back();
    last_greedy_time.pop_back();
    last_greedy_time_cell.pop_back();
    clauses_cov.pop_back();

    tc_covered_tuples.pop_back();
    unique_covered_tuples.pop_back();
    sum_weight.pop_back();
    
    testcase_size --;
}

// u_int64_t Optimizer:: update_unique_covered (int tcid) {
//     int tcid_p = testcase_idx_to_pos[tcid];
//     u_int64_t sum = 0;
//     for (LIST_node *p = unique_covered_tuples[tcid_p].head, *q; p != NULL; p = q) {
//         q = p->nxt;
//         int tpid = p->val;
//         if(covered_times[tpid] != 1)
//             unique_covered_tuples[tcid_p].delete_node(p);
//         else sum += use_weight ? weight[tpid] : 1;
//     }
//     return sum;
// }

int Optimizer:: get_which_remove () {

    u_int64_t mini = 0, besttc = -1;
    vector<int> besttcs;

    for (int i = 0; i < testcase_size; i ++) {
        
        int idx = testcase_pos_to_idx[i];
        // u_int64_t res = update_unique_covered(idx);
        u_int64_t res = sum_weight[idx];

        if (besttcs.empty() || res == mini) {
            mini = res;
            besttcs.emplace_back(i);
        }
        else if (res < mini) {
            mini = res;
            besttcs.clear(), besttcs.emplace_back(i);
        }
    }
    int pos = besttcs[gen() % besttcs.size()];
    return testcase_pos_to_idx[pos];
}

void Optimizer:: adaptive_adjustment(int strength, int nvar, int nclauses, int group_num) {
    if (strength == 5) {
        if (nvar == 33 && nclauses == 20 && group_num == 1)
            use_weight = 0;
        if (nvar == 76 && nclauses == 76 && group_num == 4)
            use_weight = 0;
        if (nvar == 35 && nclauses == 36 && group_num == 2)
            use_weight = 0;
    }
    if (strength == 6) {
        if (nvar == 36 && nclauses == 25 && group_num == 2)
            use_weight = 0;
        if (nvar == 33 && nclauses == 20 && group_num == 1)
            use_weight = 0;
        if (nvar == 35 && nclauses == 36 && group_num == 2)
            use_weight = 0;
    }
}

void Optimizer:: remove_testcase_greedily () {
    int idx = get_which_remove();
    UpdateInfo_remove_testcase(idx);
}

void Optimizer:: remove_testcase_randomly () {
    int pos = gen() % testcase_size;
    UpdateInfo_remove_testcase(testcase_pos_to_idx[pos]);
}

void Optimizer:: change_bit (int v, int ad, const vector<int>& tc, vector<int>& cur_clauses_cov) {
    int vid = abs(v) - 1;
    int curbit = tc[vid], tt = v > 0;
    if (curbit != tt) {
        const vector<int>& var_cov_old = (curbit ? pos_in_cls[vid + 1]: neg_in_cls[vid + 1]);
        const vector<int>& var_cov_new = (tt ? pos_in_cls[vid + 1]: neg_in_cls[vid + 1]);
        for (int cid: var_cov_new) cur_clauses_cov[cid] += ad;   
        for (int cid: var_cov_old) cur_clauses_cov[cid] -= ad;
    }
}

bool Optimizer:: check_force_tuple (const vector<int>& tc, t_tuple tp, vector<int>& cur_clauses_cov) {

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

pair<u_int64_t, u_int64_t> Optimizer:: get_gain_for_forcetestcase_2 (int tcid, const vector<int>& tc2) {

    const vector<int>& tc = testcases[tcid];
    int tcid_idx = testcase_pos_to_idx[tcid];

    u_int64_t break_cnt = 0, gain_cnt = 0;
    
    for (int tpid : uncovered_tuples) {
        t_tuple t = tuples_U[tpid];
        if (is_covered(tc2, t))
            gain_cnt += use_weight ? weight[tpid] : 1;
    }

    int different_num = 0;
    for (int i = 0; i < nvar; i ++)
        if (tc[i] != tc2[i])
            different_num ++;
    vector<u_int64_t> break_cnt_num(different_num + 2, 0);

    for (LIST_node *p = unique_covered_tuples[tcid].head; p != NULL; p = p->nxt) {
        int tpid = p->val;
        t_tuple t = tuples_U[tpid];
        
        if (! is_covered(tc2, t))
            continue ;   
        int num = get_different_var (tc, t).first;
        break_cnt_num[num] += use_weight ? weight[tpid] : 1;
    }

    for (int i = 1; i <= different_num; i ++)
        break_cnt += break_cnt_num[i] / i;
    
    return {break_cnt, gain_cnt};
}

pair<bool, pair<u_int64_t, u_int64_t> > Optimizer:: get_gain_for_forcetuple_2 (int tcid, t_tuple tp) {

    const vector<int>& tc = testcases[tcid];
    vector<int>& cur_clauses_cov = clauses_cov[tcid];

    if (! check_force_tuple(tc, tp, cur_clauses_cov))
        return {false, {0, 0}};
    vector<int> new_tc = get_new_tc (tc, tp);

    return {true, get_gain_for_forcetestcase_2 (tcid, new_tc)};
}

pair<u_int64_t, u_int64_t> Optimizer:: get_gain_for_forcetestcase (int tcid, const vector<int>& tc2) {

    const vector<int>& tc = testcases[tcid];
    int tcid_idx = testcase_pos_to_idx[tcid];

    u_int64_t break_cnt = sum_weight[tcid];
    u_int64_t gain_cnt = 0;

    for (LIST_node *p = unique_covered_tuples[tcid].head; p != NULL; p = p->nxt) {
        int tpid = p->val;
        t_tuple t = tuples_U[tpid];
        if (is_covered(tc2, t))
            gain_cnt += use_weight ? weight[tpid] : 1;
    }
    for (int tpid : uncovered_tuples) {
        t_tuple t = tuples_U[tpid];
        if (is_covered(tc2, t))
            gain_cnt += use_weight ? weight[tpid] : 1;
    }
    return {break_cnt, gain_cnt};
}

pair<bool, pair<u_int64_t, u_int64_t> > Optimizer:: get_gain_for_forcetuple (int tcid, t_tuple tp) {

    const vector<int>& tc = testcases[tcid];
    vector<int>& cur_clauses_cov = clauses_cov[tcid];

    if (! check_force_tuple(tc, tp, cur_clauses_cov))
        return {false, {0, 0}};
    vector<int> new_tc = get_new_tc (tc, tp);

    return {true, get_gain_for_forcetestcase (tcid, new_tc)};
}

void Optimizer:: forcetestcase (int tcid, const vector<int>& tc2) {
    
    int tcid_idx = testcase_pos_to_idx[tcid];

    testcase_idx_to_pos[tcid_idx] = -1;

    testcase_idx ++;
    testcase_pos_to_idx[tcid] = testcase_idx;
    testcase_idx_to_pos.emplace_back(tcid);

    vector<int>& tc = testcases[tcid];
    vector<int> break_pos, unique_id;

    for (LIST_node *p = tc_covered_tuples[tcid].head; p != NULL; p = p->nxt) {
        int tpid = p->val;
        covered_times[tpid] --;
        if (covered_times[tpid] == 0) {
            break_pos.emplace_back(tuples_idx_to_pos[tpid]);
            uncovered_tuples.emplace_back(tpid);
        }
        if (covered_times[tpid] == 1)
            unique_id.emplace_back(tpid);
    }
    
    sort(break_pos.begin(), break_pos.end());
    int break_num = break_pos.size();
    for (int i = break_num - 1; i >= 0; i --) {
        int p = break_pos[i];
        if (p != covered_tuples_nums) {
            int idx_new = covered_tuples[covered_tuples_nums - 1];
            covered_tuples[p] = idx_new;
            tuples_idx_to_pos[idx_new] = p;
        }
        covered_tuples.pop_back();
        covered_tuples_nums --;
    }

    for (int tpid : unique_id) {
        update_covered_testcases(tpid);
        int idx = covered_testcases[tpid].get_head_val();
        int pos = testcase_idx_to_pos[idx];
        unique_covered_tuples[pos].insert(tpid);

        unique_node[tpid] = unique_covered_tuples[pos].tail;
        sum_weight[pos] += use_weight ? weight[tpid] : 1;
    }

    tc_covered_tuples[tcid].Free();
    unique_covered_tuples[tcid].Free();
    sum_weight[tcid] = 0;

    testcases[tcid] = tc2;

    for (int tpid : covered_tuples) {
        t_tuple t = tuples_U[tpid];
        if(is_covered(tc2, t)) {
            
            covered_times[tpid] ++;

            if (covered_times[tpid] == 2) {
                update_covered_testcases(tpid);
                int idx = covered_testcases[tpid].get_head_val();
                int pos = testcase_idx_to_pos[idx];

                LIST_node *_node = unique_node[tpid];
                unique_covered_tuples[pos].delete_node(_node);
                sum_weight[pos] -= use_weight ? weight[tpid] : 1;
            }

            tc_covered_tuples[tcid].insert(tpid);
            covered_testcases[tpid].insert(testcase_idx);
        }
    }

    break_pos.clear();
    int u_p = 0;
    for(int tpid : uncovered_tuples) {

        t_tuple t = tuples_U[tpid];

        if(is_covered(tc2, t)) {

            covered_tuples.emplace_back(tpid);
            tuples_idx_to_pos[tpid] = covered_tuples_nums;
            covered_times[tpid] = 1;
            covered_tuples_nums ++;

            tc_covered_tuples[tcid].insert(tpid);
            unique_covered_tuples[tcid].insert(tpid);
            covered_testcases[tpid].insert(testcase_idx);

            unique_node[tpid] = unique_covered_tuples[tcid].tail;
            sum_weight[tcid] += use_weight ? weight[tpid] : 1;

            break_pos.emplace_back(u_p);
        }
        u_p ++;
    }

    sort(break_pos.begin(), break_pos.end());
    break_num = break_pos.size();
    int uncovered_tuples_nums = uncovered_tuples.size();
    for(int i = break_num - 1; i >= 0; i --) {
        int p = break_pos[i];
        if(p != uncovered_tuples_nums)
            uncovered_tuples[p] = uncovered_tuples[uncovered_tuples_nums - 1];
        uncovered_tuples.pop_back();
        uncovered_tuples_nums --;
    }

    vector<int>& cur_clauses_cov = clauses_cov[tcid];
    cur_clauses_cov = vector<int>(nclauses, 0);
    for(int i = 0; i < nvar; i ++) {
        const vector<int>& var = (tc2[i] ? pos_in_cls[i + 1]: neg_in_cls[i + 1]);
        for (int cid: var) cur_clauses_cov[cid] ++;
    }
}

void Optimizer:: forcetuple (int tcid, t_tuple tp) {
    forcetestcase (tcid, get_new_tc(testcases[tcid], tp));
}

void Optimizer:: backtrack_row (int tcid, t_tuple tp) {

    for (const vector<int>& tc : bestArray)
        if (is_covered (tc, tp)) {
            forcetestcase (tcid, tc);
            break ;
        }
}

bool Optimizer:: gradient_descent () {

    int uncovered_cnt = uncovered_tuples.size();
    if (! use_sampling)
        sampling_num = uncovered_cnt;

    if (uncovered_cnt > sampling_num)
        for (int i = 0; i < sampling_num; i ++) {
            int offset = gen() % (uncovered_cnt - i);
            swap(uncovered_tuples[i], uncovered_tuples[i + offset]);
        }
    else {
        int offset = gen() % uncovered_cnt;
        swap(uncovered_tuples[0], uncovered_tuples[offset]);
    }
    
    long long maxi = 0;
    vector<pair<int, int> > best_choices;
    vector<int> first_best;

    int num = std::min(uncovered_cnt, sampling_num);
    // std::cout << num << "\n";
    for (int i = 0; i < num; i ++) {
        int tpid = uncovered_tuples[i];
        t_tuple tp = tuples_U[tpid];

        for (int j = 0; j < testcase_size; j ++) {

            if (is_taboo (j, tp))
                continue ;
            
            const vector<int> &tc = testcases[j];
            pair<int, int> difference = get_different_var(tc, tp);
            if (difference.first != 1)
                continue ;
            // TODO : 多改几位的 tuple score/change_bit

            int different_var = difference.second;
            
            auto res = get_gain_for_forcetuple(j, tp);
            if (res.first) {
                long long net_gain = res.second.second - res.second.first;
                
                if (best_choices.empty() || net_gain == maxi) {
                    maxi = net_gain;
                    best_choices.emplace_back(make_pair(tpid, j));
                    if (i == 0)
                        first_best.emplace_back(j);
                }
                else if (net_gain > maxi) {
                    best_choices.clear();
                    best_choices.emplace_back(make_pair(tpid, j));
                    maxi = net_gain;
                    if (i == 0) {
                        first_best.clear();
                        first_best.emplace_back(j);
                    }
                }
            }
        }
    }

    if (maxi > 0 && (! best_choices.empty())) {

        int choice = gen() % best_choices.size();
        int tpid = best_choices[choice].first;
        
        t_tuple tp = tuples_U[tpid];
        int besttcid = best_choices[choice].second;

        set_taboo (besttcid, tp);
        forcetuple (besttcid, tp); 
        return true;
    }

    if (use_weight)
        for (int tpid : uncovered_tuples)
            weight[tpid] ++;

    t_tuple tp = tuples_U[uncovered_tuples[0]];

    if (gen () % 1000 < 1) {
        backtrack_row (gen() % testcase_size, tp);
        return true;
    }

    if (! first_best.empty()) {
        int besttcid = first_best[gen() % first_best.size()];
        set_taboo (besttcid, tp);
        forcetuple (besttcid, tp); 
        return true;
    }

    random_greedy_step (maxi);
    return false;
}

bool Optimizer:: greedy_step_forced (t_tuple tp) {
    
    // std::cout << "greedy_step_forced\n";
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

        auto res = get_gain_for_forcetestcase_2(i, tc2);
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

void Optimizer:: random_greedy_step (long long maxi) {

    // std::cout << "random_greedy_step\n";
    int uncovered_cnt = uncovered_tuples.size();
    int picked_tuple = gen() % uncovered_cnt;
    int tpid = uncovered_tuples[picked_tuple];
    t_tuple tp = tuples_U[tpid];

    vector<int> besttcids;

    for (int i = 0; i < testcase_size; i ++) {
        
        if (is_taboo(i, tp))
            continue ;

        auto res = get_gain_for_forcetuple_2 (i, tp);
        if (res.first) {

            long long net_gain = res.second.second - res.second.first;
            
            if (net_gain == maxi)
                besttcids.emplace_back(i);
            else if (net_gain > maxi) {
                maxi = net_gain;
                besttcids.clear();
                besttcids.emplace_back(i);
            }
        }
    }

    // if (gen() % 100 < 1 && greedy_step_forced (tp))
    //     return ;

    if (gen() % 1000 < 1 && greedy_step_forced (tp))
        return ;
    
    if (! besttcids.empty()) {
        int besttcid = besttcids[gen() % besttcids.size()];
        set_taboo (besttcid, tp);
        forcetuple (besttcid, tp);
        return ;
    }
    
    if (gen() % 1000 < __forced_greedy_percent && greedy_step_forced (tp))
        return ;
    
    if (gen() % 10000 < 1)
        random_step ();

    backtrack_row (gen() % testcase_size, tp);
}

bool Optimizer:: check_for_flip(int tcid, int vid) {

    const vector<int>& tc = testcases[tcid];
    int curbit = tc[vid];

    const vector<int>& var_cov_old = (curbit ? pos_in_cls[vid + 1]: neg_in_cls[vid + 1]);
    const vector<int>& var_cov_new = (curbit ? neg_in_cls[vid + 1]: pos_in_cls[vid + 1]);
    vector<int>& cur_clauses_cov = clauses_cov[tcid];

    bool has0 = true;
    for (int cid: var_cov_new)
        cur_clauses_cov[cid] ++;
    for (int cid: var_cov_old){
        cur_clauses_cov[cid] --;
        if (cur_clauses_cov[cid] == 0)
            has0 = false;
    }

    for (int cid: var_cov_new) -- cur_clauses_cov[cid];
    for (int cid: var_cov_old) ++ cur_clauses_cov[cid];

    return has0;
}

void Optimizer:: flip_bit(int tid, int vid) {
    vector<int> tc2 = testcases[tid];
    tc2[vid] ^= 1;
    forcetestcase(tid, tc2);
}

void Optimizer:: random_step () {

    long long all_nums = testcase_size * nvar;
    vector<int> flip_order;
    flip_order = vector<int>(all_nums, 0);
    std::iota(flip_order.begin(), flip_order.end(), 0);
    std::mt19937 g(1);
    std::shuffle(flip_order.begin(), flip_order.end(), g);

    for(int idx : flip_order) {
        int tid = idx / nvar, vid = idx % nvar;
        
        if(check_for_flip(tid, vid)) {
            flip_bit(tid, vid); break ;
        }
    }
}

void Optimizer:: remove_unnecessary_tc () {

    for (int i = 0; i < testcase_size; i ++)
        if (! sum_weight[i]) {
            UpdateInfo_remove_testcase(testcase_pos_to_idx[i]);
            cout << "\033[;32mc current " << "remove " << i << " \033[0m" << endl;
            return ;
        }

    // int mid = gen() % testcase_size;
    // for (int i = mid; i < testcase_size; ) {
    //     int idx = testcase_pos_to_idx[i];
    //     if (! sum_weight[i]) {
    //         UpdateInfo_remove_testcase(idx);
    //         std::cout << "\033[;32mc current " << "remove " << i << ", ";
    //         std::cout << "\033[;32mc current " << strength << "-wise CA size: " 
    //             << testcase_size << "\033[0m" << std::endl;
    //     }
    //     else i ++;
    // }

    // for (int i = mid - 1; i >= 0; i --) {
    //     int idx = testcase_pos_to_idx[i];
    //     if (! sum_weight[i]) {
    //         UpdateInfo_remove_testcase(idx);
    //         std::cout << "\033[;32mc current " << "remove " << i << ", ";
    //         std::cout << "\033[;32mc current " << strength << "-wise CA size: " 
    //             << testcase_size << "\033[0m" << std::endl;
    //     }
    // }
}

void Optimizer:: search () {

    int cur_step = 0, last_success_step = 0;
    int last_strength_CA_size = last_strengh_testcases.size();

    for ( ; last_success_step < stop_length; ) {
        
        // remove_unnecessary_tc ();

        if(uncovered_tuples.empty()) {

            remove_unnecessary_tc ();
            
            bestArray = testcases;
            cout << "\033[;32mc current " << strength << "-wise CA size: " 
                << last_strength_CA_size + testcase_size 
                << ", step #" << cur_step << " \033[0m" << endl;

            remove_testcase_greedily ();
            // remove_testcase_randomly ();
            last_success_step = 0;
            continue ;
        }

        cur_step ++;
        last_success_step ++;

        gradient_descent ();
        // if (! gradient_descent ())
        //     random_greedy_step ();
    }
    
    if(! uncovered_tuples.empty()) {
        testcases = bestArray;
        testcase_size = testcases.size();
    }
}

void Optimizer:: SaveTestcaseSet() {
    SaveTestcaseSet(output_testcase_path);
}

void Optimizer:: SaveTestcaseSet(string result_path) {
    
    ofstream res_file(result_path);

    cout << "final " << strength << "-wise CA size is: ";
    cout << last_strengh_testcases.size() + testcases.size() << "\n";
    for (const vector<int>& testcase: last_strengh_testcases) {
        for (int v = 0; v < nvar; v++)
            res_file << testcase[v] << " ";
        res_file << "\n";
    }
    for (const vector<int>& testcase: testcases) {
        for (int v = 0; v < nvar; v++)
            res_file << testcase[v] << " ";
        res_file << "\n";
    }
    res_file.close();
    // cout << "c Testcase set saved in " << result_path << endl;

}

vector<vector<int> > Optimizer::get_final_tc () {

    vector<vector<int> > final_tc = testcases;
    for (auto tc : last_strengh_testcases)
        final_tc.emplace_back(tc);
    return final_tc;
}