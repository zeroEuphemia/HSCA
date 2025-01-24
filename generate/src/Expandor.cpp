#include "Expandor.h"

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
#include <pthread.h>

using std::vector;
using std::string;
using std::pair;
using std::mt19937_64;

Expandor:: Expandor (const Argument &params) {
    
    need_all_valid_tuples = false;
    
    // 设置参数
    strength = params.init_strength;
    input_cnf_path = params.input_cnf_path;
    init_CA_file_path = params.init_CA_file_path;
    output_testcase_path = params.output_testcase_path;
    thread_num = params.thread_num;

    gen.seed(seed = params.seed);

    candidate_set_size = params.candidate_set_size;
    use_cnf_reduction = params.use_cnf_reduction;

    use_cache = true;
    if (params.use_cache == 0)
        use_cache = false;
    use_invalid_expand = true;
    if (params.use_invalid_expand == 0)
        use_invalid_expand = false;
    // 设置 cnf 和求解器
    string real_input_cnf_path = input_cnf_path;
    if (use_cnf_reduction) {
        mt19937_64 tmp_gen(params.seed);
        string reduced_cnf_path = "/tmp/" + std::to_string(getpid()) + std::to_string(tmp_gen()) + "_reduced.cnf";
        string cmd = "./bin/coprocessor -enabled_cp3 -up -subsimp -no-bve -no-bce -no-dense -dimacs=" + 
            reduced_cnf_path + " " + params.input_cnf_path;
        
        int _ = system(cmd.c_str());
        real_input_cnf_path = reduced_cnf_path;
    }
    // std::cout << real_input_cnf_path << std::endl;
    CNFInfo cnf(real_input_cnf_path);
    cnf.dump(nvar, nclauses, clauses, pos_in_cls, neg_in_cls);
    
    // group 设置
    group_id.resize(nvar);
    for (int i = 0; i < nvar; i ++)
        group_id[i] = -1;

    // read group_file
    ifstream fin(params.group_file_path.c_str());
    fin >> group_num;
    member.resize(group_num);
    group_flag.resize(group_num);
    group_flag_array.resize(thread_num);
    for (int i = 0; i < thread_num; i ++)
        group_flag_array[i].resize(group_num);
    
    for (int i = 0; i < group_num; i ++) {
        int num; fin >> num;
        for (int j = 0; j < num; j ++) {
            int x; fin >> x; x --;
            group_id[x] = i;
            member[i].emplace_back(x);
        }
    }
    fin.close();

    // for (int i = 0; i < nvar; i ++)
    //     std::cout << "(" << i << ", " << group_id[i] << ")\n";

    // cdcl_solver = new CDCLSolver::Solver;
    // cdcl_solver->read_clauses(nvar, clauses);

    // 设置 2-wise CA
    FILE *in_ca = fopen(init_CA_file_path.c_str(), "r");
    MyBitSet tmp(nvar);
    vector<int> tep; tep.resize(nvar, 0);
    for (int c, p = 0; (c = fgetc(in_ca)) != EOF; ) {
        if (c == '\n') {
            CA.emplace_back(tmp);
            new_testcase.emplace_back(tep);
            p = 0;
        }
        else if (isdigit(c)){
            if(c == '1')
                tmp.set(p), tep[p] = 1;
            else tmp.unset(p), tep[p] = 0;
            p ++;
        }
    }
    
    count_each_var_uncovered[0].resize(nvar + 2, 0);
    count_each_var_uncovered[1].resize(nvar + 2, 0);

    combnum = new u_int64_t *[mx_strength + 1];
    for (int i = 2; i <= mx_strength; i ++)
        combnum[i] = new u_int64_t[nvar + 1]{0};
    for (int i = 2; i <= nvar; i ++) {
        combnum[2][i] = 1ll * i * (i - 1) / 2;
        for (int j = 3; j <= mx_strength; j ++)
            combnum[j][i] = combnum[j][i - 1] + combnum[j - 1][i - 1];
    }

    // long long M = combnum[strength - 1][nvar] * (1ll << strength) + 2;
    // std::cout << "M = " << M << std::endl;
    // covered_now_strength_bitmap = new MyBitSet[M + 5];
    // for (int i = 0; i <= M; i ++)
    //     covered_now_strength_bitmap[i] = MyBitSet(nvar);
    
    // for (int i = 0; i < group_num; i ++)
    //     group_flag[i] = 0;
    // for (int i = 0; i < thread_num; i ++)
    //     for (int j = 0; j < group_num; j ++)
    //         group_flag_array[i][j] = 0;
    
    // // set_covered_now_strength_bitmap
    // std::cout << "begin to set covered_now_strength_bitmap" << std::endl;
    // for (const MyBitSet &tc : CA)
    //     set_covered_now_strength_bitmap(tc, 1, -1, 0);

    long long M = (strength > 2 ? (combnum[strength - 1][nvar] * (1ll << strength)) : (nvar << 2)) + 2;
    // std::cout << "M = " << M << std::endl;
    covered_last_strength_bitmap = nullptr;
    covered_now_strength_bitmap = new MyBitSet[M + 5];
    for (int i = 0; i < M; i ++)
        covered_now_strength_bitmap[i] = MyBitSet(nvar);
    
    if (strength > 2) {
        for (int i = 0; i < group_num; i ++)
                group_flag[i] = 0;
        for (int i = 0; i < thread_num; i ++)
            for (int j = 0; j < group_num; j ++)
                group_flag_array[i][j] = 0;
        // std::cout << "begin to set covered_now_strength_bitmap" << std::endl;
        for (const MyBitSet &tc : CA)
            set_covered_now_strength_bitmap(tc, 1, -1, 0);
        // std::cout << "end\n";
    }
    
    else {
        for (auto testcase : CA)
            for (int i = 0; i < nvar; i ++) {
                if (testcase.get(i)) {
                    covered_now_strength_bitmap[(nvar + i) << 1 | 1] |= testcase;
                    covered_now_strength_bitmap[(nvar + i) << 1] %= testcase;
                }
                else {
                    covered_now_strength_bitmap[i << 1 | 1] |= testcase;
                    covered_now_strength_bitmap[i << 1] %= testcase;
                }
            }
    }
        

    //多线程设置
    cdcl_solver_array.resize(thread_num + 1);
    for (int i = 0; i < thread_num; i ++) {
        cdcl_solver_array[i] = new CDCLSolver::Solver;
        cdcl_solver_array[i]->read_clauses(nvar, clauses);
    }
    all_tuples = new u_int64_t[thread_num + 1]{0};
    invalid_nums = new u_int64_t[thread_num + 1]{0};
    covered_nums = new u_int64_t[thread_num + 1]{0};
    uncovered_tuples_array.resize(thread_num);
    if (need_all_valid_tuples)
        vaild_tuples_array.resize(thread_num);

    // cdcl_sampler = new ExtMinisat::SamplingSolver(nvar, clauses, seed, true, 0);
    cdcl_sampler_array.resize(thread_num + 1);
    for (int i = 0; i < thread_num; i ++)
        cdcl_sampler_array[i] = new ExtMinisat::SamplingSolver(nvar, clauses, seed + i, true, 0);
}

Expandor:: ~Expandor() {
    
}

vector<vector<int> > Expandor:: get_final_tc () {
    vector<vector<int> > final_tc;
    for (auto tc : old_testcase)
        final_tc.emplace_back(tc);
    for (auto tc : new_testcase)
        final_tc.emplace_back(tc);
    return final_tc;
}

void Expandor:: set_covered_now_strength_bitmap (const MyBitSet &tc, int idx, int last, t_tuple tuple) {

    if (idx == strength) {
        u_int64_t offset = TupleToIndex(strength - 1, tuple);
        covered_now_strength_bitmap[offset << 1 | 1] |= tc;
        covered_now_strength_bitmap[offset << 1] %= tc;
        return ;
    }
    for (int i = last + 1; i < nvar - (strength - idx) + 1; i ++) {   
        int gid = group_id[i];
        if (gid == -1 || (! group_flag[gid])) {
            if (~gid) {
                if (! tc.get(i))
                    continue ;
                group_flag[gid] = 1;
            }
            tuple.v[idx - 1] = tc.get(i) ? (i + 1) : -(i + 1);
            set_covered_now_strength_bitmap(tc, idx + 1, i, tuple);

            if (~gid)
                group_flag[gid] = 0;
        }
    }

}
bool Expandor:: check_part_invalid(t_tuple tp) {

    for (int ban = 0; ban < strength; ban ++) {
        t_tuple tmp(strength - 1);
        for (int i = 0, j = 0; i < strength; i ++) {
            if (i == ban)
                continue ;
            
            j ++;
            if (j == strength - 1) {
                int p = abs(tp.v[i]) - 1, v = tp.v[i] > 0;
                u_int64_t state = TupleToIndex(strength - 2, tmp);
                if (! covered_last_strength_bitmap[state << 1 | v].get(p))
                    return 1;
            }
            else
                tmp.v[j - 1] = tp.v[i]; 
        }
    }
    return 0;

}

bool Expandor:: check_clauses_invalid(t_tuple tp) {

    for (t_tuple t : t_clauses) {
        bool flag = true;
        for (int i = 0; i < strength; i ++)
            if (t.v[i] != -tp.v[i]) {
                flag = false; break ;
            }
        if (flag) return 1;
    }
    return 0;
}

void Expandor:: get_remaining_valid_tuples (int thread_id, int value, int idx, int last, int _ed, t_tuple tuple) {
    
    if (idx == strength) {
        
        for (int i = strength - 2, j = value; i >= 0; i --, j >>= 1)
            if (! (j & 1))
                tuple.v[i] = -tuple.v[i];
        
        u_int64_t offset = TupleToIndex(strength - 1, tuple);
        for (int i = last + 1; i < nvar; i ++) {
            
            int gid = group_id[i];
            if ((~gid) && group_flag_array[thread_id][gid])
                continue ;

            for (int v = gid == -1 ? 0 : 1; v <= 1; v ++) {
                
                all_tuples[thread_id] ++;
                tuple.v[idx - 1] = v ? (i + 1) : -(i + 1);

                if (covered_now_strength_bitmap[offset << 1 | v].get(i)) {
                    covered_nums[thread_id] ++;

                    if (need_all_valid_tuples)
                        vaild_tuples_array[thread_id].emplace_back(tuple);
                }
                else {
                    
                    if (use_invalid_expand && check_part_invalid(tuple))
                        invalid_nums[thread_id] ++;
                    else if (use_invalid_expand && check_clauses_invalid(tuple))
                        invalid_nums[thread_id] ++;
                    else {
                        for (int _ = 0; _ < strength; _ ++) {
                            int j = abs(tuple.v[_]) - 1, vj = tuple.v[_] > 0;
                            cdcl_solver_array[thread_id]->add_assumption(j, vj);
                        }

                        bool res = cdcl_solver_array[thread_id]->solve();
                        cdcl_solver_array[thread_id]->clear_assumptions();
                        
                        if (res) {
                            t_tuple tep(strength);
                            for (int i = 0; i < strength; i ++)
                                tep.v[i] = tuple.v[i];
                            uncovered_tuples_array[thread_id].emplace_back(tep);
                            if (need_all_valid_tuples) 
                                vaild_tuples_array[thread_id].emplace_back(tep);
                        }

                        else invalid_nums[thread_id] ++;
                    }
                }
            }
        }
        return ;
    }

    int ed = std::min(nvar - (strength - idx) + 1, _ed);
    for (int i = last + 1; i < ed; i ++) {
        
        int gid = group_id[i];
        if (gid == -1 || (! group_flag_array[thread_id][gid])) {
            if (~gid) {
                if (! (value & (1 << (strength - idx - 1))))
                    continue ;
                group_flag_array[thread_id][gid] = 1;
            }
            tuple.v[idx - 1] = i + 1;
            get_remaining_valid_tuples(thread_id, value, idx + 1, i, nvar, tuple);

            if (~gid)
                group_flag_array[thread_id][gid] = 0;
        }

    }
}

void Expandor:: thread_work(int thread_id, int st_value, int ed_value, int st_dfs, int ed_dfs) {

    for (int value = st_value; value < ed_value; value ++)
        get_remaining_valid_tuples(thread_id, value, 1, st_dfs - 1, ed_dfs, t_tuple(strength));
    
}

void Expandor:: cdcl_sampler_thread_work(int thread_id, const vector<pair<int, int> > &prob, int st, int ed) {

    for (int i = st; i < ed; i ++) {
        cdcl_sampler_array[thread_id]->set_prob(prob);
        cdcl_sampler_array[thread_id]->get_solution(candidate_testcase_set_[i]);
    }

}

void Expandor:: Expand() {

    strength ++;
    
    t_clauses.clear();
    for (const vector<int> &c : clauses)
        if (c.size() == strength) {
            t_tuple tep(c);
            t_clauses.emplace_back(tep);
        }

    for (vector<int> tc : new_testcase)
        old_testcase.emplace_back(tc);
    new_testcase.clear();

    delete [] covered_last_strength_bitmap;
    covered_last_strength_bitmap = covered_now_strength_bitmap;

    num_combination_all_possible_ = combnum[strength][nvar] * (1ll << strength);
    std::cout << "num_combination_all_possible_ = " << num_combination_all_possible_ << std::endl;
    
    long long M = combnum[strength - 1][nvar] * (1ll << strength) + 2;
    // std::cout << "M = " << M << std::endl;
    covered_now_strength_bitmap = new MyBitSet[M + 5];
    for (int i = 0; i <= M; i ++)
        covered_now_strength_bitmap[i] = MyBitSet(nvar);
    
    for (int i = 0; i < group_num; i ++)
        group_flag[i] = 0;
    for (int i = 0; i < thread_num; i ++)
        for (int j = 0; j < group_num; j ++)
            group_flag_array[i][j] = 0;
    
    // set_covered_now_strength_bitmap
    // std::cout << "begin to set covered_now_strength_bitmap" << std::endl;
    for (const MyBitSet &tc : CA)
        set_covered_now_strength_bitmap(tc, 1, -1, 0);
    
    // 枚举
    // std::cout << "begin to get remaining valid tuples" << std::endl;
    uncovered_tuples.clear();
    if (need_all_valid_tuples) {
        vaild_tuples_array.clear();
        vaild_tuples_array.resize(thread_num);
    }
    // get_remaining_valid_tuples(1, -1, 0);
    for (int i = 0; i < thread_num; i ++) {
        all_tuples[i] = invalid_nums[i] = covered_nums[i] = 0;
        uncovered_tuples_array[i].clear();
        if (need_all_valid_tuples)
            vaild_tuples_array[i].clear();
    }
    
    Threads.clear();
    int offset = 0, mx_value = (1 << (strength - 1)) - 1;
    if (thread_num <= mx_value + 1) {
        int every = (mx_value + thread_num) / thread_num;
        for (int i = 0; i < thread_num; i ++) {
            int ed_value = i == thread_num - 1 ? mx_value + 1 : offset + every;
            Threads.emplace_back(&Expandor:: thread_work, this, i, offset, ed_value, 0, nvar);
            offset += every;
        }
    }
    else {
        int need_div = (thread_num + mx_value) / (mx_value + 1);
        for (int i = 0, thread_id = 0; i <= mx_value; i ++) {
            int every = nvar / need_div;
            if (every < 1)
                every = 1;
            for (int j = 0, k = 0; k < need_div; j += every, k ++) {
                Threads.emplace_back(&Expandor:: thread_work, this, 
                thread_id ++, i, i + 1, j, k == need_div - 1 ? nvar : j + every);
            }
        }
    }

    for (auto &thread : Threads)
        thread.join();
    
    // std::cout << "meow\n";
    uncovered_tuples.clear();
    vaild_tuples.clear();
    for (int i = 0; i <= nvar; i ++) {
        count_each_var_uncovered[0][i] = 0;
        count_each_var_uncovered[1][i] = 0;
    }
    for (int i = 0; i < thread_num; i ++) {
        for (t_tuple t : uncovered_tuples_array[i]) {
            uncovered_tuples.emplace_back(t);
            for (int _ = 0; _ < strength; _ ++) {
                int j = abs(t.v[_]) - 1, vj = t.v[_] > 0;
                count_each_var_uncovered[vj][j] ++;
            }
        }
        if (need_all_valid_tuples)
            for (t_tuple t : vaild_tuples_array[i])
                vaild_tuples.emplace_back(t);
    }

    if (need_all_valid_tuples)
        std::vector<vector<t_tuple> >().swap(vaild_tuples_array);

    tuples_U = uncovered_tuples;
    
    u_int64_t sum_covered_nums = 0, sum_invalid_nums = 0, sum_all_tuples = 0;
    for (int i = 0; i < thread_num; i ++) {
        sum_covered_nums += covered_nums[i];
        sum_invalid_nums += invalid_nums[i];
        sum_all_tuples += all_tuples[i];
    }
    
    // std::cout << "success" << std::endl;
    std::cout << "covered valid " << strength <<"-wise tuple nums: " << sum_covered_nums << "\n";
    std::cout << "uncovered valid " << strength << "-wise tuple nums: " << (uncovered_nums = uncovered_tuples.size()) << "\n";
    std::cout << "all valid " << strength <<"-wise tuple nums: " << sum_covered_nums + uncovered_nums << "\n";
    std::cout << "invalid " << strength << "-wise tuple nums: " << sum_invalid_nums << "\n";
    std::cout << "all tuples: " << sum_all_tuples << std::endl;

    // 其他设置
    if (use_cache) {
        candidate_testcase_set_.resize(2 * candidate_set_size);
        gain.resize(2 * candidate_set_size);
        vector<pair<int, int> > prob;
        prob.reserve(nvar);
        for (int i = 0; i < nvar; i ++) {
            int v1 = count_each_var_uncovered[1][i];
            int v2 = count_each_var_uncovered[0][i];
            prob.emplace_back(make_pair(v1, v2));
        }

        // for (int i = 0; i < candidate_set_size; i ++) {
        //     cdcl_sampler->set_prob(prob);
        //     cdcl_sampler->get_solution(candidate_testcase_set_[i]);
        // }

        Threads.clear();
        int every = (candidate_set_size + thread_num - 1) / thread_num;
        for (int i = 0, j = 0; i < thread_num; i ++, j += every) {
            Threads.emplace_back(&Expandor:: cdcl_sampler_thread_work, this, 
                i, prob,
                j, std::min(j + every, candidate_set_size));
        }
        for (auto &thread : Threads)
            thread.join();
    }
    else {
        gain.resize(candidate_set_size);
        candidate_testcase_set_.resize(candidate_set_size);
    }

    // 生成CA
    GenerateCoveringArray();
}

void Expandor::GenerateCandidateTestcaseSet() {
    
    vector<pair<int, int> > prob;
    prob.reserve(nvar);

    for (int i = 0; i < nvar; i ++) {
        int v1 = count_each_var_uncovered[1][i];
        int v2 = count_each_var_uncovered[0][i];
        prob.emplace_back(make_pair(v1, v2));
    }

    if (use_cache) {
        // for (int i = candidate_set_size; i < 2 * candidate_set_size; i ++) {
        //     cdcl_sampler->set_prob(prob);
        //     cdcl_sampler->get_solution(candidate_testcase_set_[i]);
        // }

        Threads.clear();
        int every = (candidate_set_size + thread_num - 1) / thread_num;
        for (int i = 0, j = candidate_set_size; i < thread_num; i ++, j += every) {
            Threads.emplace_back(&Expandor:: cdcl_sampler_thread_work, this, 
                i, prob,
                j, std::min(j + every, 2 * candidate_set_size));
        }
        for (auto &thread : Threads)
            thread.join();
    }
    else {
        // for (int i = 0; i < candidate_set_size; i ++) {
        //     cdcl_sampler->set_prob(prob);
        //     cdcl_sampler->get_solution(candidate_testcase_set_[i]);
        // }

        Threads.clear();
        int every = (candidate_set_size + thread_num - 1) / thread_num;
        for (int i = 0, j = 0; i < thread_num; i ++, j += every) {
            Threads.emplace_back(&Expandor:: cdcl_sampler_thread_work, this, 
                i, prob,
                j, std::min(j + every, candidate_set_size));
        }
        for (auto &thread : Threads)
            thread.join();
    }
    
}

int Expandor:: get_gain(const vector<int> &testcase) {
    
    int score = 0;
    for (t_tuple t : uncovered_tuples)
        if (is_covered(testcase, t))
            score ++;
    return score;

}

void Expandor:: get_gain_thread(int st, int ed) {

    for (int i = st; i < ed; i ++) {
        gain[i].first = get_gain(candidate_testcase_set_[i]);
        gain[i].second = i;
    }

}

int Expandor:: SelectTestcaseFromCandidateSetByTupleNum() {
    
    if(use_cache) {

        // for (int i = 0; i < 2 * candidate_set_size; i ++) {
        //     gain[i].first = get_gain(candidate_testcase_set_[i]);
        //     gain[i].second = i;
        // }
        int every = (2 * candidate_set_size + thread_num - 1) / thread_num;
        Threads.clear();
        for (int i = 0; i < 2 * candidate_set_size; i += every)
            Threads.emplace_back(&Expandor:: get_gain_thread, this, 
                i, std::min(i + every, 2 * candidate_set_size));
        for (auto &thread : Threads)
            thread.join();

        sort(gain.begin(), gain.end());

        if(gain[2 * candidate_set_size - 1].first <= 0)
            return -1;
        
        vector<vector<int> > tmp;
        tmp.resize(candidate_set_size);
        for (int i = candidate_set_size; i < 2 * candidate_set_size; i ++)
            tmp[i - candidate_set_size] = candidate_testcase_set_[gain[i].second];
        for(int i = 0; i < candidate_set_size; i ++)
            candidate_testcase_set_[i] = tmp[i];
        return candidate_set_size - 1;
    }

    else {
        // int mx = get_gain(candidate_testcase_set_[0]), mxi = 0;
        // for (int i = 1; i < candidate_set_size; i ++) {
        //     int score = get_gain(candidate_testcase_set_[i]);
        //     if(score > mx)
        //         mx = score, mxi = i;
        // }
        // return mx > 0 ? mxi : -1;

        int every = (candidate_set_size + thread_num - 1) / thread_num;
        Threads.clear();
        for (int i = 0; i < candidate_set_size; i += every)
            Threads.emplace_back(&Expandor:: get_gain_thread, this, 
                i, std::min(i + every, candidate_set_size));
        for (auto &thread : Threads)
            thread.join();

        int mx = gain[0].first, mxi = gain[0].second;
        for (int i = 1; i < candidate_set_size; i ++)
            if (gain[i].first > mx)
                mx = gain[i].first, mxi = gain[i].second;
        return mx > 0 ? mxi : -1;
    }
}

int Expandor:: GenerateTestcase() {
    GenerateCandidateTestcaseSet();
    return SelectTestcaseFromCandidateSetByTupleNum();
}

void Expandor:: Update_t_TupleInfo(const vector<int> &testcase, bool setmap) {
    
    vector<t_tuple> tep;
    for (t_tuple t : uncovered_tuples) {
        
        if (! is_covered(testcase, t)) {
            tep.emplace_back(t);
            continue ;
        }

        for (int i = 0; i < strength; i ++) {
            int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
            count_each_var_uncovered[vi][pi] --;
        }
    
        u_int64_t state = TupleToIndex(strength - 1, t);
        int p = abs(t.v[strength - 1]) - 1, v = t.v[strength - 1] > 0;
        if (setmap)
            covered_now_strength_bitmap[state << 1 | v].set(p);
    }
    uncovered_tuples = tep;
    uncovered_nums = uncovered_tuples.size();

}

void Expandor::GenerateCoveringArray() {

    int old_size = old_testcase.size();
    for (int num_generated_testcase_ = 1; ; num_generated_testcase_ ++) {
        
        int idx = GenerateTestcase();
        if(idx == -1)
            break ;
        
        new_testcase.emplace_back(candidate_testcase_set_[idx]);
        Update_t_TupleInfo(candidate_testcase_set_[idx], true);

        std::cout << "\033[;32mc current test suite size: " 
            << old_size + num_generated_testcase_ 
            << ", current uncovered valid " << strength << "-wise tuple nums: " 
            << uncovered_nums 
            << " \033[0m" << std::endl;
    }
    ReplenishTestCase();

    // 更新 CA
    for (const vector<int> &tc : new_testcase) {
        MyBitSet tep(nvar);
        for (int i = 0; i < nvar; i ++)
            if (tc[i])
                tep.set(i);
        CA.emplace_back(tep);
    }
}

void adaptive_adjustment () {
    
}

void Expandor:: Update_t_TupleInfo(int st, const vector<int> &testcase, const vector<int> &sidx) {

    int sz = uncovered_tuples.size();
    for (int p = st; p < sz; p ++) {

        t_tuple t = uncovered_tuples[sidx[p]];

        if (! is_covered(testcase, t))
            continue ;
        
        u_int64_t state = TupleToIndex(strength - 1, t);
        int pl = abs(t.v[strength - 1]) - 1, vl = t.v[strength - 1] > 0;
        if (covered_now_strength_bitmap[state << 1 | vl].get(pl))
            continue ;

        uncovered_nums --;
        covered_now_strength_bitmap[state << 1 | vl].set(pl);
        for (int i = 0; i < strength; i ++) {
            int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
            count_each_var_uncovered[vi][pi] --;
        }
    }

}

void Expandor:: ReplenishTestCase() {

    std::cout << "uncovered_nums: " << uncovered_nums << "\n";
    std::cout << "add new testcase: " << new_testcase.size() << "\n";

    vector<int> sidx(uncovered_nums, 0);
    std::iota(sidx.begin(), sidx.end(), 0);
    std::shuffle(sidx.begin(), sidx.end(), gen);

    int old_size = old_testcase.size(), sz = uncovered_nums;
    for (int p = 0; p < sz; p ++) {

        t_tuple t = uncovered_tuples[sidx[p]];

        if (uncovered_nums <= 0)
            break ;

        u_int64_t state = TupleToIndex(strength - 1, t);
        int pl = abs(t.v[strength - 1]) - 1, vl = t.v[strength - 1] > 0;
        if (covered_now_strength_bitmap[state << 1 | vl].get(pl))
            continue ;

        for (int i = 0; i < strength; i ++) {
            int pi = abs(t.v[i]) - 1, vi = t.v[i] > 0;
            cdcl_sampler_array[0]->add_assumption(pi, vi);
        }
        
        vector<pair<int, int> > prob;
        prob.resize(nvar);
        for (int i = 0; i < nvar; i ++) {
            int v1 = count_each_var_uncovered[1][i];
            int v2 = count_each_var_uncovered[0][i];
            prob.emplace_back(make_pair(v1, v2));
        }
        cdcl_sampler_array[0]->set_prob(prob);

        vector<int> tep(nvar, 0);
        cdcl_sampler_array[0]->get_solution(tep);
        cdcl_sampler_array[0]->clear_assumptions();

        new_testcase.emplace_back(tep);
        Update_t_TupleInfo(p, tep, sidx);

        std::cout << "\033[;32mc current test suite size: " 
        <<  old_size + new_testcase.size()
        << ", current uncovered valid " << strength << "-wise tuple nums: " << uncovered_nums
        << " \033[0m" << std::endl;
    }
    std::cout << "uncovered_nums: " << uncovered_nums << "\n";

}

void Expandor:: SaveTestcaseSet(string result_path) {

    ofstream res_file(result_path);

    // std::cout << old_testcase.size() << "\n";
    for (const vector<int>& testcase: old_testcase) {
        for (int v = 0; v < nvar; v++)
            res_file << testcase[v] << " ";
        res_file << "\n";
    }

    // std::cout << new_testcase.size() << "\n";
    for (const vector<int>& testcase: new_testcase) {
        for (int v = 0; v < nvar; v++)
            res_file << testcase[v] << " ";
        res_file << "\n";
    }
    res_file.close();

    // std::cout << "c Testcase set saved in " << result_path << endl;

}

void Expandor:: SaveTestcaseSet() {
    SaveTestcaseSet(output_testcase_path);
}

void Expandor:: Free() {
    delete [] covered_last_strength_bitmap;
    delete [] covered_now_strength_bitmap;
    for (int i = 0; i < thread_num; i ++) {
        delete cdcl_solver_array[i];
        delete cdcl_sampler_array[i];
    }
}