#include <cstring>
#include <iostream>
#include <string>
using namespace std;
using std::string;

struct Argument {

    // 基础
    string input_cnf_path;
	string reduced_cnf_file_path;
	string init_CA_file_path;
    string output_testcase_path;
    string group_file_path;

    bool flag_input_cnf_path;
    bool flag_init_CA_file_path;
    bool flag_output_testcase_path;
    bool flag_group_file_path;

    // 共用
    int seed, strength, thread_num;
    int use_cnf_reduction;

    // Expandor
    int use_cache, use_invalid_expand, init_strength, candidate_set_size;
    int use_RALS;

    // Optimizer
    int stop_length;
    int forced_greedy_percent;
    int sampling_num, taboo_method;
    int use_group, use_weight;
    int opt_method;
};