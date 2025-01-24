#include "GlobalOptimizer.h"

#include <iostream>
#include <signal.h>
#include <chrono>
#include <cstring>
#include <string>
using namespace std;

void HandleInterrupt(int sig){
    std::cout << "c" << endl;
    std::cout << "c caught signal... exiting" << std::endl;

    exit(-1);
}

void SetupSignalHandler(){
    signal(SIGTERM, HandleInterrupt);
    signal(SIGINT, HandleInterrupt);
    signal(SIGQUIT, HandleInterrupt);
    signal(SIGKILL, HandleInterrupt);
}

bool ParseArgument(int argc, char **argv, Argument &argu) {

    // 基础
    // string input_cnf_path;
	// string reduced_cnf_file_path;
	// string init_CA_file_path;
    // string output_testcase_path;
    // string group_file_path;
    // bool flag_input_cnf_path;
    // bool flag_init_CA_file_path;
    // bool flag_output_testcase_path;
    // bool flag_group_file_path;
    argu.flag_input_cnf_path = argu.flag_init_CA_file_path = false;
    argu.flag_output_testcase_path = argu.flag_group_file_path = false;
    
    // 共用
    // int seed, strength, thread_num;
    // int use_cnf_reduction;
    argu.seed = 1;
    argu.strength = -1;
    argu.thread_num = 32;
    argu.use_cnf_reduction = 1;

    // Expandor
    // int use_cache, use_invalid_expand, init_strength, candidate_set_size;
    // int use_RALS;
    argu.use_cache = argu.use_invalid_expand = argu.use_RALS = 1;
    argu.init_strength = 2;
    argu.candidate_set_size = 100;

    // Optimizer
    // int stop_length;
    // int forced_greedy_percent;
    // int sampling_num, taboo_method;
    // int use_group, use_weight;
    // int opt_method;
    argu.stop_length = 10000;
    argu.forced_greedy_percent = 90;
    argu.sampling_num = 100;
    argu.taboo_method = 1;
    argu.use_group = 1;
    argu.opt_method = 1;

    for (int i = 1; i < argc; i ++) {
        
        // 基础
        if (strcmp(argv[i], "-input_cnf_path") == 0) {
			i ++;
			if(i >= argc) return false;
			argu.input_cnf_path = argv[i];
            argu.flag_input_cnf_path = true;
			continue ;
		}

        else if (strcmp(argv[i], "-group_file_path") == 0) {
			i ++;
			if(i >= argc) return false;
			argu.group_file_path = argv[i];
            argu.flag_group_file_path = true;
			continue ;
		}

        else if (strcmp(argv[i], "-init_CA_file_path") == 0) {
			i ++;
			if(i >= argc) return false;
			argu.init_CA_file_path = argv[i];
            argu.flag_init_CA_file_path = true;
			continue ;
		}

        else if (strcmp(argv[i], "-output_testcase_path") == 0) {
			i ++;
			if(i >= argc) return false;
			argu.output_testcase_path = argv[i];
            argu.flag_output_testcase_path = true;
			continue ;
		}

        // 共用
        // int seed, strength, thread_num;
        // int use_cnf_reduction;
        else if (strcmp(argv[i], "-seed") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.seed);
			continue ;
		}

        else if (strcmp(argv[i], "-strength") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.strength);
			continue ;
		}

        else if (strcmp(argv[i], "-thread_num") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.thread_num);
			continue ;
		}

        // Expandor
        // int use_cache, use_invalid_expand, init_strength, candidate_set_size;
        // int use_RALS;
        else if (strcmp(argv[i], "-use_cache") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.use_cache);
			continue ;
		}

        else if (strcmp(argv[i], "-use_FID") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.use_invalid_expand);
			continue ;
		}

        else if (strcmp(argv[i], "-init_strength") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.init_strength);
			continue ;
		}

        else if (strcmp(argv[i], "-candidate_set_size") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.candidate_set_size);
			continue ;
		}

        else if (strcmp(argv[i], "-use_RALS") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.use_RALS);
			continue ;
		}

        // Optimizer
        // int stop_length;
        // int forced_greedy_percent;
        // int sampling_num, taboo_method;
        // int use_group, use_weight;
        // int opt_method;
        else if (strcmp(argv[i], "-L") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.stop_length);
			continue ;
		}
        else if (strcmp(argv[i], "-forced_greedy_percent") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.forced_greedy_percent);
			continue ;
		}
        else if (strcmp(argv[i], "-sampling_num") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.sampling_num);
			continue ;
		}
        else if (strcmp(argv[i], "-taboo_method") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.taboo_method);
			continue ;
		}
        else if (strcmp(argv[i], "-use_group") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.use_group);
			continue ;
		}
        else if (strcmp(argv[i], "-use_weight") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.use_weight);
			continue ;
		}
        else if (strcmp(argv[i], "-opt_method") == 0) {
			i ++;
			if(i >= argc) return false;
			sscanf(argv[i], "%d", &argu.opt_method);
			continue ;
		}
    }

    if(! argu.flag_input_cnf_path || ! argu.flag_init_CA_file_path 
        || ! argu.flag_group_file_path || argu.strength == -1) {
            
            std::cout << argu.flag_input_cnf_path << " "
                << argu.flag_init_CA_file_path << " "
                << argu.flag_group_file_path << " "
                << argu.strength << "\n";
        
            return false;
        }

    if(! argu.flag_output_testcase_path) {
        int pos = argu.input_cnf_path.find_last_of( '/' );
        string cnf_file_name = argu.input_cnf_path.substr(pos + 1);
	    cnf_file_name.replace(cnf_file_name.find(".cnf"), 4, "");

        argu.output_testcase_path = cnf_file_name + "_testcase.out";
    }
    
    return true;
}

int main(int argc, char **argv) {
    
    SetupSignalHandler();
    
    Argument argu;

	if (! ParseArgument(argc, argv, argu)) {
		std::cout << "c Argument Error!" << std::endl;
		return -1;
	}

    Expandor expandor(argu);
    std::cout << "Expandor init success" << endl;

    if (! argu.use_RALS)
        argu.opt_method = 1;
    
    std::cout << "opt_method = " << argu.opt_method << "\n";
    if (argu.opt_method == 1) {
        for (int i = argu.init_strength + 1; i <= argu.strength; i ++) {
            std::cout << "------------------------ strengh = " << i
                        << " ------------------------" << std::endl;
            expandor.Expand();
        }

        if((! argu.use_RALS) || (! expandor.new_testcase.size())) {
            expandor.SaveTestcaseSet();
            std::cout << "End" << std::endl;
            return 0;
        }
        
        // expandor.SaveTestcaseSet();

        Optimizer optimizer(expandor, argu, false);
        expandor.Free();

        optimizer.search();
        optimizer.SaveTestcaseSet();
    }

    else if (argu.opt_method == 2) {

        for (int i = argu.init_strength + 1; i <= argu.strength; i ++) {
            std::cout << "------------------------ strengh = " << i
                        << " ------------------------" << std::endl;
            expandor.Expand();
        }

        GlobalOptimizer global_optimizer(expandor, expandor.get_final_tc(), argu);
        expandor.Free();
        global_optimizer.search();
        global_optimizer.SaveTestcaseSet();
    }

    else if (argu.opt_method == 22) {
        
        for (int i = argu.init_strength + 1; i <= argu.strength; i ++) {
            std::cout << "------------------------ strengh = " << i
                        << " ------------------------" << std::endl;
            if (i == argu.strength)
                expandor.need_all_valid_tuples = true;
            expandor.Expand();
        }

        Optimizer optimizer(expandor, argu, true);
        expandor.Free();

        optimizer.search();
        optimizer.SaveTestcaseSet();
    }

    else if (argu.opt_method == 3) {
        
        for (int i = argu.init_strength + 1; i <= argu.strength; i ++) {
            std::cout << "------------------------ strengh = " << i
                        << " ------------------------" << std::endl;
            expandor.Expand();
        }

        int orgL = argu.stop_length;
        argu.stop_length = 5000;

        Optimizer optimizer(expandor, argu, false);
        expandor.Free();
        optimizer.search();
        
        vector<vector<int> > init_tc = optimizer.get_final_tc();
        optimizer.Free();

        argu.stop_length = orgL;
        GlobalOptimizer global_optimizer(expandor, init_tc, argu);
        global_optimizer.search();
        global_optimizer.SaveTestcaseSet();
    }
    
    std::cout << "End\n";
    return 0;
}