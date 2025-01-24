import os, sys

def turnCA(cnf_path, model_path, boolean_CA_path, output_path):
    # cnf_path, model_path, 

    # cnf_path = f"Input/cnf/{index_name}.cnf"
    with open(cnf_path, "r") as cnf_file:
        lines = cnf_file.readlines()
    
    name_to_index = {}
    for line in lines:
        arr = line.split()
        if len(arr) == 4 and arr[2] == "-->":
            name_to_index[arr[1]] = int(arr[3]) - 1
    
    # model_path = f"Input/3way/{index_name}.model"
    with open(model_path, "r") as model_file:
        lines = model_file.readlines()
    mvar = int(lines[1].replace("\n", ""))
    values = list(map(int, lines[2].split()))

    # boolean_CA_path = f"Results/seed{seed}/CA/{index_name}_CA.out"
    # output_path = f"Results/seed{seed}/initCA/{index_name}_CA.out"
    output_file = open(output_path, "w")
    with open(boolean_CA_path, "r") as boolean_CA_file:
        lines = boolean_CA_file.readlines()

    for line in lines:
        arr = list(map(int, line.split()))
        tc = [0 for i in range(mvar)]
        for i in range(mvar):
            if f"p{i}" in name_to_index:
                tc[i] = arr[name_to_index[f"p{i}"]]
            else:
                for v in range(values[i]):
                    idx = name_to_index[f"p{i}@v{v}"]
                    if arr[idx] == 1:
                        tc[i] = v
                        break
        for i in range(mvar):
            output_file.write(str(tc[i]) + " ")
        output_file.write("\n")
    output_file.close()

if __name__ == '__main__':

    # print("QAQ")
    if len(sys.argv) != 5:
        print("args Error !")
        exit(0)

    turnCA(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

    # python3 turnCA.py _Input/Input/bugzilla_.cnf _Input/models/benchmark_bugzilla_6.model tmp/bugzilla_t6_seed1_CA.out_boolean_CA.out tmp/bugzilla_t6.out