import os, sys

def get_info(model_path, constr_path):
    with open(model_path, "r") as model_file:
        line = model_file.readline()
        line = model_file.readline()
    nvar = int(line.strip())
    with open(constr_path, "r") as constr_file:
        line = constr_file.readline()
    ncon = int(line.strip())
    return [nvar, ncon]

if __name__ == "__main__":

    output_file = open("benchmark_information.csv", "w")
    
    output_file.write("Instance,#options,#constraints")
    output_file.write("\nreal-world:")
    for name in ["apache", "bugzilla", "gcc", "spins", "spinv"]:
        model_path = f"benchmarks/real-world/{name}_2wise.model"
        constr_path = f"benchmarks/real-world/{name}.constraints"
        info = get_info(model_path, constr_path)
        output_file.write(f"\n{name},{info[0]},{info[1]}")
    output_file.write("\nsynthetic:")
    for idx in range(1, 31):
        model_path = f"benchmarks/synthetic/Syn_{idx}_2wise.model"
        constr_path = f"benchmarks/synthetic/Syn_{idx}.constraints"
        info = get_info(model_path, constr_path)
        output_file.write(f"\nSyn_{idx},{info[0]},{info[1]}")
    output_file.close()