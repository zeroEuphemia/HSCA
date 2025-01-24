import os, sys

def get_group(cnf_path, output_path):
    
    with open(cnf_path, "r") as f:
        lines = f.readlines()
    
    member = {}
    for line in lines:
        arr = line.split()
        if len(arr) != 4:
            break
        if arr[0] != "c" or arr[1] == "@top" or arr[1].count("@") == 0:
            continue
        
        p, v = arr[1].split("@")
        p, v = int(p[1:]), int(v[1:])
        if p not in member:
            member[p] = []
        member[p].append(int(arr[3]))
    
    # for group in member:
    #     print(group, member[group])
    
    output_file = open(output_path, "w")
    output_file.write(str(len(member)) + "\n")
    for group in member:
        arr = member[group]
        output_file.write(str(len(arr)) + "\n")
        for x in arr:
            output_file.write(str(x) + " ")
        output_file.write("\n")
    output_file.close()
        
if __name__ == "__main__":
    
    # get_input("spins.cnf", "spins_group.txt")
    if len(sys.argv) != 3:
        print("Usage: python3 get_group.py cnfFile groupFile")
        exit(0)
    get_group(sys.argv[1], sys.argv[2])
