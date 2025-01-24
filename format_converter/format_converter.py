import os, sys
import secrets

def generate_random_hex(n):
    return secrets.token_hex(n)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 format_converter.py modelFile constrFile cnfFile groupFile")
        exit(0)

    modelFile = sys.argv[1]
    constrFile = sys.argv[2]
    cnfFile = sys.argv[3]
    groupFile = sys.argv[4]
    actsFile = f"tmp/{secrets.token_hex(32)}_acts.txt"
    command = f"./FormatConverter {modelFile} {constrFile} {actsFile}"
    
    os.system(command)

    ctwFile = f"tmp/{secrets.token_hex(32)}.ctw"
    command = f"python3 turn.py {actsFile} {ctwFile}"
    os.system(command)

    listfile = f"tmp/{secrets.token_hex(32)}_list.txt"
    command = f"./ctw_parser {ctwFile} {cnfFile} {listfile} 2> /dev/null"
    os.system(command)

    command = f"python3 get_group.py {cnfFile} {groupFile}"
    os.system(command)

    if os.path.exists(actsFile):
        os.remove(actsFile)
    if os.path.exists(ctwFile):
        os.remove(ctwFile)
    if os.path.exists(listfile):
        os.remove(listfile)
    
    
