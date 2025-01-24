"""
This file is used for translating CNFs (in Dimacs format) into the format
accepted by AutoCCAG, FastCA and TCA. 
"""

import os, sys, copy
import uuid

class FormatEncoder:
    def __init__(self, cnf_instance_path, t_wise_strength, model_file_path, constraints_file_path):
        self.t_wise_strength = t_wise_strength
        self.num_vars = 0
        self.num_clauses = 0
        self.clauses = []
        self.cnf_instance_path = cnf_instance_path
        self.model_file_path = model_file_path
        self.constraints_file_path = constraints_file_path
        self._reading_cnf_instance()
        return
    
    def encoding(self):
        self._generate_model_file()
        self._generate_constraints_file()
        return
    
    def clear(self):
        os.system("rm %s %s" % (self.model_file_path, self.constraints_file_path))

    def _generate_model_file(self):
        fout = open(self.model_file_path, 'w+')
        fout.write("%d\n" % self.t_wise_strength)
        fout.write("%d\n" % self.num_vars)
        line = "2 " * self.num_vars
        line = line.strip()
        fout.write("%s\n" % line)
        fout.close()
        return
    
    def _generate_constraints_file(self):
        fout = open(self.constraints_file_path, 'w+')
        fout.write("%d\n" % self.num_clauses)
        for clause in self.clauses:
            line = ""
            for literal in clause:
                if literal < 0:
                    value = 2*(-literal) - 1
                else:
                    value = 2*literal - 2
                line += "- %d " % value
            line = line.strip()
            fout.write("%d\n" % len(clause))
            fout.write("%s\n" % line)
        fout.close()
        return
    
    def _reading_cnf_instance(self):
        fin = open(self.cnf_instance_path, 'r')
        while True:
            line = fin.readline()
            if not line:
                break
            words = line.strip().split()
            if len(words) <= 0:
                continue
            elif (words[0] =='c') or (words[0][0] == 'c'):
                continue
            elif words[0] == 'p':
                self.num_vars = int(words[2])
                self.num_clauses = int(words[3])
                continue
            else:
                clause = []
                for word in words:
                    literal = int(word)
                    if literal != 0:
                        clause.append(literal)
                self.clauses.append(copy.deepcopy(clause))
                continue
        fin.close()
        return

if __name__ == "__main__":

    cnf_path, strength = sys.argv[1], int(sys.argv[2])
    model_path, constr_path = sys.argv[3], sys.argv[4]

    formatEncoder = FormatEncoder(cnf_path, strength, model_path, constr_path)
    formatEncoder.encoding()

    # python3 formatencoding.py _Input/Input/spins_.cnf 3 spins_boolean.model spins_boolean.constraints
