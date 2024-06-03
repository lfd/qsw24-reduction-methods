import sys
import os
import shutil 
from ExperimentClass import *

def executeExp(nm, nj, qis, quarkMethod, path):
    if not os.path.exists(path):
            os.makedirs(path)

    shutil.copy2("instance.py", path)
    shutil.copy2("ExperimentClass.py", path)
    shutil.copy2("ising_part_model.py", path)
    shutil.copy2("Transform_dict_arrays.py", path)


    exp = ExperimentClass("ExperimentRigClasses.csv", path, quarkMethod)        
    exp.execute(int(nm), int(nj), bool((qis.lower() in ['true', 'jep', 'certainly', 'yes', 'yeah', 'yup', 'uh-huh'])))

if __name__ == "__main__":
    nm, nj, qis = ["","",""]
    exec = True
    
    if len(sys.argv) == 3:
        nm = sys.argv[1]
        nj = sys.argv[2]
        qis = "false"
        quark = "false"
    elif len(sys.argv) == 4:
        nm = sys.argv[1]
        nj = sys.argv[2]
        qis = sys.argv[3]
        quark = "false"
    elif len(sys.argv) == 5:
        nm = sys.argv[1]
        nj = sys.argv[2]
        qis = sys.argv[3]
        quark = sys.argv[4]
    else:
        exec = False
        print("You must specify correct parameters: Number of machines, Number of Jobs, [Qiskit Transpile]: True or False, [useAllQuarkReducMethods]: True or False")
    
    
    
    if exec:    
        path = "Erg/" + str(int(nm)) + "m_" + str(int(nj)) + "j_" + str(bool((qis.lower() in ['true', 'jep', 'certainly', 'yes', 'yeah', 'yup', 'uh-huh'])))
        
        executeExp(nm, nj, qis, "stupid_but_supposedly_fastest", path + "_Stupid" + "/")
        
        if (bool((quark.lower() in ['true', 'jep', 'certainly', 'yes', 'yeah', 'yup', 'uh-huh']))):
            executeExp(nm, nj, qis, "probably_better", path + "_Better" + "/")
            executeExp(nm, nj, qis, "supposedly_best_but_slowest", path + "_Best" + "/")
        
      
