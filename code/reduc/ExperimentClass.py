import quark
import pandas as pd
import numpy as np

import os 

import time 
import scipy.special

from pathos.multiprocessing import ProcessPool as Pool

import csv
import json


import networkx as nx
import pprint

from qiskit import QuantumCircuit, transpile, Aer
from qiskit.providers.fake_provider import ConfigurableFakeBackend
from qiskit.providers.aer import AerSimulator, noise

from itertools import combinations

from ising_part_model import *

from util import *

class ExperimentClass:
    """
        class for executing the experiments
    """
    
    def __init__(self, csvName, resultsPath="", quarkReducMethod = "stupid_but_supposedly_fastest"):
        """
        Parameters
        ----------
            csvName: string 
                Name of the csv file in which the experiments data will be saved
            resultsPath: string 
                Path to where the results (csv inclusive) will be stored (must be "/" terminated)
            quarkReducMethod: string
                possible values are "stupid_but_supposedly_fastest", "probably_better", "supposedly_best_but_slowest"
        """
        self.csvName = csvName
        self.resultsPath = resultsPath
        self.quarkReducMethod = quarkReducMethod
        self.csvColumns = ["NoMachines", "NoJobs", "ElapsedTimeS", "ElapsedTimeSModelGen", "Difficulty", "NoVariablesBefore", "NoVariablesAfter",
 "NoTermsDeg1", "NoTermsDeg2", "NoTermsDeg3", "NoTermsDeg4", "EstimatedCircDepth", "QiskitDepth", "QiskitDepthTranspiled", "avgNodeDegree", "num_edges", "num_multi_edges",
 "NoTermsDegr1", "NoTermsDegr2", "NoTermsDegr3", "NoTermsDegr4", "EstimatedCircDepthr", "QiskitDepthr", "QiskitDepthTranspiledr", "avgNodeDegreer", "num_edgesr", "num_multi_edgesr"]
        
        self.pp=pprint.PrettyPrinter(indent=4)

    def execute(self, MaxNumMachines, MaxNumJobs, useQiskitTranspile = False):
        """
        execute the specified experiments up to MaxNumJobs and MaxNumMachines with MaxNumJobs >= MaxNumMachines during execution
        Parameters
        ----------
            MaxNumMachines: integer
            MaxNumJobs: integer     
            useQiskitTranspile: boolean
                transpile the circuit with Qiskit to analyze depth
        """
        completeExp = []
        
        
        argsM = []
        argsJ = []
        argsDifficulty = []
        argsQiskit = []
        
        for nM in range(1,MaxNumMachines):
            for nJ in range(nM, MaxNumJobs):
                for diffi in list(["easy", "medium", "hard"]):
                    
                    
                    argsM.insert(0, nM)
                    argsJ.insert(0, nJ)
                    argsDifficulty.insert(0, diffi)
                    argsQiskit.insert(0, useQiskitTranspile)
         
        with Pool() as p:
            completeExp = p.imap(self.executeSingle, argsM, argsJ, argsDifficulty, argsQiskit)
        
        with open(self.resultsPath + self.csvName, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(self.csvColumns)
         
            for singleExp in completeExp:   
                writer.writerow(list(singleExp))
        
        self.addDensities()
                
    def executeSingle(self, machines, jobs, difficulty, useQiskitTranspile):
        print("Currently in (m,j): ", machines, ", ", jobs)
        start=time.time()
        
        jobs_dict, machines_dict, rigs_matrix = ExperimentClass.genExample(jobs, machines, difficulty)
       
        dur = []
        for j in jobs_dict:
            dur.append(j["Duration"])
        
        A=1
        B=2*np.max(dur)
        C1=4*max(dur)**2
        lagrange_list = [A,B,C1]
        JSS_instance = JSS_model_part(jobs, machines, jobs_dict, machines_dict, rigs_matrix, lagrange_list, "min_makespan", "qubo", [],0)

        intermediate = time.time() - start
        
        quark_dict = JSS_instance.linear | JSS_instance.quadratic | JSS_instance.third | JSS_instance.fourth 
        
        polyBinary = quark.PolyBinary(quark_dict)     

        
        varBefore = len(polyBinary.variables)
        t1 = len(ExperimentClass.PolyPart(polyBinary, 1))
        t2 = len(ExperimentClass.PolyPart(polyBinary, 2))
        t3 = len(ExperimentClass.PolyPart(polyBinary, 3))
        t4 = len(ExperimentClass.PolyPart(polyBinary, 4))

        
        depth_estimate = t1 + 3*t2 + 5*t3 + 7*t4

        
        qcirc, qdepth = ExperimentClass.QiskitDepthEstimate(polyBinary)

        if(useQiskitTranspile):
            transcirc, transcircdepth = ExperimentClass.QiskitTranspile(qcirc)
        else:
            transcircdepth = -1
            
        
        G, Gdeg, Gavg, Gnum_edges, Gnum_multi_edges = ExperimentClass.genCouplingGraphFromPolynomial(polyBinary)

        
        
        reduced_poly, penaltyTerms= polyBinary.reduce(max_degree=2, reduction_variable_prefix='r', var_pair_choice=self.quarkReducMethod)
        penalty_poly = quark.PolyBinary(ExperimentClass.PenaltyDictToPolyDict(penaltyTerms))

        
        CombinedPoly = reduced_poly + penalty_poly
        varAfter = len(CombinedPoly.variables)
        t1r = len(ExperimentClass.PolyPart(CombinedPoly, 1)) 
        t2r = len(ExperimentClass.PolyPart(CombinedPoly, 2))
        t3r = len(ExperimentClass.PolyPart(CombinedPoly, 3))
        t4r = len(ExperimentClass.PolyPart(CombinedPoly, 4))

        
        qcircr, qdepthr = ExperimentClass.QiskitDepthEstimate(CombinedPoly)

        if(useQiskitTranspile):
            transcircr, transcircdepthr = ExperimentClass.QiskitTranspile(qcircr)
        else:
            transcircdepthr = -1
       
        
        depth_estimater = t1r + 3*t2r + 5*t3r + 7*t4r
        
        
        Gr, Gdegr, Gavgr, Gnum_edgesr, Gnum_multi_edgesr = ExperimentClass.genCouplingGraphFromPolynomial(CombinedPoly)
        
        end = time.time()
        elapsed = end - start
        
        
        pathExtended = self.resultsPath + str(machines) + "m_" + str(jobs) + "j_" + str(difficulty) + "_"
        
        with open(pathExtended + "polyBefore.json", "w") as pB:
            json.dump(polyBinary.__repr__(), pB) 
            
        with open(pathExtended + "polyAfter.json", "w") as pA:
            json.dump(CombinedPoly.__repr__(), pA) 
        
        
        with open(pathExtended + "GraphDegBefore.json", "w") as gdB:
            json.dump(str(Gdeg), gdB) 
            
        with open(pathExtended + "GraphDegAfter.json", "w") as gdA:
            json.dump(str(Gdegr), gdA) 
        
        
        with open(pathExtended + "GraphBefore.json", "w") as gB:
            json.dump(nx.node_link_data(G), gB) 
            
        with open(pathExtended + "GraphAfter.json", "w") as gA:
            json.dump(nx.node_link_data(Gr), gA) 
        
        
        """
        with open("7m_7j_medium_GraphBefore.json") as f:
        data = json.load(f)
        G = nx.node_link_graph(data)

        print(G.edges())

        nx.draw(G, with_labels = True)
        """
        
        return [machines, jobs, elapsed, intermediate, difficulty, varBefore, varAfter, 
               t1, t2, t3, t4, depth_estimate, qdepth, transcircdepth, Gavg, Gnum_edges, Gnum_multi_edges,
               t1r, t2r, t3r, t4r, depth_estimater, qdepthr, transcircdepthr, Gavgr, Gnum_edgesr, Gnum_multi_edgesr]
    
    def addDensities(self):
        """
        Add all possible density columns to the csv
        """
        dfp = pd.read_csv(self.resultsPath + self.csvName)
        ExperimentClass.AddDensities(dfp, "NoTermsDeg", "NoVariablesBefore", "density_deg", 4)
        ExperimentClass.AddDensities(dfp, "NoTermsDegr", "NoVariablesAfter", "density_degr", 4)
        dfp.to_csv(self.resultsPath + self.csvName)
        
    @staticmethod  
    def AddDensity(dataframe, noTermsCol, noVarsCol, nameNewCol, k):
        """
        Add a k-density column to a pandas dataframe
        """
        c1 = []
        for x in np.array(dataframe[[noTermsCol, noVarsCol]]):
            scp = (scipy.special.comb(x[1], k, exact=True))
            if (scp <= 0): 
                c1 += [0] 
            else:
                c1 += [x[0] / (scipy.special.comb(x[1], k, exact=True))]
            
        dataframe[nameNewCol] = c1
        
    @staticmethod
    def AddDensities(dataframe, noTermsCol, noVarsCol, nameNewCol, k):
        """
        Add densities up to degree k
        Expects enumerated string literals
        """
        for i in range(1,k+1):
            ExperimentClass.AddDensity(dataframe, noTermsCol + str(i), noVarsCol, nameNewCol + str(i), i)    
    
    @staticmethod
    def genExample(num_jobs, num_machines, difficulty):
        """
        Generates scalable examples
        
        Parameters
        ----------
            num_jobs: int
                Number of jobs
            num_machines: int
                Number of machines
            difficulty: string
                Valid forms are: "easy", "medium", "hard"
        """
         
        prng = np.random.default_rng(seed=42)
        
        
        
        
        rigsMatrix = np.zeros(shape=(num_jobs,num_jobs), dtype=np.int32)
        
        
        jobs = []
        
        
        machines = []
        
        
        if (difficulty == "easy"):
            rigsMatrix = np.ones(shape=(num_jobs,num_jobs), dtype=np.double)
            for i in range(0,num_jobs):
                rigsMatrix[i][i] = 0                
            
            for i in range(0, num_jobs):
                jobs.append({'Duration': 1, 'Rig': 0})
            
            for i in range(0,num_machines):
                machines.append({'Rig': 0})
                
        
        elif (difficulty == "medium"):
            rigsMatrix = prng.random((num_jobs, num_jobs))
            zerosetters = np.rint(prng.random((num_jobs, num_jobs))).astype(int)
            
            for i in range(0,num_jobs):
                rigsMatrix[i][i] = 0 
                for j in range(0, num_jobs):
                    if(zerosetters[i][j] == 0):
                        rigsMatrix[i][j] = 0
                               
            for i in range(0, num_jobs):
                jobs.append({'Duration': 1, 'Rig': i})
                
            for i in range(0,num_machines):
                machines.append({'Rig': 0})
        
        elif (difficulty == "hard"):
            job_durations = prng.random((num_jobs))
            
            rigsMatrix = prng.random((num_jobs, num_jobs))
            for i in range(0,num_jobs):
                rigsMatrix[i][i] = 0
                
            for i in range(0, num_jobs):
                jobs.append({'Duration': job_durations[i], 'Rig': i})
            
            for i in range(0,num_machines):
                machines.append({'Rig': i})
        
        else :
            print("Invalid difficulty in genExample(...)")
            
        
        rigsMatrix = (np.rint(np.array(rigsMatrix) * 20)).astype(int)
        
        for i in range(len(jobs)):
            jobs[i]['Duration'] = int(jobs[i]['Duration']*100)
            
        return jobs, machines, rigsMatrix
    
    
    
    @staticmethod
    def PolyPart(poly, degree):
        monomials = sorted(var_tuple for var_tuple in poly.keys() if len(var_tuple) == degree)
        return {var_tuple : poly.get_coefficient(*var_tuple) for var_tuple in monomials}
    
    
    @staticmethod
    def PenaltyDictToPolyDict(penaltyDict):
        out = dict()
        for key in penaltyDict:
            
            t = dict(penaltyDict[key])
            out.update(t)
        return out
    
    
    @staticmethod
    def QiskitDepthEstimate(polynomial):
        circ = QuantumCircuit(len(polynomial.variables))
        pC = polynomial.compact()
        for key in pC: 
            deg = len(key)
            if deg == 1:
                circ.rz(pC[key], key[0])
            elif deg == 2:
                circ.cx(key[0], key[1])
                circ.rz(pC[key], key[1])
                circ.cx(key[0], key[1])
            elif deg == 3:
                circ.cx(key[0], key[1])
                circ.cx(key[1], key[2])
                circ.rz(pC[key], key[2])
                circ.cx(key[1], key[2])
                circ.cx(key[0], key[1])
            elif deg == 4:
                circ.cx(key[0], key[1])
                circ.cx(key[1], key[2])
                circ.cx(key[2], key[3])
                circ.rz(pC[key], key[3])
                circ.cx(key[2], key[3])
                circ.cx(key[1], key[2])
                circ.cx(key[0], key[1])

        return circ, circ.depth()

    
    @staticmethod
    def generate_noiseModel(GateSet1, GateSet2):
        prob_1 = 0.001
        prob_2 = 0.01

        noiseM = noise.NoiseModel()
        noiseM.add_all_qubit_quantum_error(noise.depolarizing_error(prob_1,1), GateSet1)
        noiseM.add_all_qubit_quantum_error(noise.depolarizing_error(prob_2,2), GateSet2)

        return noiseM

    
    @staticmethod
    def QiskitTranspile(circ):
        cplGraphEdges = nx.complete_graph(circ.num_qubits + 1).edges()
        cplMap = [list(elem) for elem in cplGraphEdges]

        bend = AerSimulator(coupling_map=cplMap, basis_gates = ["id", "sx", "x","rz", "cx"], 
                            noise_model=ExperimentClass.generate_noiseModel(["id", "sx", "x","rz"],["cx"]))

        opt_circ = transpile(circ, backend=bend, optimization_level=3)

        return opt_circ, opt_circ.depth()
        
    @staticmethod
    def genJSSModelInstance(nm, nj, difficulty):
        jobs_dict, machines_dict, rigs_matrix = ExperimentClass.genExample(nj, nm, difficulty)
        
        dur = []
        for j in jobs_dict:
            dur.append(j["Duration"])
        
        A=1
        B=2*np.max(dur)
        C1=4*max(dur)**2
        lagrange_list = [A,B,C1]
        JSS_instance = JSS_model_part(nj, nm, jobs_dict, machines_dict, rigs_matrix, lagrange_list, "min_makespan", "qubo", [],0)
        
        return JSS_instance
        
    @staticmethod
    def genCouplingGraphFromPolynomial(polynomial):
        """
            Generates a multi graph from a given quark polynomial.
            Variables = nodes, variables in the same monomial = edge between them
        """
        
        pC = polynomial.compact()
        
        
        G = nx.MultiGraph()
        
        G.add_nodes_from(pC.variables)
        
        
        edges = []
        for monomial in pC.keys():
            edges += (list(combinations(monomial, 2)))
        
        G.add_edges_from(edges)
        
        avg_deg = 0
        for (node, deg) in G.degree():
            avg_deg += deg
        avg_deg /= len(G.nodes())
        
        num_edges = 0
        num_multi_edges = 0
        for n1 in G.nodes():
            for n2 in G.nodes():
                num_ed = G.number_of_edges(n1,n2)
                if num_ed >= 2:
                    num_multi_edges += num_ed
                num_edges += num_ed
        
        
        num_edges /= 2
        num_multi_edges /=2
        
        return G, G.degree(), avg_deg, num_edges, num_multi_edges
        
        
        
        
