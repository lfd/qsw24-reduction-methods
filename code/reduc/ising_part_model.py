import json
import numpy as np
import itertools
from pyqubo import Binary
from scipy.optimize import minimize
from qiskit import QuantumCircuit, Aer, transpile
from qiskit import QuantumCircuit
import re
from Transform_dict_arrays import *

import pprint

import warnings
warnings.filterwarnings('ignore')

np.set_printoptions(formatter={'float': '{: 0.4f}'.format})

class JSS_model_part:
    """
    A class for implementing Ising models and QUBOs for JSS problem with tsp encoding
    cases to be covered:

    - different models: minimize makespan or minimize total idle times
    
    """
    
    def __init__(self, num_jobs, num_machines, jobs_array, machine_array, rigs_array, lagrange_params, qubo_model, eval_method, optimal, opt_val, *shift_length):

        """
        Constructs all the necessary attributes for the QAOASolver object.

        Parameters
        ----------
            num_jobs : integer
            num_machines : integer
            jobs_array : dict
                dictionary with job as key and duration and gear group as entries
            machine_array : dict
                dictionary with machine number as key and initial gear group as entry
            rigs_array : np.array
                contains gear change times between gear groups. 
            lagrange_params: list of integers
                Lagrange parameters for objective terms and capacity. depends on model and encoding
            qubo_model: string
                which qubo objective model to use. Possible values are: 'min_makespan', 'min_idle'
            eval_method: string
                using qubo energy for evaluation or some other method. Possible values are: 'qubo', 'makespan'
            optimal: list of optimal solutions
                optimal solution, depends on qubo model. e.g. minimal makespan
            shift_length: integer
                length of the shift to be scheduled (optional). Only used for qubo_model = 'min_idle'.
        Returns
        ----------
            None

        """

        
        trans=Transform_Data_Jobscheduling()
        self.rigs, self.rigs_initial, self.durations = trans.rigs_to_Jobsarray(rigs_array, jobs_array, machine_array)

        
        self.num_jobs = num_jobs
        self.num_machines = num_machines
        self.optimal=optimal
        self.opt_val = opt_val
        self.qubo_model = qubo_model
        self.eval_method = eval_method

        self.A, self.B, self.C1 = lagrange_params
        
        if shift_length is not None:
            self.shift_length = shift_length
        elif shift_length is None:
            self.shift_length = sum(self.durations) + sum(sum(self.rigs)) + sum(sum(self.rigs_initial))

        opt_bitstrings=[]
        for j in range(len(optimal)):
            opt_bitstrings.append(self.schedule_to_bitstring(self.dict_to_schedule(optimal[j])))
        self.opt_bitstrings = opt_bitstrings
        
        
        (self.variables,
        self.linear,
        self.quadratic,
        self.third,
        self.fourth,
        self.var_to_index,
        self.index_to_var,
        self.model,
        self.offset) = self.create_ising_model_part1()
        self.num_qubits = len(self.variables)

        if eval_method == 'qubo':
            self.energy_penalty = self.qubo_energy_penalty_min_ms
 

    
    
    

    def create_ising_model_part1(self):

        """
        Parameters:
        num_jobs: Number of jobs (integer)  
        num_machines: Number of machines (integer)
        durations: Durations of the jobs in arbitrary time units (1-dim array of floats)
        rigs: Matrix denoting the differences in rigs between the jobs (array of integers with dimension num_jobs x num_jobs)
        rigs_initial: Matrix denoting the differences in rigs between the jobs and the initial rig of the machines (array of integers with dimension num_jobs x num_machines)

        Binary variables: x_ij = 1 if job i is processed on machine j. Solve only the distribution, the rest is going to be solved classically or with a second QUBO.
        We have N jobs and M machines. 
        Total number of variables: N*M
        
        Bitstrings: For each job, x_ij denote job i on machine j
        (x_01, x_02, ..., x_0M, x_10, ...,x_1M,..., x_NM)

        Lagrange multipliers for the different parts
        A: Minimize total duration (job durations + initial rigs. without rig changes for QUBO to have seond order)
        B: Minimize changing the rig between jobs as well as between the first job and the initial rig of the machine
        C1: Constraint that each job can run only on one machine
        """

        mach_pairs = list(itertools.combinations(range(self.num_machines), 2)) 
        job_pairs = list(itertools.combinations(range(self.num_jobs), 2)) 

        x = dict()
        for i in range(self.num_jobs):
            for j in range(self.num_machines):
                x[(i,j)] = Binary(f'x_{i}_{j}')

        
        H_obj=0
        for pair in mach_pairs:
            temp1 = sum([x[i,pair[0]] * (self.durations[i]+self.rigs_initial[i,pair[0]]) for i in range(self.num_jobs)]) + sum([x[job_pair[0],pair[0]]*x[job_pair[1],pair[0]]*self.rigs[job_pair[0],job_pair[1]] for job_pair in job_pairs])
            temp2 = sum([x[i,pair[1]] * (self.durations[i]+self.rigs_initial[i,pair[1]]) for i in range(self.num_jobs)]) + sum([x[job_pair[0],pair[1]]*x[job_pair[1],pair[1]]*self.rigs[job_pair[0],job_pair[1]] for job_pair in job_pairs])
            H_obj += (temp1-temp2)**2

        
        H_r=0
        for j in range(self.num_machines):
            H_r += sum([x[pair[0],j] * x[pair[1],j] * self.rigs[pair[0], pair[1]] for pair in job_pairs])
            H_r += sum([x[i,j] * self.rigs_initial[i, j] for i in range(self.num_jobs)])

        
        H_C1=0
        for i in range(self.num_jobs):
            H_C1 += (sum([x[i,j] for j in range(self.num_machines)]) - 1)**2

        H = self.A*H_obj + self.B*H_r + self.C1*H_C1

        model = H.compile()

        variables=[] 
        for i in range(self.num_jobs):
            for j in range(self.num_machines):
                variables.append('x_'+str(i)+'_'+str(j))

        linear, quadratic, offset = model.to_ising()

        var_to_index = dict([(n, i) for i, n in enumerate(variables)])
        index_to_var = dict([(i, n) for i, n in enumerate(variables)])

        lin_new = {}
        quadr_new = {}
        third_new = {}
        fourth_new = {}

        for var_dict in [linear, quadratic]:
            for key in var_dict.keys():
                
                nums=[int(s) for s in re.findall(r'\d+', str(key))]
                new_key=[]
                for var_count in range(int(len(nums)/2)):
                    
                    if f'x_{nums[2*var_count]}_{nums[2*var_count+1]}' not in new_key:
                        new_key.append(f'x_{nums[2*var_count]}_{nums[2*var_count+1]}')
                if len(new_key) == 1:
                    found = False
                    if tuple(new_key) in lin_new.keys():
                        lin_new[new_key[0]] += var_dict[key]
                        found = True
                    if not found:
                        lin_new[new_key[0]]=var_dict[key]
                if len(new_key) == 2:
                    found = False
                    for perm in list(itertools.permutations(tuple(new_key))):
                        
                        if perm in quadr_new.keys():
                            quadr_new[tuple(perm)]+= var_dict[key]
                            found = True
                    
                    if not found:
                        quadr_new[tuple(new_key)]=var_dict[key]
                if len(new_key) == 3:
                    found = False
                    for perm in list(itertools.permutations(tuple(new_key))):
                        
                        if perm in third_new.keys():
                            third_new[tuple(perm)]+= var_dict[key]
                            found = True
                    
                    if not found:
                        third_new[tuple(new_key)] = var_dict[key]
                if len(new_key) == 4:
                    found = False
                    for perm in list(itertools.permutations(tuple(new_key))):
                        
                        if perm in fourth_new.keys():
                            fourth_new[tuple(perm)]+= var_dict[key]
                            found = True
                    
                    if not found:
                        fourth_new[tuple(new_key)] = var_dict[key]
   
        return variables, lin_new, quadr_new, third_new, fourth_new, var_to_index, index_to_var, model, offset

    
    
    

    def create_quark_dict(self):

        quark_dict = {}
        for var_dict in [self.linear, self.quadratic]:
            for key in var_dict.keys():
                
                nums=[int(s) for s in re.findall(r'\d+', str(key))]
                new_key=[]
                
                for var_count in range(int(len(nums)/2)):
                    var = f'x_{nums[2*var_count]}_{nums[2*var_count+1]}'
                    index = self.var_to_index[var]
                    new_key.append(('x', index))
                    val = var_dict[key]
                found = False
                for perm in list(itertools.permutations(tuple(new_key))):
                    
                    if perm in quark_dict.keys():
                        quark_dict[tuple(perm)]+= val
                        found = True
                
                if not found:
                    quark_dict[tuple(new_key)]=val




    
    
    
    def qubo_energy_penalty_min_ms(self, x_val):
        """
        Method to compute the QUBO penalty and QUBO cost

        Parameters
        ----------
            x_val : bitstring

        Returns
        ----------
             cost : total QUBO energy including objective and penalty terms
             penalty : qubo penalty value
        """

        mach_pairs = list(itertools.combinations(range(self.num_machines), 2)) 
        job_pairs = list(itertools.combinations(range(self.num_jobs), 2)) 

        x = dict()
        for k in range(len(x_val)):
            variable = self.index_to_var[k]
            (i,j)=int(variable.split('_')[1]), int(variable.split('_')[2])
            x[(i,j)]=int(x_val[k])

        
        energy_obj=0
        for pair in mach_pairs:
            temp1 = sum([x[i,pair[0]] * (self.durations[i]+self.rigs_initial[i,pair[0]]) for i in range(self.num_jobs)]) + sum([x[job_pair[0],pair[0]]*x[job_pair[1],pair[0]]*self.rigs[job_pair[0],job_pair[1]] for job_pair in job_pairs])
            temp2 = sum([x[i,pair[1]] * (self.durations[i]+self.rigs_initial[i,pair[1]]) for i in range(self.num_jobs)]) + sum([x[job_pair[0],pair[1]]*x[job_pair[1],pair[1]]*self.rigs[job_pair[0],job_pair[1]] for job_pair in job_pairs])
            energy_obj += (temp1-temp2)**2
        energy_obj *= self.A

        
        energy_r=0
        for j in range(self.num_machines):
            energy_r += sum([x[pair[0],j] * x[pair[1],j] * self.rigs[pair[0], pair[1]] for pair in job_pairs])
            energy_r += sum([x[i,j] * self.rigs_initial[i, j] for i in range(self.num_jobs)])
        energy_r *= self.B
        
        
        energy_C1=0
        for i in range(self.num_jobs):
            energy_C1 += (sum([x[i,j] for j in range(self.num_machines)]) - 1)**2
        energy_C1 *= self.C1


        penalty = energy_C1
        cost = penalty + energy_obj + energy_r
        return cost, penalty

    

    
    
    
   
   
    def dict_to_schedule(self, dict):
        
        schedule=np.zeros((self.num_machines, self.num_jobs))
        for i in range(self.num_machines):
            for l in range(len(dict[str(i)])):
                    schedule[i,l]=dict[str(i)][l]+1
        return schedule

    def schedule_to_bitstring(self, schedule):
        
        bitstring_new=np.zeros((self.num_machines*self.num_jobs))

        for j in range(self.num_machines):
            for i in range(self.num_jobs):
                if schedule[j,i] > 0:
                    bitstring_new[int(schedule[j,i]-1)*self.num_machines+j] = 1

        return bitstring_new


    def bitstring_to_schedule(self, bitstring):

        schedule=np.zeros((self.num_machines, self.num_jobs))
        for j in range(self.num_machines):
            count=0
            for i in range(self.num_jobs):
                if bitstring[i*self.num_machines+j] > 0:
                    schedule[j,count] = i+1
                    count+=1
        return schedule


    def opt_dict_to_bitstring_schedule(self, opt_sols):
        """
        get simplified schedule from optimal solutions dict: All jobs are just sorted in ascending order. 
        get bitstrings in part1-encoding corresponding to the optimal solutions
        To be compared with solutions for the part 1 problems, which include no ordering
        """
        
        opt_schedules=[]
        opt_bitstrings=[]
        for opt_sol in opt_sols:

            schedule=np.zeros((self.num_machines, self.num_jobs))
            for i in range(self.num_machines):
                job_list= np.sort(opt_sol[str(i)])
                for l in range(len(job_list)):
                    schedule[i,l]=job_list[l]+1
            opt_schedules.append(schedule)

        for schedule in opt_schedules:
            bitstring_new=np.zeros((self.num_machines*self.num_jobs))

            for j in range(self.num_machines):
                for i in range(self.num_jobs):
                    if schedule[j,i] > 0:
                        bitstring_new[int(schedule[j,i]-1)*self.num_machines+j] = 1
            opt_bitstrings.append(bitstring_new)

        return opt_bitstrings, opt_schedules
    

    def check_optimal_bitstring(self, bitstring):
        opt_bitstrings, _ = self.opt_dict_to_bitstring_schedule(self.optimal)
        opt = False
        for opt_bitstring in opt_bitstrings:
            if np.sum(np.abs(np.array(opt_bitstring) - np.array(bitstring))) == 0:
                opt= True
            return opt