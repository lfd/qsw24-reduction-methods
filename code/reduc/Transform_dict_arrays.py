import numpy as np

class Transform_Data_Jobscheduling():
    def __init__(self) -> None:
        pass  
    def rigs_to_Jobsarray(self, rigs_array, jobs_array, machine_array):
        """
        Transforms dictionary with rigs to rigs array:duration of gearchange to 
        an array with jobs in the index and the durations in the matrix.

        Transforms jobsdictionary with job as key and duration and gear group as entries
        to durations array.

        Transforms machinedict with m as key and initial gear group as entry to
        initial_rigs array with machines as lines and jobs as columns.

        Output: rigs, rigs_initial, durations
        """
        num_jobs = len(jobs_array)
        num_machines = len(machine_array)

        durations = np.zeros(num_jobs).astype(int)
        rigs_initial = np.zeros((num_jobs, num_machines))
        rigs = np.zeros((num_jobs, num_jobs))
        

        
        for i in range(num_jobs):
            durations[i] = jobs_array[i]["Duration"]

        
        for j in range(num_machines):
            rj = machine_array[j]["Rig"]
            for i in range(num_jobs):
                ri = jobs_array[i]["Rig"]
                rigs_initial[i,j] = rigs_array[rj][ri]


        for i1 in range(num_jobs):
            r1 = jobs_array[i1]["Rig"]
            for i2 in range(num_jobs):
                r2 = jobs_array[i2]["Rig"]
                rigs[i1,i2] = rigs_array[r1][r2]

        return rigs, rigs_initial, durations

