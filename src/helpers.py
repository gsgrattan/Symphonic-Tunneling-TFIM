import numpy as np
import copy
import sys

qubit_pairs = { #which spins to chose for the ZZ correlator
    "(2,5)" : (0,7),
    "(2,6)": (0,9),
    "(2,7)" :(0,10),
    "(2,8)" : (0,12),
    "(2,9)" : (0,13),
    "(2,10)" : (0,15),
    "(2,11)" : (0,16),
    "(2,12)" : (0,18),
    "(2,13)" : (0,19),
    "(2,14)" : (0,21),
    "(2,15)" : (0,22),

    "(3,4)": (0,7),
    "(3,5)": (0,7),
    "(3,6)": (0,10),
    "(3,7)": (0,10),
    "(3,8)": (0,13),
    "(3,9)": (0,13),
    "(3,10)": (0,16),

    "(4,4)" : (0,10),
    "(4,5)" : (0,10),
    "(4,6)" : (0,14),
    "(4,7)" : (0,14),

    "(5,5)" : (0,12),
    "(5,6)" : (0,17)}

parameter_dict_skeleton = { #This is just a skeleton for what we pass into each simulation
    'lx':0,
    'ly':0,
    'p0':0,
    'p1':0,
    'kappa':0,
    'alpha' : 0,
    'alpha_coef': 0,
    'freq': 0,
    'freq_coef': 0,
    'ramptime': 0,
    'plateau_time' : 0,
    'dt' : 0,
    'num_samples': 1000
}


#Creates an array of possible parameter combinations to feed into the Parallel package, should be easy to generalize to arbitrary parameter choices
def generate_params(steps,kappa_vals, lattice_arr, alpha_coef_vals, freq_coef, dt_coef, ramptime_divisor, sim_time_coef, num_samples):
    
    param_array = []
    for lx,ly in lattice_arr:
        for alpha_coef in alpha_coef_vals:
            for kappa in kappa_vals:
                for step in range(steps):
                    param = copy.deepcopy(parameter_dict_skeleton)
                
                
     
                    #skeleton for 
                    param['lx'] = lx
                    param['ly'] = ly
                    n = lx*ly
                    
                    param['num_qubits'] = n
                
                    param['kappa'] = kappa
                    frequency = freq_coef*np.log(n)
                    param['freq'] = frequency
                    param['freq_coef'] = freq_coef

                    param['alpha'] = alpha_coef * np.log(n)/np.log(8) #Log(8) is legacy

                    param['alpha_coef'] = alpha_coef

                    param['dt'] = dt_coef/frequency
                    param['plateau_time'] = sim_time_coef*(0.8/frequency + (n*0.25 - 0.8/frequency)*step/(steps-1))
                    ramptime = n/ramptime_divisor

                    param['ramptime'] = ramptime

                    param['num_samples'] = num_samples
                    
                    if lx == 1:                    
                        p0 = 0
                        p1 = n//2
                    elif (lx in {2, 3 ,4 ,5}): 
                        lattice = "({},{})".format(lx,ly)
                        p0,p1 = qubit_pairs.get(lattice)

                    param['p0'] = n-1-p0 #reindex for qulacs
                    param['p1'] = n-1-p1
            
                    param_array.append(param)

    return param_array


def generate_filename(filename_prefix, param):
    lx,ly = (param['lx'], param['ly'])
    alpha_coef = str(int(param['alpha_coef']*100))
    freq_coef = param['freq_coef']
    kappa = str(int(param['kappa']*100))
    
    base_filename = filename_prefix + '_d{}x{}_k{}_a{}_f{}.pkl' 
    return base_filename.format(lx,ly,kappa,alpha_coef,freq_coef)
