import qulacs
from qulacs import *
from qulacs.gate import *
import numpy as np
import math
import random
import pickle
import copy
import networkx as nx

from joblib import Parallel, delayed
import multiprocessing

#Plotting and Data Manipulation
import matplotlib.pyplot as plt
import pandas as pd

import time

# uniform X layer -- note it applies 2*ang to be consistent with the conventions in my earlier cirq code, feel free to change as needed

def uniform_TF_layer(n,ang):
    circuit = QuantumCircuit(n)
    for i in range(n):
        circuit.add_RX_gate(i,2*ang) # argument order is reversed compared to cirq
    return circuit

# uniform AC Y layer

def uniform_AC_Y_layer(n,T,freq,alpha,ang):
    circuit = QuantumCircuit(n)
    amp = 2*ang*alpha*np.sin(2*np.pi*freq*T)
    for i in range(n):
        circuit.add_RY_gate(i,amp) # argument order is reversed compared to cirq
    return circuit

# 1d ZZ term

def ZZ_1d_PBC(n,ang):
    circuit = QuantumCircuit(n)
    for i in range(n-1):
        circuit.add_gate(PauliRotation([i,i+1],[3,3],2*ang))
    circuit.add_gate(PauliRotation([0,n-1],[3,3],2*ang)) # enforce PBCs
    return circuit

# measure ZZ from a vector (note reversed qubit order in the output of qulacs when choosing p1, p2)
# TODO: implement this same thing with samples instead 

def measure_ZZ_statevec(n,p1,p2,statevec):
    res = 0
    binaryformat = '{0:0' + str(n) + 'b}'
    for i in range(len(statevec)):
        bitstring=binaryformat.format(i)
        res = res + (-1)**(int(bitstring[p1])+int(bitstring[p2]))*abs(statevec[i])**2
    return res


def measure_ZZ_samples(n,p1,p2,samples):
    res = 0
    binaryformat = '{0:0' + str(n) + 'b}'
    for i in range(len(samples)):
        bitstring=binaryformat.format(samples[i])
        res += (-1)**(int(bitstring[p1])+int(bitstring[p2]))
    return (res/len(samples))

def measure_PFM_samples(n,samples):
    return [list(samples).count(0)/len(samples),list(samples).count(2**n-1)/len(samples)]
    



# full TFIM circuit, assuming J = 1
# this is done in terms of time instead of layers

def TFIM_evolve_1d_with_measurement(n,kappa,freq,alpha,ramptime,plateau_time,dt): # 
    res = []
    ZZvals = []
    psi = QuantumState(n)
    psi.set_zero_state() # initializes in zero
    Tseg = 0
    T = 0
    baseval = - 2*np.pi*dt # note - sign absorbed here
    # ramp up evolution
    while(Tseg<ramptime): 
        Tseg = Tseg + dt
        T = T + dt
        ramping_value = np.sin(np.pi*Tseg/(2*ramptime)) # ramping profile
        exponent = ramping_value * kappa * baseval # factor of 2 was absorbed into my uniform TF layer
        uniform_TF_layer(n,exponent).update_quantum_state(psi) # apply transverse field
        uniform_AC_Y_layer(n,T,freq,alpha,baseval * ramping_value).update_quantum_state(psi) # VHF AC terms (note: ramped)
        ZZ_1d_PBC(n,baseval).update_quantum_state(psi) # Ising Hamiltonian
        
    # now main hold evolution
    
    Tseg = 0
    while(Tseg<plateau_time):
        Tseg = Tseg + dt
        T = T + dt
        ramping_value = 1 # hold time
        exponent = ramping_value * kappa * baseval # factor of 2 was absorbed into my uniform TF layer
        uniform_TF_layer(n,exponent).update_quantum_state(psi) # apply transverse field
        uniform_AC_Y_layer(n,T,freq,alpha,baseval * ramping_value).update_quantum_state(psi) # VHF AC terms
        ZZ_1d_PBC(n,baseval).update_quantum_state(psi) # Ising Hamiltonian
        # do the ZZ measurement at each timestep and store it as a list
        # this is hard coded here to measure 0, L/2 but can be easily tweaked
        ZZvals.append((Tseg,measure_ZZ_statevec(n,0,math.floor(n/2),psi.get_vector())))
    
    Tseg = 0
    while(Tseg<ramptime): 
        Tseg = Tseg + dt
        T = T + dt
        ramping_value = np.cos(np.pi*Tseg/(2*ramptime)) # ramping profile
        exponent = ramping_value * kappa * baseval # factor of 2 was absorbed into my uniform TF layer
        uniform_TF_layer(n,exponent).update_quantum_state(psi) # apply transverse field
        uniform_AC_Y_layer(n,T,freq,alpha,baseval * ramping_value).update_quantum_state(psi) # VHF AC terms (note: ramped)
        ZZ_1d_PBC(n,baseval).update_quantum_state(psi) # Ising Hamiltonian
        
    sv = psi.get_vector() # get the state vector back
    # return P0, P1, ZZ correlator
    return [abs(sv[0])**2,abs(sv[2**n - 1])**2,ZZvals]

def ZZ_graph(n,problem_graph,ang):
    circuit=QuantumCircuit(n)
    for edge_index, edge in enumerate(problem_graph):
            i, j = edge
            circuit.add_gate(PauliRotation([i,j],[3,3],2*ang))
    return circuit



def TFIM_evolve_2d_sampling(**param_dict): 
    lx = param_dict['lx']
    ly = param_dict['ly']
    p0 = param_dict['p0']
    p1 = param_dict['p1']
    kappa = param_dict['kappa']
    freq = param_dict['freq']
    alpha = param_dict['alpha']
    ramptime = param_dict['ramptime']
    plateau_time = param_dict['plateau_time']
    dt = param_dict['dt']
    num_samples = param_dict['num_samples']
    
    res = []
    ZZvals = dict()
    
    n = lx * ly
    
    graph = nx.grid_2d_graph(n=lx, m=ly, periodic=True)
        # Convert indices from adj matrix form A[i,j] --> single index V[i]
    graph = nx.convert_node_labels_to_integers(graph)
    problem_graph = graph.edges
    
    psi = QuantumState(n)
    psi.set_zero_state() # initializes in zero
    

    baseval = - 2*np.pi*dt # note - sign absorbed here

    num_ramping_layers = ramptime // dt
    
    # ramp up transverse field and Symphonic drive
    Tseg = 0
    trotter_index = 0  #Trotter layer index
    T = 0
    while(trotter_index < num_ramping_layers): 
        trotter_index+=1
        T = T + dt
        Tseg += dt
        ramping_value = np.sin(np.pi*Tseg/(2*ramptime)) # ramping profile
        exponent = ramping_value * kappa * baseval # factor of 2 was absorbed into my uniform TF layer
        ZZ_graph(n,problem_graph,baseval).update_quantum_state(psi) # Ising Hamiltonian

        uniform_TF_layer(n,exponent).update_quantum_state(psi) # apply transverse field
        uniform_AC_Y_layer(n,T,freq,alpha,baseval * ramping_value).update_quantum_state(psi) # VHF AC terms (note: ramped)
                
    # now main hold evolution
    
    Tseg = 0

    num_evolution_layers = plateau_time//dt

    #get the timesteps to sample ZZ
    sample_trotter_steps = np.linspace(start = num_ramping_layers, stop = num_ramping_layers+num_evolution_layers, num=40, dtype = int) #set the sample timesteps

    while(trotter_index< (num_evolution_layers+ num_ramping_layers)):
        trotter_index+=1
        T = T + dt
        Tseg += dt
        ramping_value = 1 # hold time
        exponent = ramping_value * kappa * baseval # factor of 2 was absorbed into my uniform TF layer
        ZZ_graph(n,problem_graph,baseval).update_quantum_state(psi) # Ising Hamiltonian
        uniform_TF_layer(n,exponent).update_quantum_state(psi) # apply transverse field
        uniform_AC_Y_layer(n,T,freq,alpha,baseval * ramping_value).update_quantum_state(psi) # VHF AC terms
        
        # do the ZZ measurement at each timestep and store it as a list
        # this is hard coded here to measure 0, L/2 but can be easily tweaked
        if trotter_index in sample_trotter_steps:

            ZZvals[Tseg] = psi.sampling(num_samples)
       
    #Ramp Down transverse Field and Symphonic Drive
    Tseg = 0
    while(trotter_index<(num_evolution_layers + 2*num_ramping_layers)):
        trotter_index+=1 
        Tseg += dt
        T = T + dt
        ramping_value = np.cos(np.pi*Tseg/(2*ramptime)) # ramping profile
        exponent = ramping_value * kappa * baseval # factor of 2 was absorbed into my uniform TF layer
        ZZ_graph(n,problem_graph,baseval).update_quantum_state(psi) # Ising Hamiltonian

        uniform_TF_layer(n,exponent).update_quantum_state(psi) # apply transverse field
        uniform_AC_Y_layer(n,T,freq,alpha,baseval * ramping_value).update_quantum_state(psi) # VHF AC terms (note: ramped)
                
    # return P0, P1, ZZ correlator
    ZZ_mean_array = [measure_ZZ_samples(n, p0, p1, samples) for time, samples in ZZvals.items()] 

    ZZ_mean = np.mean(ZZ_mean_array)


    PFM_vals = measure_PFM_samples(n,psi.sampling(num_samples))
    data = {"zeros_overlap" : PFM_vals[0],
            "ones_overlap" : PFM_vals[1],
            "ZZ_mean": ZZ_mean,
            "ZZ_sample_dict" : [ZZvals]
            }
    
    return {**param_dict, **data}


def TFIM_evolve_2d_continous_sampling(**param_dict): 
    lx = param_dict['lx']
    ly = param_dict['ly']
    p0 = param_dict['p0']
    p1 = param_dict['p1']
    kappa = param_dict['kappa']
    freq = param_dict['freq']
    alpha = param_dict['alpha']
    ramptime = param_dict['ramptime']
    plateau_time = param_dict['plateau_time']
    dt = param_dict['dt']
    num_samples = param_dict['num_samples']
   

    ZZvals = dict()
    
    n = lx * ly
    
    graph = nx.grid_2d_graph(n=lx, m=ly, periodic=True)
        # Convert indices from adj matrix form A[i,j] --> single index V[i]
    graph = nx.convert_node_labels_to_integers(graph)
    problem_graph = graph.edges

    
    psi = QuantumState(n)
    psi.set_zero_state() # initializes in zero
    



    baseval = -2*np.pi*dt # note - sign absorbed here


    num_ramping_layers = ramptime // dt
    num_evolution_layers = plateau_time //dt

    # ramp up transverse field and Symphonic drive
    Tseg = 0
    trotter_index = 0  #Trotter layer index
    T = 0
    while(trotter_index < num_ramping_layers): 
        trotter_index+=1
        T = T + dt
        Tseg += dt
        ramping_value = np.sin(np.pi*Tseg/(2*ramptime)) # ramping profile

        exponent = ramping_value * kappa * baseval # factor of 2 was absorbed into my uniform TF layer
        ZZ_graph(n,problem_graph,baseval).update_quantum_state(psi) # Ising Hamiltonian

        uniform_TF_layer(n,exponent).update_quantum_state(psi) # apply transverse field
        uniform_AC_Y_layer(n,T,freq,alpha,baseval * ramping_value).update_quantum_state(psi) # VHF AC terms (note: ramped)
                



    # Plateau
    Tseg = 0
    while(trotter_index < (num_evolution_layers + num_ramping_layers)):
        trotter_index+=1
        T = T + dt
        Tseg += dt

        ramping_value = 1. # hold time
        exponent = ramping_value * kappa * baseval # factor of 2 was absorbed into my uniform TF layer

        ZZ_graph(n,problem_graph,baseval).update_quantum_state(psi) # Ising Hamiltonian
        uniform_TF_layer(n,exponent).update_quantum_state(psi) # apply transverse field
        uniform_AC_Y_layer(n,T,freq,alpha,baseval * ramping_value).update_quantum_state(psi) # VHF AC terms 
        ZZvals[str(T)] = psi.sampling(num_samples)
        # do the ZZ measurement at each timestep and store it as a list
        # this is hard coded here to measure 0, L/2 but can be easily tweaked
       
    #Ramp Down transverse Field and Symphonic Drive
    Tseg = 0

    while(trotter_index < (num_evolution_layers + 2*num_ramping_layers)):
        trotter_index+=1 
        Tseg += dt
        T = T + dt
        ramping_value = np.cos(np.pi*Tseg/(2*ramptime)) # ramping profile
        exponent = ramping_value * kappa * baseval # factor of 2 was absorbed into my uniform TF layer
        ZZ_graph(n,problem_graph,baseval).update_quantum_state(psi) # Ising Hamiltonian

        uniform_TF_layer(n,exponent).update_quantum_state(psi) # apply transverse field
        uniform_AC_Y_layer(n,T,freq,alpha,baseval * ramping_value).update_quantum_state(psi) # VHF AC terms (note: ramped)
                
    # return P0, P1, ZZ correlator
    ZZ_mean_array = [measure_ZZ_samples(n, p0, p1, samples) for time, samples in ZZvals.items()] 

    ZZ_mean = np.mean(ZZ_mean_array)


    PFM_vals = measure_PFM_samples(n,psi.sampling(num_samples))

    data = {"zeros_overlap" : PFM_vals[0],
            "ones_overlap" : PFM_vals[1],
            "ZZ_mean" : ZZ_mean,
            "ZZ_sample_dict" : [ZZvals]
            }
    return {**param_dict, **data}