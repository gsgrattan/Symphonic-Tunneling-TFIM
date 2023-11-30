# Symphonic_Tunneling-TFIM
<h3>Enviornment and dependencies:</h3>


Here is a brief overview of the enviornment that last run this notebook
<p>

    Python version: 3.9.6 (default, May  7 2023, 23:32:44) 
    [Clang 14.0.3 (clang-1403.0.22.14.1)]
    numpy==1.26.2
    pandas==2.1.3
    scipy==1.11.4
    qulacs==0.6.2
    networkx==3.2.1

</p>
See requirements.txt for a full list of dependencies

<h3>main.ipynb: </h3>

This contains driver code for the simulations of varying lattice x dimensions for the Transverse Field Ising Model. It  creates a quick rabi oscillation plot, and 
<p>

    - Uses the joblib.Parallel() function to parallelize the simulations over each trajectory.
       NOTE: It is programmed all cores on the machine by default, can be changed by altering 
       `num_parallel_jobs`

    - By default we run 48 different plateau times to collect the Rabi Oscillation Data, this is 
    arbitrary to the machine used to collect the data and the number of cores it had. this can be
     changed but it should remain relatively large to limit the impacts of aliasing.

    - I have commented out the section that saves the data, comment it back if you want

    - When turned on, the output datafiles are pickled pandas dataframe objects and are stored by 
     default in a directory with the path `Symphonic_Tunneling-TFIM/data/`, an empty directory `data` 
     has been provided for convenience
</p>

<h3> `src/tfim_st.py` </h3>
Contains the various functions that construct the circuit and run the simulations

There are two main functions to note:
<p>

    - TFIM_evolve_2d_sampling(**param_dict)
        Samples `<ZZ>` at 40 evenly spaced times during the plateau
    - TFIM_evolve_2d_continous_sampling(**param_dict)
        Sample `<ZZ>` at every trotter step in the plateau
</p>

<h3>`src/helpers.py` </h3>
Contains various functions for creating the data inputs to enabling parallel simulation and 
generating filenames for the data files
