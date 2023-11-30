# Symphonic_Tunneling-TFIM
<h3>Enviornment and dependencies:</h3>
<p>
    See requirements.txt for the dependencies
</p>

<h3>main.ipynb: </h3>
<p>
    - Uses the joblib.Parallel() function to parallelize the simulations over each trajectory.
       NOTE: It is programmed all cores on the machine by default, can be changed by altering `num_parallel_jobs`

    - The output datafiles are pickled pandas dataframe objects and should be stored in a directory with the path `Symphonic_Tunneling-TFIM/data/`

</p>
<h3> `src/tfim_st.py` </h3>

<p>
contains the various functions that perform the simulation

there are two main functions to note:
    TFIM_evolve_2d_sampling(**param_dict)
        Samples `<ZZ>` only at the end of the plateau
    TFIM_evolve_2d_continous_sampling(**param_dict)
        Sample `<ZZ>` at every trotter step in the plateau
</p>

<h3>`src/helpers.py` </h3>
<p>
contains various functions for creating the data inputs to enabling parallel simulation and generating filenames for the data files
</p>



