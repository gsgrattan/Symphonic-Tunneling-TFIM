# Symphonic_Tunneling-TFIM
<h3>Enviornent: </h3>
The enviornment and packages in stenv were installed and running on an Apple M2 chip,
the venv may need to be reconfigured for other systems.

src/helpers.py contains various functions for creating the data inputs to enabling parallel simulation and generating filenames for the data files

src/tfim_st.py contains the various functions that perform the simulation

there are two main functions to note:
    TFIM_evolve_2d_sampling(**param_dict)
        Samples <ZZ> only at the end of the plateau
    TFIM_evolve_2d_continous_sampling(**param_dict)
        Sample <ZZ> at every trotter step in the plateau

main.ipynb:
    The output datafiles are pickled pandas dataframe objects and should be stored in a directory in Symphonic_Tunneling-TFIM/data/

