# Overview of the parameters in `parameters.py`.

## IBM Credentials
- run_on_quantum_computer : Boolean. If you want to run the code on a quantum computer, change this to 'True'. Make sure your IBM credentials are correctly configured below.
- max_circuit_per_job : Integer that states how many circuits should a single job have. The code splits the workload on multiple jobs, hence why this parameter may be important.
- channel : Channel to connect to. To connect to IBM quantum computers, input 'ibm_quantum'
- token : Refer to your IBM account for your token.
- backend_device : Name of the device to run the simulation on. The list of available backends is on the IBM website. The link is provided [here](https://quantum-computing.ibm.com/services/resources?tab=yours).
- custom_qubits : Inside this list, specify the qubits to be used on the quantum computer. If left empty, it will try to find the best qubits, although the code might make mistakes and also not work if the circuit is insufficiently big.
- optimisation_level : How much runtime is optimized. See qiskit's documentation for more info on this parameter.
- resilience_level : How much the execution resists noise. See qiskit's documentation for more info on this parameter.
max_execution_time : Maximum runtime of the job. Leave it to None to not apply any limits.
- execution (dict) : Dictionnary containing some important parameters, like the number of shots of the circuit. Shots are the number of time a single circuit is ran for a single measure. More shots mean less noise, but longer runtimes.
- recover_jobs: Boolean that dictates if the code should try to recover jobs that were already ran on IBM quantum.
- job_ids : Dictionnary that contains the jobs ids of the jobs to be recovered. Format : {'job0':'id1', 'job3':'id3',...}. Downside is that you have to know in what order the jobs were sent.


## Physic parameters
- N : Number of sites in the lattice.
- t : Value of the hopping energy. Usually equal to -1.
- U : Value of the potential energy. Usually equal to 2.
- mu : Value of the chemical potential energy. Usually equal to half of U.
- spin_green : The spin to be used to calculate the Green function. This dictates the spins of the operators when calculating the elements of the H+,H-,S+ and S- matrices.
- spin_gs : This parameter is for telling the code which ground state should be use in preference. This is important, because some ground states may have the same ground state energy, and thus the code doesn't know which ground state to use to perform the calculations. With this parameter, you can tell the code to have a preference for lower total spin ground states of higher total spin ground states. With this preference in mind, the code will always try to choose the ground state with the total spin closest to 0.  


## Code parameters
- use_qcm Boolean : Calls or doesn't call the QCM module in the code. This is because the QCM installation may be difficult for some, thus this benchmark is optionnal.

- force_custom_lattice : Boolean. This will tell the code to use the custom lattice models inputted in `parameters.py` instead of the automatically generated lattice. The code will always try to generate a lattice that forms a grid. If the code cannot make the number of sites inputted into a grid shapped lattice, the code will HAVE to use the custom lattice you inputted EVEN if this is set to False. The code cannot generate grid shapped lattices if the inputted number of sites is a prime number. In other words, if N == prime_number, then force_custom_lattice = True.
- force_custom_circuit : Boolean. By default, the code will initialize the qubits of the circuit with an exact ground state vector. If you wish to test your own custom circuit, input "True".
- decompose_and_print_circuit : Boolean that states if the circuit used in the code should be decomposed by qiskits into more elementary gates. By decomposing the circuit, said circuit will be printed in the CMD for you to see, and this decomposed circuit will also be used in the code. You will quickly notice that it is not very practical to have this on, because calculation times are way longer and the circuit itself becomes enourmous.
- produce_latex_circuit : Boolean that states if the latex commands to create an image of the used circuit should be printed in the cmd.

- generate_npy : Boolean. This is for saving, or not, the matrices as `.npy` files. This should always be left on 'True', because `graph.py` uses these files as input for the program.
- generate_matrix (H+,H-,S+,S-,ALL) : By default, this is on 'ALL'. It means that everytime you run the code, all the matrices are going to be calculated again. If you only want to calculate one, then you change this parameter to reflect that choice.    
- excit_document : This is the name of the `excitation.def` file to be used by the program. Said file needs to be in the same directory as `parameters.py`.

- parrallelize_observable_calculation : Boolean that states if the code should be parallelized for optimal run times.

## Noisy simulation parameters
- noisy_simulation : Boolean. Noise can be added to the simulation. Input 'True' to have noise. Note that qiskit does not parallelize this simulation, thus the calculation can take days to complete. Not very practical. 


# Graph parameters
graph_for_ibmq: Boolean that dictates if `combined_graph.py` should graph the ibmq spectrum.
graph_for_fock: Boolean that dictates if `combined_graph.py` should graph the fock benchmark spectrum.
graph_for_qcm: Boolean that dictates if `combined_graph.py` should graph the qcm spectrum.


# Custom Circuits
With this code, you can use your custom circuits to calculate the Green function. However, be aware that qiskit has a very specific notation to use. Indeed, each qubits correspond to these positions in fock's basis: |... 5 4 3 2 1 0 >. In other words, the first qubit starts to the right. This is important because in order to create, for example, state #8, you need to apply circuit.x(3) and NOT circuit.x(0). 
