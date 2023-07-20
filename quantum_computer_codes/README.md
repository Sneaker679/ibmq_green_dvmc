# Overview of the parameters in `parameters.py`.

## IBM Credentials
- run_on_quantum_computer (Y/N) : If you want to run the code on a quantum computer, change this to 'Y'. Make sure your IBM credentials are correctly configured below.
- channel : Channel to connect to. To connect to IBM quantum computers, input 'ibm_quantum'
- token : Refer to your IBM account for your token.
- backend_device : Name of the device to run the simulation on. The list of available backends is on the IBM website. The link is provided [here](https://quantum-computing.ibm.com/services/resources?tab=yours).


## Physic parameters
- N : Number of sites in the lattice.
- t : Value of the hopping energy. Usually equal to -1.
- U : Value of the potential energy. Usually equal to 2.
- mu : Value of the chemical potential energy. Usually equal to half of U.
- spin_green : The spin to be used to calculate the Green function. This dictates the spins of the operators when calculating the elements of the H+,H-,S+ and S- matrices.
- spin_gs : This parameter is for telling the code which ground state should be use in preference. This is important, because some ground states may have the same ground state energy, and thus the code doesn't know which ground state to use to perform the calculations. With this parameter, you can tell the code to have a preference for lower total spin ground states of higher total spin ground states. With this preference in mind, the code will always try to choose the ground state with the total spin closest to 0.  


## Code parameters
- use_qcm (Y/N) : Calls or doesn't call the QCM module in the code. This is because the QCM installation may be difficult for some, thus this benchmark is optionnal.

- force_custom_lattice (Y/N): This will tell the code to use the custom lattice models inputted in `parameters.py` instead of the automatically generated lattice.
- hopping_matrix_for_qiskit_lattice (Y/N) (NO LONGER A PARAMETER SINCE A SIMPLER SOLUTION WAS FOUND): If force_custom_lattice == 'Y', the current implementation can initialize a qiskit lattice in 2 ways. The first one is by using a hopping matrix, defined in `parameters.py`. The second option is to initialize the lattice using qiskit's Lattice model classes and work from there. Highly custom lattices can be created from it, which could be very useful for the user. However, by using this second method, it is possible that the spectrums of the benchmarks don't overlap with the spectrum of qiskit, because the labels of the sites are different depending on the source code of the packages. This is why, as the user, you get to choose between using a hopping matrix to initialize the lattice, or using the lattice classes of qiskit. 

- force_custom_circuit (Y/N): By default, the code will initialize the qubits of the circuit with an exact ground state vector. If you wish to test your own custom circuit, input 'Y'.
- decompose_and_print_circuit (Y/N): The circuit used in the code can be decomposed by qiskits into more elementary gates. By decomposing the circuit, said circuit will be printed in the CMD for you to see, and this decomposed circuit will also be used in the code. You will quickly notice that it is not very practical to have this on, because calculation times are way longer and the circuit itself becomes enourmous.

- generate_npy (Y/N): This is for saving, or not, the matrices as `.npy` files. These should always be left on 'Y', because `graph.py` uses these files as input for the program.
- generate_matrix (H+,H-,S+,S-,ALL) : By default, this is on 'ALL'. It means that everytime you run the code, all the matrices are going to be calculated again. If you only want to calculate one, then you change this parameter to reflect that choice.    
- excit_document : This is the name of the `excitation.def` file to be used by the program. Said file needs to be in the same directory as `parameters.py`.


## Noisy simulation parameters
- noisy_simulation (Y/N): Noise can be added to the simulation. Input 'Y' to have noise. Note that qiskit does not parallelize this simulation, thus the calculation can take days to complete. Not very practical. 


# Custom Circuits
With this code, you can use your custom circuits to calculate the Green function. However, be aware that qiskit has a very specific notation to use. Indeed, each qubits correspond to these positions in fock's basis: |... 5 4 3 2 1 0 >. In other words, the first qubit starts to the right. This is important because in order to create, for example, state #8, you need to apply circuit.x(3) and NOT circuit.x(0). 

