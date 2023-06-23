# ibmq_green_dvmc
## Table of Contents
- [Installation](#Installation)
- [Description](#Description)
- [How it works](#how-it-works)
- [Execution](#execution)
- [Todo](#todo)

## Installation
This is a Python code. Use the latest version of language as of 15/06/2023 to avoid issues.
To run this code, you will need to install the following packages for your Python. The links to their documentation is also provided :
- qiskit-->  (https://qiskit.org/documentation/getting_started.html)
- qiskit_nature--> (https://qiskit.org/ecosystem/nature/getting_started.html)
- pyqcm-->  (https://qcm-wed.readthedocs.io/en/latest/intro.html#installation)
- matplotlib--> (https://matplotlib.org/stable/users/getting_started/index.html#installation-quick-start)
- numpy--> (https://numpy.org/install/)

Note that the "qcm" package is not manditory as its only purpose is to benchmark the calculations using the quantum computer.

Also note that in its current form, the "qcm" package needs a custom modification so it can work with the current code. In the "qcm" package, locate `pyqcm/_spectral.py` and find in it the "*cluster_spectral_function*" function. Remove the hashtags preceding the "*return*" at the end. You will need to rebuild the package and perform all the installation steps again, as described by the documentation.

## Description
This code is meant to calculate the Green Function for a specified system of low amount of sites. Instead of calculating it using conventional binary based computers, we simulate the calculation on a quantum computer using the qiskit package. 

The code itself is accompanied by a few supplementary benchmarks codes. These benchmark codes come in two types : benchmarks that use a custom second quantization base implementation, and benchmarks that use the "qcm" package. The usage of those benchmarks is explained in the following section.

## How it works
The primary use of this code is to calculate the Green Function. I will walktrough you in calculating it using the codes in `quantum_computer_codes/`.

All the files located in `quantum_computer_codes/` are useful to the calculation of the Green Function for any systems. Here is the uses of each individual files in this folder.

---
- The `parameters.py` file is where you can modify the parameters and variables to be used by all the files in the folder. It is also where you need to modify your lattice structures and circuit configurations if you are opting for a customized model. Specifications on the inner workings of the lattice and circuit objects are located in the "qiskit" package documentation. As for the qcm lattice, its usage is specified in its own documentation. One thing to note though it the HOPPING MATRIX. If the code detects it cannot create a grid lattice with the specified number of sites, you need to change manually the hopping matrix.
- The `generate.py` is meant to generate the specified matrices, print them in the CMD and save them as `.npy` files in `output/`.
- The `graph.py` is to make a graph with the generated matrices. The code will use the right set of matrices depending if you are generating a graph from the quantum computer matrices or the fock benchmark matrices. The spectrum is outputted as a pdf in the same folder as `parameter.py`, and its file name corresponds to which matrices were used (ibmq or fock). This code also outputs a `local_dos.dat` file in `output/`.
- The `hamiltonian_circuit.py` contains the hidden code to generate automatically grid shapped lattices for all implementations of the green function calculation (qiskit, qcm and fock).
---

- The `qcm_benchmark.py`, similarly to `graph.py`, will generate a graph of the Green Function using the "qcm" package. Said graph will be called `qcm_spectrum.pdf`. The file also generates a `local_dos_qcm.dat` file in `output/`. 

---

- The `fock_benchmark`, similarly to `graph.py`, will generate a graph of the Green Function using my implementation of the second quantization base. Said graph will be called `fock_spectrum.pdf`. The file also generates a `local_dos_fock.dat` file in `output/`. 


---

- The `combined_graphs.py` file, using the two or three `.dat` files in the `output/` directory, generates a combined graph of the Green Function calculated by 3 methods (ibmq,fock and qcm) or only 2 methods (ibmq and fock). This is meant to directly compare these methods.

---

- The `excitationNsites.def`, located in `excitation_files/`, is the file containing all the excitation types of the system. If this file is to be changed, it needs to respect the exact formatting of the current files.

---

All codes should be ran from the cmd. Also note that this code supports and needs the use of arguments: ```$ python3 <a_code.py> <arg1> ```. The argument is to be used when you want to run the code using the parameters stored in the provided examples. In other words, if you want to run, for example, the calculation of the function for a 2 sites system using the provided example, you would need to do the following:

- ```$ python3 generate.py 2``` 
- ```$ python3 graph.py 2```

The outputted graph will be stored in `examples/2sites/`.


## Execution

After setting up your parameters and, if applicable, your custom lattices and your custom circuit in `parameters.py`, refer to the following sections to run the code properly.

### Running the quantum computer simulation
Run these the first time:
- `generate.py`
- `graph.py`

Your `ibmq_spectrum.pdf` is stored in the same folder as `parameters.py`. As long as `graph.py` is provided matrices, the code can work idependantly or, in other words, without running `generate.py` before. 

### Running the qcm benchmark
Run this:
- `qcm_benchmark.py`

Your `qcm_spectrum.pdf` is stored in the same folder as `parameters.py`. Note that in the event the the ground_state is located in a block where the total spin is not of 0, you will need to modify this code to target the proper block for the calculation. However, this occurence is highly unlikely.

### Running the fock benchmark
Run this:
- `fock_benchmark.py`

Your `fock_spectrum.pdf` is stored in the same folder as `parameters.py`.

### Combining the graphs
Run this:
- `combined_graphs.py`

Your `combined_spectrum.pdf` is stored in the same folder as `parameters.py`.

### REMINDER
Don't forget the argument when you are running this code for the provided examples. For example, calculating for 4 sites using the provided example:
- `$ python combined_graph.py 4`


## TODO
- Add a code to generate excitation.def files automatically.
- Commenting all the codes, dear god.
- Optimisation
