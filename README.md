# ibmq_green_dvmc
## Table of Contents
- [Installation](#installation)
- [Description](#Description)
    - [Usage](#usage)
    - [How it works](#how-it-works)
    - [Execution](#execution)
- [Todo](#todo)

## Installation
This is a Python code. Use the latest version of language as of 15/06/2023 to avoid issues.
To run this code, you will need to install the following packages for your Python. The links to their documentation is also provided :
- qiskit-->  (https://qiskit.org/documentation/getting_started.html)
- qcm-->  (https://qcm-wed.readthedocs.io/en/latest/intro.html#installation)

Note that the "qcm" package is not manditory as its only purpose is to benchmark the calculations using the quantum computer.

Also note that in its current form, the "qcm" package needs a custom modification so it can work with the current code. In the "qcm" package, locate `pyqcm/_spectral.py` and find in it the "*cluster_spectral_function*" function. Remove the hashtags preceding the "*return*" at the end. You will need to rebuild the package and perform all the installation steps again, as described by the documentation.

## Description
This code is meant to calculate the Green Function for a specified system of low amount of sites. Instead of calculating it using conventional binary based computers, we simulate the calculation on a quantum computer using the qiskit package. 

The code itself is accompanied by a few supplementary benchmarks codes. These benchmark codes come in two types : benchmarks that use a custom second quantization base implementation, and benchmarks that use the "qcm" package. All the codes related to this custom second quantization base implementation are located in the *second_quantization_codes/* directory, while those from the "qcm" package are located in `quantum_computer_codes/...` and its subsequent directories. This is because the quantum computer is directly benchmarked with "qcm".

## Usage
The primary use of this code is to calculate the Green Function. I will walktrough you into using the included preset for a 2 sites system. The general idea is the same for other sites. Note that the usage of the second quantization benchmarks will not be explained.

### How it works
All the files located in `quantum_computer_codes/2sites/` are useful to the calculation of the Green Function for 2 sites systems. Here is the uses of each individual files:

---
- The `parameters.py` file is where you can modify the parameters and variables to be used by all the files in the folder.
- The `generate.py` is meant to generate the specified matrices, print them in the CMD and save them as `.npy` files. This file needs to be modified depending on the wanted structure of the lattice. Specifications on the inner workings of the lattice objects are located in the "qiskit" package documentation. The quantum circuit used in the code may also need to be modified depending on your needs.
- The `gs.py` file calculates the ground state energy of the system using the quantum computer.
- The `graph.py` file will use the `.npy` files and the ground state energy generated by `gs.py` to calculate the Green Function and make a graph with it. Said graph will be called `ibmq_spectrum.pdf`. This file also generates a ```local_dos.dat``` file.
- The `generate_graph.py` file will run `generate.py` and then `graph.py` for an easy usage of the code. Note that those last two files can be ran independantly with no problems.

---

- The `qcm_benchmark.py`, similarly to `graph.py`, will generate a graph of the Green Function using the "qcm" package. Said graph will be called `qcm_spectrum.pdf`. The file also generates a ```local_dos_qcm.dat``` file. 
- The `qcm_model_2D_C4.py` file is used to build the lattice of the system. This needs to be manually modified depending on the wanted structure of the lattice so the the "qcm" package can properly benchmark the quantum computer implementation.

---

- The `combined_graphs.py` file, using the two `.dat` files in the directory, generates a combined graph of the Green Function calculated by both methods (ibmq and qcm). This is meant to directly compare these methods.

---

- The `excitation2sites.def` is the file containing all the excitation types of the system. If this file is to be changed, it needs to respect the exact formatting of the current files.

---
### Execution
For a basic utilisation of the code, first go in this subdirectory : `quantum_computer_codes/...`. From there, directories should already be there for a few specific systems. These directories will be named like this : (# of sites)sites/. Enter one of them, then open the `parameters.py` file. Modify the parameters and variables as you see fit, then save the file and run the `generate_graph.py` code. Wait for the calculation to end. Finally, open the `ibmq_spectrum.pdf` for the result.

### Custom systems
***This section is unfinished.
To create and calculate the green function for custom systems, you will need to copy one of the examples and modify the following files to your needs :
- `parameters.py` : change N and N_min.
- `generate.py` : First modify the lattice structure to match your desired system. Refer to qiskit's documentation for this. Also make sure the file has been copied from `4sites/` and not `2sites/`, or else you have to modify the quantum circuit to your needs.
- `gs.py` : Make sure the lattice structure used in this file is the same as `generate.py`.
- `qcm_model_*.py` : Make sure the proper lattice structure is setup for your needs. Refer to the qcm documentation.
- `excitation.def` : Make sure the excitation types there match the number of sites and the number of excitations. The number of excitations is calculated using the following formula : N_exc = 2 + N_min ( N_min + 1 ).

## TODO
- Implementation of the "t" matrix instead of just a number for the second quantization codes. Note that because this is not implemented yet, the second quantization codes cannot benchmark the quantum computer's calculation for more than 2 sites unless the fermions can jump to *any* site.
- Implementation of a custom quantum circuit for 4 sites instead of calculating its ground state using numpy.
- Testing for 3, 5 and more sites and confirming the steps described in the "Custom_systems" section.
- Add a code to generate excitation.def files automatically.
- Have "qcm" figure out the location of the ground state instead of assuming it is at N:Nc,S:0.
- Add to `parameters.py` options to create a lattice and have it match in both qcm and qiskit. 
- Optimisation of the code in general, including only having to generate the hamiltonian once instead of twice.

