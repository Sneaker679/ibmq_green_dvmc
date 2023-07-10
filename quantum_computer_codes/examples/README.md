# Examples
Currently, there are two configurations of lattice in this directory. 

## Basic execution
To run each of them, make sure you are in the directory corresponding to the desired configuration (`2sites/` or `4sites/`), and run the following lines to run the entirety of this code and get all the results:
```bash
$ python3 ../../generate.py
$ python3 ../../graph.py
$ python3 ../../fock_benchmark.py
$ python3 ../../combined_graphs.py
```

In your current working directory, where `parameters.py` is located, new files should have appeared, as well a a `output/` directory. You can ignore the previous directory as the actual results are the PDFs that appeared in your current directory. In fact, there will be 3 of them: `ibmq_spectrum.pdf`, `fock_spectrum.pdf` and `combined_spectum.pdf`, the latter being a mix of both previous ones.


## Basic execution - with pyqcm
If qcm is installed and correctly modified according to the procedure described in the main README, change use_qcm in `parameters.py` for 'Y'. Then, you can also run this:
```bash
$ python3 ../../qcm_benchmark.py
$ python3 ../../combined_graphs.py
```

A new file (PDF) will appear in your current working directory: `qcm_benchmark.pdf`. Furthermore, `combined_spectrum.pdf` will be updated to include the qcm benchmark as well.
