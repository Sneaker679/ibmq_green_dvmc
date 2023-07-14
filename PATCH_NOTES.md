# PATCH NOTE

## Bug list and their fix
Bug: The 3 ways for calculating the Green Function were not agreeing at certain mu values.
The issue was: The ground state fed to the program didn't take into account spin symmetry. Thus, depending on the block containing the ground state, said ground state would have reverse spins compared to the operators used to calculating the Green Function. These would yield trivial values and other unusualresults.
The fix: To fix this, all benchmarks were extensively modified so that everytime, they would select the exact same ground state, from the exact same block of the hamiltonian.

Bug: When QCM benchmark, despite being our most reliable source for benchmarking, was yielding different results compared to the other 2 benchmarks at specific mu values.
The issue was: When QCM finds out that the ground state is located in many blocks of the hamiltonian, it somehow decides to combine these results to calculate the Green Function. In other words, it would take part of one ground state and part of another ground state to build the spectrum. This behavior is not implemented at all with the other benchmarks.
The fix: When QCM finds multiple blocks for a single ground state energy, the ground state is recalculated with only one of those blocks. Priority is given to the block with the spin closest to Sz=0.

Bug: Multiprocessing is slower than singleprocessing.
The issue was: Passing down the hamiltoninan to all processes was very ressource intensive.
The fix: Share the memory for the hamiltonian using the mpire module.

Bug: Peaks were missing from the spectrum, even in cases where it should be symmetrical with the spectrum of another mu that doesn't have missing peaks.
The issue: The graph.py program was doing too much filtering with the e_ac_r,u_ac_r,e_ca_r and u_ca_r variables.
The fix: Have the program use directly the unfiltered e_ac and ... variables. Now, all the peaks are there.
