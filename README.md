##README 


### What is this repository for?

* Designed to compute geodesics in Lie groups for quantum computation 

* Includes a templated Nelder Mead Optimiser



### How do I get set up?
 
* Requires the Armadillo matrix library, which are header files only. For downloads and installation instructions see http://arma.sourceforge.net/docs.html

* Source code for solvers is in /solvers 

* /src contains templated Nelder Mead solver, code to generate Pauli matrices, and an implementation of the Nelder Mead for the unitary group 

* Make files for serial and parallel exist in /serialSolver and /segmentedSolver respectively. Build and run in these files. Makefiles build code found in /solvers. 

### Running the code?

* "./<exename> <Max No. Nelder Mead steps> <No. of steps in Cayley Integrator> <tolerance>  <Max No. Leap-frog steps> "

### Contribution guidelines


### Who do I talk to?

* meswaddle@gmail.com 
