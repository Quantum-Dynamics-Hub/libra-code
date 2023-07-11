# Instructions 

  1. Create a directory "res"

  2. Run the calculations:  python run.py

  3. Plot the PES and population:   gnuplot pes.plt   and gnuplot pops.plt

  4. Copy files loop.plt and start.plt  into res folder and run gnuplot: gnuplot start.plt


# What to explore

## Explore various ways of setting up the initial wavefunction on the grid:

  1.  Case = 0  HO ground state wavefunction, |0>

  2.  Case = 1  HO 1-st excited state wavefunction, |1>

      In these examples the HO eigenstates correspond to the potential in which the present system evolves
      Therefore, there should be no overall motion of the wavepackets



## Create the initial state as a superposition of the HO eigenstates

  1.  Case = 2  a superposition of the first 2 HO eigenstates: ~ |0> + |1> (up to the overall normalization)

      Note that in this example, we will need to have 2 values of each parameter: x0, px0, alpha, initial electronic, initial vibrational state

      Because the superposition is no longer an eigenstate of the Hamiltonian, we should in principle observe some evolution


## Make the wavepackets moving 

  1.  Case = 3  the ground state wavefunction is now has a momentum, so we'll see some movement


## Displace the potential

  1.  Case = 4  we now change the minimum of the potential, so (unless we change the parameters of the 
      wavefunctions) accordingly, the initialized Gaussian is no longer an eigenfunction of the (new) Hamiltonian.
      So, we get some motion again



