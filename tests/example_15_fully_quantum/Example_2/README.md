# Instructions 

  1. Create a directory "res"

  2. Run the calculations:  python run.py

  3. Plot the PES and population:   gnuplot pes.plt   and gnuplot pops.plt

  4. Copy files loop.plt and start.plt  into res folder and run gnuplot: gnuplot start.plt


# What to explore

## Explore three model potentials/initial conditions (which counts to many more models):

  1.  Case = 1 *Double well potential*

  2.  Case = 2 *Decay of a metastable state*

  3.  Case = 3 *Modified decay of a metastable state*

      This is a model similar to the case 2, but without infinitely deel well to the right of the barrier
      It is still sufficiently deep so that the wavefunction would have little probability of going back



## Explore two options for measuring the rate of quantum transitions:

  1.  opt = 0  *Measuring the amount of probability dissapering in the complex absorbing potential*

      Pros:
      * can use smaller grid, which can help reduce the costs

      Cons:
      * the total population and energy of the system are not conserved
      * the dynamics may be delayed because the wavepacket would need to travel to the absorbing 
        potential (maybe much later than the moment of actual barrier crossing)
      * a partial reflection from the absorbing potential is possible, such that the dynamics may be
        perturbed


  2.  opt = 1  *Measuring the amount of probability passing through a point (dividing surface)*

      Pros:
      * the total population and energy of the system are conserved
      * the dynamics is not delayed because we measure exactly what we need 
        (observing the moment of actual barrier crossing)
      * no reflection from the absorbing potential, the dynamics is not perturbed

      Cons:
      * may need larger grid to avoid reflections/boundary passing
