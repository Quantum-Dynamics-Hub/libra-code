# Libra

[![Build Status](https://travis-ci.org/Quantum-Dynamics-Hub/libra-code.svg?branch=master)](https://travis-ci.org/Quantum-Dynamics-Hub/libra-code)

This is the main page of the computational chemistry methodology discovery library, Libra

## About

* [extensive tutorials, example, and documentation](https://github.com/compchem-cybertraining/Tutorials_Libra).
This resource is actively developed, but the older tutorials may not fully reflect the current stage of the code itself, so the 
tutorials demonstrating the older features may be failing - please let us know if you need to use any of those and run into problems.

* [Read-the-Docs documentation](https://libra-documentation.readthedocs.io/en/latest/)
This is far from being complete. The work in progress...

* [old program website](https://quantum-dynamics-hub.github.io/libra/index.html)
Not maintained for long time, but may still contain some useful info on some topics

## Community and Support

* Click below to join our community:
👉 [Join Slack](https://join.slack.com/t/quantumdynamicshub/shared_invite/zt-mjbhjssx-GGhsbYHxeBMvhmumK_j7LA)

* Open an Issue - to ask questions or report a problem

* The following [public forum](https://groups.google.com/forum/#!forum/quantum-dynamics-hub) exists, but
  hasn't been used for a while - use the Slack workspace instead


## Installation 

Please see [here](INSTALLATION.md)

## Developers and Contributors

  * Dr. Alexey Akimov (University at Buffalo, [link](https://akimovlab.github.io/index.html) )  
      The main developer and maintainer of the code
  
  * Dr. Daeho Han (University at Buffalo)
      Implementation of the exact-factorization (XF) methods, quantum trajectories surface hopping (QTSH),
      validating and testing Ehrenfest dynamics and other internals, implementing some model Hamiltonians 

  * Mr. Brendan Smith (University at Buffalo) 
      Entangled trajectories Hamiltonian, NA-MD with spin-orbit coupling, NBRA workflows, BL-LZ NA-MD tutorials and examples, 
      Libra/DFTB+, Libra/QE, Libra/ErgoSCF, Libra/CP2K, and Libra/Gaussian interfaces
      
  * Mr. Mohammad Shakiba (Shahid Bahonar University of Kerman, Iran)
      Cube file processing scripts, Libra/CP2K and Libra/Gaussian, Libra/Libint2 interfaces

  * Mrs. Story Temen (University at Buffalo)
      Implementation and testing of the HEOM codes

  * Dr. Wei Li (Hunan Agricultural University)
      NA-MD with spin-orbit coupling

  * Dr. Kosuke Sato (Toyota Research Lab) 
      State reordering scripts, Libra/GAMESS interface (Libra-X)

  * Dr. Ekadashi Pradhan (York University) 
      Libra/QE interface, delta-SCF NA-M (Libra-X)

  * Dr. Amber Jain (Indian Institute of Technology Bombay, India)
      Implementation and testing of the HEOM codes

  * Dr. Xiang Sun (NYU Shanghai, China)
      Implementation and testing of the FGR codes

  * Dr. Sophya Garashchuk (University of South Carolina)
      QTAG theory development

  * Dr. Matthew Dutra (University of South Carolina)
      Implementation and testing of the QTAG codes

## References

This code is provided in the hope it will be useful. 

  If you use Libra in your research, please cite the following paper:

  ### The "main" Libra papers 
  * [More recent overview of Libra's capabilities](https://doi.org/10.1016/j.simpa.2022.100445)
  Shakiba, M.; Smith, B.; Li, W.; Dutra, M.; Jain, A.; Sun, X.; Garashchuk, S.; Akimov, A.V.* "Libra: A modular software library for quantum
  nonadiabatic dynamics" *Software Impacts* **2022**  14, 100445
   
  * [The initial implementation](http://onlinelibrary.wiley.com/doi/10.1002/jcc.24367/full)
  Akimov, A. V. "Libra: An open-Source 'methodology discovery' library for quantum and classical dynamics simulations" 
  *J. Comput. Chem.*  **2016**  37, 1626-1649

  If you use any of the Libra's methods or implementation, please cite as appropriate:

  ### Papers that describe various developments in Libra
  * [FISH: Fully-integrated surface hopping](https://doi.org/10.1021/acs.jpclett.5c01509)
  Han, D.; Shakiba, M.; Akimov, A. V. "Fully-Integrated Surface Hopping as Quantum Decoherence Correction in Nonadiabatic Dynamics"
  *J. Phys. Chem. Lett.* **2025** 16, 7168-7176

  * [Multiple-state QTSH](https://doi.org/10.1021/acs.jctc.4c01751)
  Han, D.; Martens, C. C.; Akimov, A. V. "Generalization of Quantum-Trajectory Surface Hopping to Multiple Quantum States"
  *J. Chem. Theory Comput.* **2025** 20, 5022-5042

  * [F-tracking of excited states](https://pubs.acs.org/doi/10.1021/acs.jpclett.4c02909)
   Akimov, A. V. "State Tracking in Nonadiabatic Molecular Dynamics Using Only Forces and Energies"
   *J. Phys. Chem. Lett.* **2024** 15, 11944-11953
    
  * [FSSH-3, modified FSSH-2, implmenentation of GFSH](https://www.tandfonline.com/doi/full/10.1080/00268976.2024.2376893)
  Akimov, A. V. "The Fewest Switches Surface Hopping as An Optimization Problem" *Mol. Phys.* **2024** e2376893

  * [Exact Factorization methods: SHXF, MQCXF, MFXF](https://doi.org/10.1021/acs.jctc.4c00343)
   Han, D.; Akimov, A.V. "Nonadiabatic Dynamics with Exact Factorization: Implementation and Assessment."
   *J. Chem. Theory Comput.* **2024** 20, 5022–5042

  * [ML Hamiltonian mapping approach](https://doi.org/10.1021/acs.jctc.4c00008)
  Shakiba, M.; Akimov, A.V. "Machine-Learned Kohn–Sham Hamiltonian Mapping for Nonadiabatic Molecular Dynamics."
  *J. Chem. Theory Comput.* **2024** 20, 2992–3007

  * [TC-NBRA](https://doi.org/10.1021/acs.jpclett.3c03029)
  Akimov, A.V. "Energy-Conserving and Thermally Corrected Neglect of Back-Reaction Approximation Method for Nonadiabatic Molecular Dynamics"
  *J. Phys. Chem. Lett.* **2023** 14, 11673-11683

  * [QTAG](https://doi.org/10.1002/qua.27078)
  Dutra, M.; Garshchuk, S.; Akimov, A. "The Quantum Trajectory-guided Adaptive Gaussian Methodology in the Libra Software Package"
  *Int. J. Quntum Chem.* **2023.** 123, e27078

  * [(generalized) Local diabatization](https://doi.org/10.1007/s00214-023-03007-7)
  Shakiba, M.; Akimov, A.V. "Generalization of the Local Diabatization Approach for Propagating Electronic Degrees of Freedom in Nonadiabatic Dynamics"
  *Theor. Chem. Acc.* **2023** 142, 68

  * [xTB from cp2k/Libra interface for NA-MD in large-scale systems, k-point formulation, spin-adaptation](https://doi.org/10.1021/acs.jctc.2c00297)
  Shakiba, M.; Stippel, E.; Li, W.; Akimov, A. V. "Nonadiabatic Molecular Dynamics with Extended Density Functional Tight-Binding:
  Application to Nanocrystals and Periodic Solids" *J. Chem. Theory Comput.* **2022** 18, 5157-5180

 * [many-body effects, NA-MD with TD-DFT states](https://pubs.acs.org/doi/10.1021/acs.jctc.0c01009)
  Smith, B.; Shakiba, M.; Akimov, A. V. "Nonadiabatic Dynamics in Si and CdSe Nanoclusters: Many-Body vs. Single-Particle
  Treatment of Excited States" *J. Chem. Theory. Comput.* **2021** 17, 678-693

  * [HEOM implementation](https://pubs.acs.org/doi/10.1002/qua.26373)
  Temen, S.; Jain, A.; Akimov, A. V. "Hierarchical equations of motion in the Libra software package"
  *Int. J. Quant. Chem.* **2020** 120, e26373

  * [Belyaev-Lebedev-Landau-Zener Surface Hopping within the Neglect of Back-Reaction Approximation](https://pubs.acs.org/doi/10.1021/acs.jpclett.9b03687)
  Smith, B.; Akimov, A. V. "Hot Electron Cooling in Silicon Nanoclusters via Landau-Zener Non-Adiabatic Molecular Dynamics: 
  Size Dependence and Role of Surface Termination" *J. Phys. Chem. Lett.* **2020** 11, 1456-1465 

  * [Phase correction, Ehrenfest dynamics details, basis transformations (see the SI)](https://doi.org/10.1021/acs.jpclett.8b02826)
  Akimov, A. V.; "A Simple Phase Correction Makes a Big Difference in Nonadiabatic Molecular Dynamics"
  *J. Phys. Chem. Lett.*  **2018** 9, 6096-6102



  You may find the following papers useful examples
  ### Papers that utilize Libra

  * [Formulation of a fragment-based NA-MD](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00955)
  Akimov, A. V. "Nonadiabatic Molecular Dynamics with Tight-Binding Fragment Molecular Orbitals" 
  *J. Chem. Theory Comput.*  **2016** 12, 5719-5736

  * [Quasi-stochastic Hamiltonian for longer NA-MD](http://pubs.acs.org/doi/abs/10.1021/acs.jpclett.7b02185)
  Akimov, A. V.; "Stochastic and Quasi-Stochastic Hamiltonians for Long-Time Nonadiabatic Molecular Dynamics"
  *J. Phys. Chem. Lett.*  **2017** 8, 5190-5195

  * [Entrangled-trajectories Hamiltonian dynamics to capture quantum effects of nuclei](https://doi.org/10.1063/1.5022573)
  Smith, B. A.; Akimov, A. V. "Entangled trajectories Hamiltonian dynamics for treating quantum nuclear effects" 
  *J. Chem. Phys.* **2018** 148, 144106

  * [Inclusion of the Spin-orbit coupling in NA-MD](https://doi.org/10.1021/acsenergylett.8b01226)
  Li, W.; Zhou, L.; Prezhdo, O. V.; Akimov, A. V. "Spin-Orbit Interactions Greatly Accelerate Nonradiative Dynamics
  in Lead Halide Perovskites" *ACS Energy Lett.* **2018** 3, 2159-2166




