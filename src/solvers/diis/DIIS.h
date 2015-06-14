#ifndef DIIS_H
#define DIIS_H

#include "Mathematics.h"



class DIIS{

public:

  DIIS(int,int);

  void update_diis_coefficients();
  void add_diis_matrices(MATRIX*, MATRIX*);

  void extrapolate_matrix(MATRIX*);




  int N_diis_max;                // Length of DIIS history (size of the lists)  

  int N_diis;                    // current # of matrices stored
  int N_diis_eff;                // effective size of the DIIS matrix (such that it is full-rank)
  vector<MATRIX*> diis_X;        // diis iteration of objective matrices (typically Fock matrices)
  vector<MATRIX*> diis_err;      // diis error matrices
  vector<double>  diis_c;        // diis extrapolation coefficients


};


#endif // DIIS_H 