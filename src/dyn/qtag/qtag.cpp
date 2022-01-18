/*********************************************************************************
* Copyright (C) 2022 Matthew Dutra, Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file qtag.cpp
  \brief The file implements various auxiliary functions for QTAG method
    
*/

#include "qtag.h"
#include "../gwp/libgwp.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libdyn namespace
namespace libdyn{

using namespace libgwp;

/// libqtag namespace
namespace libqtag{


CMATRIX qtag_psi(MATRIX q, MATRIX& q1, MATRIX& p1, MATRIX& alp1, MATRIX& s1, CMATRIX& Coeff){
/**
  Returns the (complex) wavefunction value *wf* at a given point *q* for all quantum states
  calculated using the TBF parameters stored in *q1, p1, alp1, s1* and coefficients *Coeff*.

  Args:
      q - MATRIX(ndof, 1) - point at which to evaluate the wavefunction
      q1 - MATRIX(ndof, ntraj) - coordinates of trajectories
      p1 - MATRIX(ndof, ntraj) - momenta of trajectories
      alp1 - MATRIX(ndof, ntraj) - Gaussian width parameters of trajectories
      s1 - MATRIX(ndof, ntraj) - Gaussian width parameters of trajectories
      Coeff - CMATRIX(nstates, ntraj) - amplitudes of all trajectories on all states

  Returns:
      wfc - CMATRIX(nstates, 1) - complex value of the wavefunction at q for all states

*/
  int idof, itraj, istate;
  int ndof = q1.n_rows;
  int ntraj = q1.n_cols;
  int nstates = Coeff.n_rows;

  CMATRIX wfc(nstates, 1);

  for(itraj<0; itraj<ntraj; itraj++){

    complex<double> prod(1.0, 0.0);  
    double sum1, sum2; sum1 = sum2 = 0.0;

    for(int idof=0; idof<ndof; idof++){
      double Q = q1.get(idof, itraj);
      double P = p1.get(idof, itraj);
      double A = alp1.get(idof, itraj);
      double S = s1.get(idof, itraj);
      double _q = q.get(idof, 0);
      double dq = _q - Q;

      prod *= pow((A/M_PI), 0.25);
      sum1 += A * dq * dq;
      sum2 +=  (P * dq + S);

    }// for idof

    double cs = cos(sum2);
    double si = sin(sum2);
    prod *= exp(-0.5 * sum1) * complex<double>(cs, si); 
    
    for(istate=0; istate<nstates; istate++){
      wfc.add(istate, 0,   Coeff.get(istate, itraj) * prod);
    }
  }// for itraj

  return wfc;

}


CMATRIX qtag_overlap_elementary(MATRIX q, MATRIX& p, MATRIX& alp, MATRIX& s){
/**
  Returns a ntraj x ntraj matrix of the GBFs, independent of their active states

  In the gwp, we have the phase parameters first, then the Gaussian width parameters 
  also, keep in mind the 0.5 factor 
    
*/

  MATRIX tmp(alp); tmp *= 0.5;

  return gwp_overlap_matrix(q, p, s, tmp, q, p, s, tmp );

}

CMATRIX qtag_kinetic_elementary(MATRIX q, MATRIX& p, MATRIX& alp, MATRIX& s, MATRIX& invM){
/**
  Returns a ntraj x ntraj matrix of the kinetic energy for GBFs, independent of their active states

  In the gwp, we have the phase parameters first, then the Gaussian width parameters 
  also, keep in mind the 0.5 factor 
    
*/

  MATRIX tmp(alp); tmp *= 0.5;

  return gwp_kinetic_matrix(q, p, s, tmp, q, p, s, tmp, invM );

}


CMATRIX qtag_overlap(vector<int>& active_states, CMATRIX& ovlp, int nstates){
/**
  Calculates the single-surface Hamiltonian matrix elements H_ij=<gi|KE+V|gj>, 
  computed using the basis parameters stored in *qpas*. This requires the single-surface 
  overlap matrix *ov*. Returns the single-surface Hamiltonian *H*.

  Args:
    active_states - list[ ntraj x int] - state indices for all trajectories
    ovlp - CMATRIX(ntraj, ntraj) - couplings of all trajectories as if they are on the same state

  Returns:
    S - CMATRIX(nstates x ntraj, nstates x ntraj) - The many-surface Hamiltonian matrix.

*/
  int i, j;
  int ntraj = active_states.size();

  CMATRIX S(nstates * ntraj, nstates * ntraj);

  for(int itraj=0; itraj<ntraj; itraj++){
    //for(int i=0; i<nstates; i++){
    i = active_states[itraj];

      int indx1 = itraj * nstates + i;

      for(int jtraj=0; jtraj<ntraj; jtraj++){
        //for(int j=0; j<nstates; j++){
          j = active_states[jtraj];
 
          int indx2 = jtraj * nstates + j;

          complex<double> val(0.0, 0.0);
          if(active_states[itraj] == active_states[jtraj]){
            val = ovlp.get(itraj, jtraj);
          }

          S.set(indx1, indx2, val);
        
//        }// j
      }// jtraj

//    }// for i
  }// for itraj

  return S;
}


CMATRIX qtag_hamiltonian(MATRIX q, MATRIX& p, MATRIX& alp, MATRIX& s, CMATRIX& Coeff,
                         vector<int>& active_states, CMATRIX& ovlp, CMATRIX& kin,
                         MATRIX& invM, nHamiltonian& ham, bp::object compute_ham_funct,  
                         bp::dict& compute_ham_params){
/**
  Calculates the single-surface Hamiltonian matrix elements H_ij=<gi|KE+V|gj>, 
  computed using the basis parameters stored in *qpas*. This requires the single-surface 
  overlap matrix *ov*. Returns the single-surface Hamiltonian *H*.

  Args:
    q - MATRIX(ndof, ntraj) - coordinates of the system
    p - MATRIX(ndof, ntraj) - momenta of the system
    alp - MATRIX(ndof, ntraj) - widths
    s - MATRIX(ndof, ntraj) - phases
    active_states - list[ ntraj x int] - state indices for all trajectories
    ovlp - CMATRIX(ntraj, ntraj) - couplings of all trajectories as if they are on the same state
    kin - CMATRIX(ntraj, ntraj) - kinetic energy of all trajectories as if they are on the same state
    invM - MATRIX(ndof, 1) - inverse masses for all DOFs
    ham - the Hamiltonian object that knows how to compute interactions for independent trajectories
    compute_ham_funct (PyObject) - the Python function that computes the model Hamiltonian
    compute_ham_params ( dict ) - the dictionary for the function above.

  Returns:
    H - CMATRIX(nstates x ntraj, nstates x ntraj) - The many-surface Hamiltonian matrix.

*/

  int ndof = q.n_rows;
  int nstates = Coeff.n_rows;
  int ntraj = active_states.size();
  int sz = nstates * ntraj;
  int i, j;

  CMATRIX H(sz, sz);

  //==================== Kinetic energy ================
  // Same ordering scheme as for super-overlap

  for(int itraj=0; itraj<ntraj; itraj++){    
    //for(int i=0; i<nstates; i++){
      i = active_states[itraj];

      int indx1 = itraj * nstates + i;

      for(int jtraj=0; jtraj<ntraj; jtraj++){
        //for(int j=0; j<nstates; j++){
          j = active_states[jtraj];
 
          int indx2 = jtraj * nstates + j;


          // Kinetic energy terms are always added
          H.set(indx1, indx2, kin.get(itraj, jtraj));
        
//        }// j
      }// jtraj

//    }// for i
  }// for itraj

  //==================== Potential energy/Couplings ================

  for(int itraj=0; itraj<ntraj; itraj++){    
    for(int i=0; i<nstates; i++){    
      int indx1 = itraj * nstates + i;

      for(int jtraj=0; jtraj<ntraj; jtraj++){
        for(int j=0; j<nstates; j++){ 
          int indx2 = jtraj * nstates + j;

          // Kinetic energy terms are always added
          H.set(indx1, indx2, kin.get(itraj, jtraj));

          // Potential term is only added for the same states
          complex<double> val(0.0, 0.0);
          if(active_states[itraj] == active_states[jtraj]){
            val = ovlp.get(itraj, jtraj);

            double pot = 0.0;
            // Here we need to add the call of a function to compute the 
            // potential = matrix elements for the mid-center of 2 trajectories
            // or otherwise

            H.add(indx1, indx2, pot*val);
          }
        
        }// j
      }// jtraj

    }// for i
  }// for itraj


  return H;


/*  TO BE IMPLEMENTED

//	qvals,pvals,avals,svals=MATRIX(qpas[0]),MATRIX(qpas[1]),MATRIX(qpas[2]),MATRIX(qpas[3])

	m_term=0.0
	for i in range(ndof):
		m_term+=(-1.0**2/(2.0*mass[i]))

	for i in range(ntraj):
		q1,p1,a1,s1=qvals.col(i),pvals.col(i),avals.col(i),svals.col(i)
		for j in range(ntraj):
			q2,p2,a2,s2=qvals.col(j),pvals.col(j),avals.col(j),svals.col(j)

			ke=m_term*gwp_kinetic(q1,p1,s1,a1/2,q2,p2,s2,a2/2)

			qpasi=[q1,p1,a1,s1];qpasj=[q2,p2,a2,s2]
			v=vcalc(ndof,qpasi,qpasj,model_params,nsurf,pot)
			H.set(i,j,ke+v*ov.get(i,j))
*/

}



CMATRIX qtag_momentum(MATRIX& q, MATRIX& p, MATRIX& alp, MATRIX& s, CMATRIX& Coeff){
/**
  Returns the momentum *mom* calculated for each basis function according to p=Im(grad(psi)/psi). 
  Also returns the corresponding real component *r*, which can be used in updates of the basis phase parameter *s*.

*/


  int i,j, idof;
  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  complex<double> dzt;
  double dq, nrm;

  CMATRIX mom(ndof,ntraj); // real = mom, imag = r
  CMATRIX dz(ndof,1);
  complex<double> one(0.0, 1.0);

  for(i=0; i<ntraj; i++){
    complex<double> z(0.0,0.0);
    dz *= 0.0;

    for(j=0; j<ntraj; j++){
      complex<double> nrm(1.0, 0.0);

      for(idof=0; idof<ndof; idof++){
          dq = q.get(idof, i) - q.get(idof, j);
          double argg = (p.get(idof, j) * dq + s.get(idof, j) );
          complex<double> tmp( cos(argg), sin(argg) );
          nrm *= exp( -0.5 * alp.get(idof, j) * dq*dq ) * tmp * pow( (alp.get(idof, i)/M_PI), 0.25);
     
      }// for idof

      z += Coeff.get(j) * nrm;

      for(idof=0; idof<ndof; idof++){
        dq = q.get(idof, i) - q.get(idof, j);
        dzt = -(alp.get(idof, j) * dq - one * p.get(idof, j)) * Coeff.get(j) * nrm;
        dz.add(idof, 0, dzt);
      }// for idof
      
    }// for j

    for(idof=0; idof<ndof; idof++){
      mom.set(i, idof,  dz.get(idof)/z );
    }// for idof

  }// for i

  return mom;

}



}// namespace libqtag
}// namespace libdyn
}// liblibra


