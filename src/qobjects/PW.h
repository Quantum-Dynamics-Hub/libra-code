/*********************************************************************************
* Copyright (C) 2013-2106 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file PW.h
  \brief The file describes data types for operating on plane waves (with general k-points)
    
*/

#ifndef PW_H
#define PW_H

#include "../math_linalg/liblinalg.h"
#include "../math_specialfunctions/libspecialfunctions.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libspecialfunctions;

namespace libqobjects{

// Some useful (potentially) macros
#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define SWAP_4(x) ( ((x) << 24) | (((x) << 8) & 0x00ff0000) | \
         (((x) >> 8) & 0x0000ff00) | ((x) >> 24) )
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_LONG(x) (*(unsigned *)&(x) = SWAP_4(*(unsigned *)&(x)))
#define FIX_FLOAT(x) FIX_LONG(x)

template<class T>
void get_value(T& val,char* memblock,int& pos)
{ memcpy((void*)&val, &memblock[pos], sizeof(T)); pos += sizeof(T); }


class PW{   
/**
  This class describes the data type for storing general MO in the plane wave basis

  psi_{i,k}(r) = 1/sqrt(V) *sum c_i(G) * exp(-i*(G+k)*r), where we run over all G points in the expansion
                             G

  this wavefunction accounts both the band index (i) and the k-point index (k)
  
*/

public:
  double energy;      ///< energy of this orbital (eigenvalue)
  double gweight;     
  double fweight;

  int iband;          ///< band (MO) index

  int kx,ky,kz;       ///< indexes of this kpoint
  VECTOR k_point;     ///< actual k-vector

  int npw;                         ///< number of plane waves in MO expansion
  vector< complex<double> > coeff; ///< coefficients of the plane waves: c(G) = exp(-i*(g+k)*r), where we run over all G points


  // Constructor
  PW(int _i, int _kx, int _ky, int _kz, int _npw){ 
    energy = 0.0; gweight = 0.0; fweight = 0.0;
    iband = _i; kx = _kx; ky = _ky; kz = _kz;
    npw = _npw;
    coeff = vector<complex<double> >(npw,complex<double>(0.0,0.0));    
  }

  // Copy constructor
  PW(const PW& plw){
    energy = plw.energy; gweight = plw.gweight; fweight = plw.fweight;
    iband = plw.iband;
    kx = plw.kx; ky = plw.ky; kz = plw.kz; k_point = plw.k_point; 
    npw = plw.npw; coeff = plw.coeff;
  }
  //PW& operator=(const PW&);

  // Destructor
  ~PW(){ if(coeff.size()>0){ coeff.clear(); } npw = 0; }


  // Methods
  PW conj();
  void normalize();
  void complete();  // add complex conjugate part of the wavefunction


  // Operators
  PW operator-();                 // Negation;
  PW operator+(const PW& m);
  PW operator-(const PW& m);
  void operator+=(const PW& m);
  void operator-=(const PW& m);
  PW operator/(double num);
  PW operator/(complex<double> num);

  // Friends
  friend PW operator*(const double& f, const PW& m1);  // Multiplication of PW and double;
  friend PW operator*(const PW& m1, const double& f);  // Multiplication of PW and double;
  friend PW operator*(const float& f, const PW& m1);   // Multiplication of PW and float;
  friend PW operator*(const PW& m1, const float& f);   // Multiplication of PW and float;
  friend PW operator*(const complex<double>& f, const PW& m1);  // Multiplication of PW and complex<double>;
  friend PW operator*(const PW& m1, const complex<double>& f);  // Multiplication of PW and complex<double>;
  friend PW operator*(const complex<float>& f, const PW& m1);   // Multiplication of PW and complex<float>;
  friend PW operator*(const PW& m1, const complex<float>& f);   // Multiplication of PW and complex<float>;

  friend complex<double> operator*(const PW& m1, const PW& m2);  // Multiplication of PW and PW


  friend int operator == (const PW& pw1, const PW& pw2){
    int res = ((pw1.iband==pw2.iband) && (pw1.energy==pw2.energy) && (pw1.kx==pw2.kx)
              && (pw1.ky==pw2.ky) && (pw1.kz==pw2.kz) && (pw1.npw==pw2.npw)            );

    if(res){ 
      for(int i=0;i<pw1.npw;i++){  res *= (pw1.coeff[i]==pw2.coeff[i]); }
    }
    return res;
  }



};

typedef std::vector<PW> PWList; ///< This is the data type for representing vector of PW objects


complex<double> I_1D(double kx, double kxp, double gx, double gxp);
complex<double> I_3D(VECTOR& k, VECTOR& kp, VECTOR& g, VECTOR& gp);
CMATRIX pw_overlap(VECTOR& k1, VECTOR& k2, CMATRIX& coeff1, CMATRIX& coeff2, vector<VECTOR>& grid1, vector<VECTOR>& grid2);

CMATRIX QE_read_acsii_wfc(std::string filename, int kpt, vector<int>& act_space, int verbose);
MATRIX QE_read_acsii_grid(std::string filename, int verbose);
vector<CMATRIX> compute_Hprime(CMATRIX& wfc, MATRIX& grid, MATRIX& reci);

/*
class wfc{

  void aux_line2vec(string line,vector<double>& a);

public:
  // Info
  int nspin;                // type of the spin-polarization used in calculations
  int gamma_only;           // gamma trick:  =1 (true), =0 (false)
  int natoms;               // number of atoms
  double tpiba;             // units of the lattice vectors (reciprocal)
  double alat;              // units of the lattice vectors (real)
  double omega;             // volume of the unit cell
  double efermi;            // fermi energy
  std::string cell_units;   // units of the simulation cell
  std::string energy_units; // units of the energy (eigenvalues)
  vector<double> a1,a2,a3;  // lattice vectors
  vector<double> b1,b2,b3;  // reciprocal lattice vectors

  // Data
  int nkpts;                // number of k-points
  int nbands;               // number of bands - the same for all k-points
  int npw;                  // number of plane waves in each mo - same for all mos
  int is_allocated;         // flag that says which level of memory is allocated 
  vector<K_point> kpts;     // array of k-points
  vector<vector<int> > grid;// internal vector<int> consists of 3 integers

  

  // Constructor
  wfc(){ nkpts = 0; nbands = 0; npw = 0; is_allocated = 0; }
  wfc(int _nkpts){ nkpts = _nkpts; nbands = 0; npw = 0; kpts = vector<K_point>(nkpts,K_point()); is_allocated = 1; }
  wfc(int _nkpts,int _nbnds){ nkpts = _nkpts; nbands = _nbnds; npw = 0; kpts = vector<K_point>(nkpts,K_point(nbands)); is_allocated = 2; }
  wfc(int _nkpts,int _nbnds,int _npw){nkpts = _nkpts; nbands = _nbnds; npw = _npw; kpts = vector<K_point>(nkpts,K_point(nbands,npw)); is_allocated = 3; }
  wfc(wfc& wfc1,int min1,int max1, wfc& wfc2,int min2,int max2);

  // Destructor
  ~wfc(){  if(kpts.size()>0){ kpts.clear(); } if(grid.size()>0){ grid.clear();} nkpts = 0; nbands = 0; npw = 0; is_allocated = 0; }

  // Copy constructor
  //wfc(const wfc&);

  // Methods:
  // QE methods
  void QE_read_binary_wfc(std::string filename,int,int,int);
  void QE_read_acsii_wfc(std::string filename);
  void QE_read_acsii_grid(std::string filename);
  void QE_read_acsii_index(std::string filename);

  // Common methods
  void complete();
  void normalize();
  void transform(int k,matrix&);
  void restore(int k1,int do_complete);
  void set_latt_vectors(boost::python::list);
  void set_reci_vectors(boost::python::list);
  void compute_Hprime(int minband,int maxband,std::string filename);
};


// Functions using the arguments of wfc type
void overlap(wfc& wfc1,int k1,int minband,int maxband,std::string filename);
void energy(wfc& wfc1,int k,int minband,int maxband,std::string filename);
void nac(wfc& wfc1,wfc& wfc2,int k1,int k2,int minband,int maxband,double dt,std::string filename);
void ham(wfc& wfc1,wfc& wfc2,int k1,int k2,int minband,int maxband,double dt,std::string filename);

*/


}// namespace libqobjects
}// namespace liblibra


#endif // PW_H

