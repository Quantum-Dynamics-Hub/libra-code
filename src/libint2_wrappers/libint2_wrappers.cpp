/*********************************************************************************
* Copyright (C) 2021 Mohammad Shakiba and Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
* This file uses Libint library (https://github.com/evaleev/libint) for computing 
* the atomic orbitals overlaps for two different geometries. This can be used to
* compute the molecular orbital overlaps for non-adiabatic molecular dynamics.
* The code uses Pybid11 to interface the code with Python.
*
* This code maily uses the test files in 
* (https://github.com/evaleev/libint/tree/master/tests/hartree-fock) with 
* modifications to be used for computing the atomic orbitals overlap.
*********************************************************************************/

/**
  \file libint2_wrappers.cpp
  \brief The file describes the classes and functions that interface libint2 functions and types with those of Libra and expose them to Python
    
*/


#include "libint2_wrappers.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libspecialfunctions;


namespace liblibint2_wrappers{


using libint2::Shell;
using libint2::Engine;
using libint2::Operator;
using std::cout;
using std::endl;



// We first need to initailize a libint2::Shell. Here is how we do it.
std::vector<libint2::Shell> initialize_shell(int l_val, bool is_spherical, 
 const std::vector<double>& exponents, const std::vector<double>& coeff, VECTOR& coords){

  std::vector<Shell> shells;

  // The same as above explanation for exponents and contraction coefficients
  libint2::svector<double> q;
  libint2::svector<double> p;
  for (int i = 0; i < exponents.size(); i++) {
       q.push_back(exponents[i]);
       p.push_back(coeff[i]);
      }

   // Append to shell
  shells.push_back(
                       {
                         q, // exponents of primitive Gaussians
                         {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                           {l_val, is_spherical, p}
                         },
                         {{coords.x, coords.y, coords.z}}   // origin coordinates
                       }
                   );

  return shells;
}




// This function adds the basis sets for each atomic type to a libint2::Shell variable
void add_to_shell(std::vector<libint2::Shell>& shells, 
  int l_val, bool is_spherical , const std::vector<double>& exponents, const std::vector<double>& coeff, VECTOR& coords){

    // We need libint2::svector (take a look at shell.h file in libint library folder)
	// q for exponents
    libint2::svector<double> q;
    // p for contraction coefficients (libint will automatically converts them to refer to 
	// normalization free primitives before computing integrals) (Shell::renorm in shell.h)
    libint2::svector<double> p;
	
    // Appending the exponents and contraction coefficients from Python to q and p respectively.
    for (int i = 0; i < exponents.size(); i++) {
        q.push_back(exponents[i]);
        p.push_back(coeff[i]);
    }

    // The spherical/cartesian flag for shell
    shells.push_back(
                     {
                       q, // exponents of primitive Gaussians
                       {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                         {l_val, is_spherical, p}
                       },
                       {{coords.x, coords.y, coords.z}}   // origin coordinates
                     }
                    );
    // Return the shell with new data
    return shells;

}

// This function prints the shells. It is useful to check the results (this is also interfaced with Python using Pybind11)
void print_shells(std::vector<libint2::Shell>& shells){

    std::cout << "\n\tShells are:\n";
    for(auto s=0; s < shells.size(); s++){
      std::cout << shells[s] << std::endl;
      // The shell size
      std::cout << "\n The shell size is:\n" << shells.size() << std::endl;
    }
}



// Below are the functions that we use for starting the engine for computing the overlap integrals.
size_t nbasis(const std::vector<libint2::Shell>& shells){
  size_t n = 0;
  for(int i=0; i<shells.size(); i++){
    n += shells[i].size();
  }

  return n;
}

size_t max_nprim(const std::vector<libint2::Shell>& shells){
  size_t n = 0;
  for(int i=0; i<shells.size(); i++){
    n = std::max(shells[i].nprim(), n);
  }

  return n;
}


int max_l(const std::vector<libint2::Shell>& shells){
  int l = 0;
  for(int i=0; i<shells.size(); i++){
    for(int j=0; j<shells[i].contr.size(); j++){
      l = std::max(shells[i].contr[j].l, l);
    }
  }

  return l;
}

std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells){
  std::vector<size_t> result;
  result.reserve(shells.size());

  size_t n = 0;
  for(int i=0; i<shells.size(); i++){
    result.push_back(n);
    n += shells[i].size();
  }

  return result;
}



/*

// Computing the integrals between two shells serial
MATRIX compute_1body_ints(const std::vector<libint2::Shell>& shells_1, const std::vector<libint2::Shell>& shells_2,libint2::Operator obtype)
{

  const auto n_1 = nbasis(shells_1);
  const auto n_2 = nbasis(shells_2);
  //std::cout << "nbasis (shells_1)" << n_1 << "\n";
  //std::cout << "nbasis (shells_2)" << n_2 << "\n";
  Matrix result(n_1,n_2);
  int max_n;
  int max_lval;
  if (max_nprim(shells_1)>max_nprim(shells_2))
  {
    max_n = max_nprim(shells_1);
  }
  else
  {
    max_n = max_nprim(shells_2);
  }
  if (max_l(shells_1)>max_l(shells_2))
  {
    max_lval = max_l(shells_1);
  }
  else
  {
    max_lval = max_l(shells_2);
  }

  Engine engine(obtype, max_n, max_lval, 0);

  auto shell2bf_1 = map_shell_to_basis_function(shells_1);
  auto shell2bf_2 = map_shell_to_basis_function(shells_2);

  // buf[0] points to the target shell set after every call  to engine.compute()
  const auto& buf = engine.results();

  // loop over unique shell pairs, {s1,s2} such that s1 >= s2
  // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
  //std::cout << shells_1.size();
  for(auto s1=0; s1!=shells_1.size(); ++s1) {

    auto bf1 = shell2bf_1[s1]; // first basis function in this shell
    //std::cout << "Flag bf1: " << bf1 << "\n";
    auto n1 = shells_1[s1].size();

    //std::cout << "shells_1[s1]" << shells_1[s1].size() << "\n";

    //for(auto s2=0; s2<=s1; ++s2) {
    for(auto s2=0; s2!=shells_2.size(); ++s2) {

      auto bf2 = shell2bf_2[s2];
      //std::cout << "Flag bf2: " << bf2 << "\n";
      auto n2 = shells_2[s2].size();

      //std::cout << "shells_2[s2]" << shells_2[s2].size() << "\n";

      // compute shell pair
      engine.compute(shells_1[s1], shells_2[s2]);
      //std::cout << "Flag after engine.compute" << "\n";
      // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
      Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
      //std::cout << "Flag after buf_mat" << "\n" << "buf_mat.size: \n" << buf_mat.size() << "\n";
      //std::cout << "Flag after buf_mat" << "\n" << "buf_mat  \n" << buf_mat << "\n";
      //std::cout << "Flag after buf_mat" << "\n" << "bf1.size" << bf1.size() << "\n";
      //std::cout << "Flag after buf_mat" << "\n" << "bf2.size" << bf2.size() << "\n";
      result.block(bf1, bf2, n1, n2) = buf_mat;
      //std::cout << "Flag after result.block" << "\n";
      if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
      result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
      //std::cout << "Flag after buf_mat.transpose" << "\n";
    }
  }
  //std::cout << "Flag before MATRIX res" << "\n";
  MATRIX res(n_1, n_2);
  for(int i=0;i<n_1;i++){
    for(int j=0; j<n_2;j++){
      res.set(i,j, result(i,j));
    }
  }


  return res;
}

*/


// Computing the integrals between two shells in parallel using OpenMP (This is also adopted from libint test files with some modifications)
// In the main file we can also use other operators such as kinetic, nuclear, etc but since we only need overlap operator here we omit those parts.
MATRIX compute_1body_ints_parallel(const std::vector<libint2::Shell>& shells_1, const std::vector<libint2::Shell>& shells_2,
libint2::Operator obtype, int nthreads){
  //using libint2::nthreads;

  const auto n_1 = nbasis(shells_1);
  const auto n_2 = nbasis(shells_2);

  // The AO overlap matrix n x n
  Matrix result(n_1,n_2);
  int max_n;
  int max_lval;
  if (max_nprim(shells_1)>max_nprim(shells_2))
  { 
    max_n = max_nprim(shells_1);
  }
  else
  { 
    max_n = max_nprim(shells_2);
  }
  if (max_l(shells_1)>max_l(shells_2))
  { 
    max_lval = max_l(shells_1);
  }
  else
  {
    max_lval = max_l(shells_2);
  }

  // Initializing different engines based on the nthreads
  std::vector<libint2::Engine> engines(nthreads);
  // Initializing the first engine. Others will be the same as this one (as is shown below)
  engines[0] = libint2::Engine(obtype, max_n, max_lval, 0);
  // This part is for other operators, I'll keep it but not necessary.
  //  engines[0].set_params(oparams);
  // Other engines
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  // Mapping the shells into basis sets
  auto shell2bf1 = map_shell_to_basis_function(shells_1);
  auto shell2bf2 = map_shell_to_basis_function(shells_2);
  // Compute for each thread_id
  auto compute = [&](int thread_id) {

    const auto& buf = engines[thread_id].results();

    for (auto s1 = 0l, s12 = 0l; s1 != shells_1.size(); ++s1) {
//    for (auto s1 = 0; s1 < shells_1.size(); ++s1) {
      auto bf1 = shell2bf1[s1];     // first basis function in this shell
      auto n1 = shells_1[s1].size();
      //std::cout << "shells_1[" << s1 << "]" << shells_1[s1].size() << "\n";   
      auto s1_offset = s1 * (s1+1) / 2;
      //for (auto s2 = 0; s2 < shells_2.size(); ++s2) {
      for (auto s2=0; s2!= shells_2.size(); ++s2) {
        //std::cout << "shells_2[" << s2 << "]" << shells_2[s2].size() << "\n";
        auto s12 = s1_offset + s2;
        if (s12 % nthreads != thread_id) continue;
          auto bf2 = shell2bf2[s2];   // first basis function in this shell
          auto n2  = shells_2[s2].size();

        auto n12 = n1 * n2;
		// Make compute for each engine
        engines[thread_id].compute(shells_1[s1], shells_2[s2]);
        // The results block 
        Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
        result.block(bf1, bf2, n1, n2) = buf_mat;
        //if (s1 != s2)  // if s1 >= s2, copy {s1,s2} to the corresponding
        //               // {s2,s1} block, note the transpose!
        //result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
        }
      }
    };
  //std::cout << "Done";
  MATRIX res(n_1,n_2);
  // Now for compute, do parallel
  parallel_do(compute, nthreads);
  for(int i=0;i<n_1;i++){
    for(int j=0; j<n_2;j++){
      res.set(i,j, result(i,j));
    }
  }

  // Return the AO overlap values
  return res;                         
}



// The main function for computing the AO overlaps. We will use this in python.
MATRIX compute_overlaps(const std::vector<libint2::Shell>& shells_1, const std::vector<libint2::Shell>& shells_2, int number_of_threads) {

   //int nthreads;


  //  using libint2::Operator;    
  //  using libint2::nthreads;

    auto nthreads_cstr = getenv("LIBINT_NUM_THREADS");
    int nthreads = number_of_threads;
    if (nthreads_cstr && strcmp(nthreads_cstr, "")) {
       std::istringstream iss(nthreads_cstr);
       iss >> nthreads;
       if (nthreads > 1 << 16 || nthreads <= 0) nthreads = 1;
    }

#if defined(_OPENMP)
      omp_set_num_threads(nthreads);
#endif
//     ;; //std::cout << "Will scale over " << nthreads
//#if defined(_OPENMP)
//     ;;//          << " OpenMP"
//#else
//     ;;//          << " C++11"
//#endif
//     ;;//          << " threads" << std::endl;
    // Initialize Libint    
    libint2::initialize();
	// Compute the AO overlap matrix
    auto S = compute_1body_ints_parallel(shells_1, shells_2, Operator::overlap, nthreads);
    //std::cout << "\n\tFinished computing overlap integral\n";
	// End of AO matrix calculation
    libint2::finalize(); // done with libint

    return S;
}


/*
MATRIX compute_overlaps_serial(const std::vector<libint2::Shell>& shells_1, const std::vector<libint2::Shell>& shells_2) {
    // Initialize Libint    
    libint2::initialize();
        // Compute the AO overlap matrix
    auto S = compute_1body_ints(shells_1, shells_2, Operator::overlap);
    //std::cout << "\n\tFinished computing overlap integral\n";
    libint2::finalize(); // done with libint
    return S;
}

*/


}// namespace liblibint2_wrappers
}// namespace liblibra

