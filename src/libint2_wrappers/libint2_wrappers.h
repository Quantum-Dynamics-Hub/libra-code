/*********************************************************************************
* Copyright (C) 2021 Mohammad Shakiba and Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libint2_wrappers.h
  \brief The file describes the classes and functions that interface libint2 functions and types with those of Libra and expose them to Python
    
*/

#ifndef LIBINT2_WRAPPERS_H
#define LIBINT2_WRAPPERS_H


#if defined(USING_PCH)
#include "../pch.h"
#else

// standard C++ headers
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <thread>

// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#endif 


// Libint Gaussian integrals library
#include <libint2.hpp>


#if !LIBINT2_CONSTEXPR_STATICS
#  include <libint2/statics_definition.h>
#endif


// OpenMP library
#if defined(_OPENMP)
#include <omp.h>
#endif


#include "../math_linalg/liblinalg.h"
#include "../math_specialfunctions/libspecialfunctions.h"




/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libspecialfunctions;


namespace liblibint2_wrappers{



// fires off \c nthreads instances of lambda in parallel
template <typename Lambda>
void parallel_do(Lambda& lambda, int nthreads) {
#ifdef _OPENMP
#pragma omp parallel
{
  auto thread_id = omp_get_thread_num();
  lambda(thread_id);
}
#else  // use C++11 threads
std::vector<std::thread> threads;
for (int thread_id = 0; thread_id != nthreads; ++thread_id) {
  if (thread_id != nthreads - 1)
    threads.push_back(std::thread(lambda, thread_id));
  else
    lambda(thread_id);
}  // threads_id
for (int thread_id = 0; thread_id < nthreads - 1; ++thread_id)
  threads[thread_id].join();
#endif
}




std::vector<libint2::Shell> initialize_shell(int l_val, bool is_spherical, 
 const std::vector<double>& exponents, const std::vector<double>& coeff, VECTOR& coords);

void add_to_shell(std::vector<libint2::Shell>& shells, 
  int l_val, bool is_spherical , const std::vector<double>& exponents, const std::vector<double>& coeff, VECTOR& coords);

void print_shells(std::vector<libint2::Shell>& shells);


size_t nbasis(const std::vector<libint2::Shell>& shells);
size_t max_nprim(const std::vector<libint2::Shell>& shells);
int max_l(const std::vector<libint2::Shell>& shells);
std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells);

using real_t = libint2::scalar_type;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>  Matrix;  

//MATRIX compute_1body_ints(const std::vector<libint2::Shell>& shells_1, const std::vector<libint2::Shell>& shells_2,libint2::Operator obtype);
MATRIX compute_1body_ints_parallel(const std::vector<libint2::Shell>& shells_1, const std::vector<libint2::Shell>& shells_2,libint2::Operator obtype);
MATRIX compute_overlaps(const std::vector<libint2::Shell>& shells_1, const std::vector<libint2::Shell>& shells_2, int number_of_threads);
//MATRIX compute_overlaps_serial(const std::vector<libint2::Shell>& shells_1, const std::vector<libint2::Shell>& shells_2);



typedef std::vector< libint2::Shell > libint2_ShellList;  ///< Data type that holds a vector of libint2::Shell objects



}// namespace liblibint2_wrappers
}// namespace liblibra



#endif // LIBINT2_WRAPPERS_H
