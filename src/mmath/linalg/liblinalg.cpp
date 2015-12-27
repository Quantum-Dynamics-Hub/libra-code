/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "liblinalg.h"

using namespace boost::python;

//using namespace libmmath::liblinalg;
//using libmmath::VECTOR;
//using libmmath::MATRIX;
//using libmmath::CMATRIX;
//using libmmath::MATRIX3x3;
//using libmmath::QUATERNION;
//using libmmath::DATA;



namespace libmmath{
namespace liblinalg{


void export_linalg_objects(){


  class_< intList >("intList")
      .def(vector_indexing_suite< intList >())
  ;

  class_< floatList >("floatList")
      .def(vector_indexing_suite< floatList >())
  ;

  class_< doubleList >("doubleList")
      .def(vector_indexing_suite< doubleList >())
  ;

  class_< complexList >("complexList")
      .def(vector_indexing_suite< complexList >())
  ;


  class_< intMap >("intMap")
      .def(vector_indexing_suite< intMap >())
  ;

  class_< floatMap >("floatMap")
      .def(vector_indexing_suite< floatMap >())
  ;

  class_< doubleMap >("doubleMap")
      .def(vector_indexing_suite< doubleMap >())
  ;

  class_< complexMap >("complexMap")
      .def(vector_indexing_suite< complexMap >())
  ;






void    (VECTOR::*cross1)(VECTOR&,VECTOR&) = &VECTOR::cross;
//VECTOR  (*cross2)(const double, const VECTOR& ,const VECTOR&) = &cross;


double (MATRIX::*get1)(int)            = &MATRIX::get;
double (MATRIX::*get2)(int,int)        = &MATRIX::get;
void   (MATRIX::*set1)(int,double)     = &MATRIX::set;
void   (MATRIX::*set2)(int,int,double) = &MATRIX::set; 
void   (MATRIX::*expt_Inverse_v1)(MATRIX*)    = &MATRIX::Inverse;
void   (MATRIX::*expt_Inverse_v2)(MATRIX&)    = &MATRIX::Inverse;
double (MATRIX::*expt_dot_product_v1)(MATRIX*)    = &MATRIX::dot_product;
double (MATRIX::*expt_dot_product_v2)(MATRIX&)    = &MATRIX::dot_product;
void   (MATRIX::*expt_show_matrix_v1)() = &MATRIX::show_matrix;
void   (MATRIX::*expt_show_matrix_v2)(char*) = &MATRIX::show_matrix;
void   (MATRIX::*expt_exp_v1)(MATRIX&) = &MATRIX::exp;
//void   (MATRIX::*expt_exp_v2)(const MATRIX&) = &MATRIX::exp;
void (MATRIX::*expt_JACOBY_EIGEN_v1)(MATRIX&, MATRIX&) = &MATRIX::JACOBY_EIGEN;
void (MATRIX::*expt_JACOBY_EIGEN_v2)(MATRIX&, MATRIX&,double) = &MATRIX::JACOBY_EIGEN;





complex<double> (CMATRIX::*get3)(int)            = &CMATRIX::get;
complex<double> (CMATRIX::*get4)(int,int)        = &CMATRIX::get;
void   (CMATRIX::*set3)(int,double,double)       = &CMATRIX::set;
void   (CMATRIX::*set4)(int,int,double,double)   = &CMATRIX::set; 
void   (CMATRIX::*set5)(int,complex<double>)       = &CMATRIX::set;
void   (CMATRIX::*set6)(int,int,complex<double>)   = &CMATRIX::set; 


void (CMATRIX::*tridiagonalize1)(CMATRIX& T)              = &CMATRIX::tridiagonalize;
void (CMATRIX::*tridiagonalize2)(CMATRIX& T,CMATRIX& H)   = &CMATRIX::tridiagonalize;


//----------- MATHEMATICS.h -------------------------

  class_<VECTOR>("VECTOR",init<>())
      .def(init<double,double,double>())
      .def(init<const VECTOR&>())
      .def("__copy__", &generic__copy__<VECTOR>) 
      .def("__deepcopy__", &generic__deepcopy__<VECTOR>)
      .def_readwrite("x",&VECTOR::x)
      .def_readwrite("y",&VECTOR::y)
      .def_readwrite("z",&VECTOR::z)
      .def_readwrite("is_transposed",&VECTOR::is_transposed)
      .def(self+self)
      .def(self+double())
      .def(self-self)
      .def(self*double())
      .def(double()*self)
      .def(self/double())
      .def(self*self)
//      .def(self*MATRIX())
      .def(MATRIX()*self)
      .def(self==other<VECTOR>())
      .def(self!=other<VECTOR>())

      .def("length",&VECTOR::length)
      .def("length2",&VECTOR::length2)
      .def("unit",&VECTOR::unit)
      .def("normalize",&VECTOR::normalize)
      .def("cross",cross1)
//      .def("cross",cross2)
      
      
  ;

  class_< VECTORList >("VECTORList")
      .def(vector_indexing_suite< VECTORList >())
  ;

  class_< VECTORMap >("VECTORMap")
      .def(vector_indexing_suite< VECTORMap >())
  ;




  class_<MATRIX>("MATRIX",init<>())      
      .def(init<int,int>())
      .def(init<const MATRIX&>())
      .def("__copy__", &generic__copy__<MATRIX>)
      .def("__deepcopy__", &generic__deepcopy__<MATRIX>)

      .def("get",get1)
      .def("get",get2)
      .def("set",set1)
      .def("set",set2)

      .def_readwrite("num_of_cols",&MATRIX::num_of_cols)
      .def_readwrite("num_of_rows",&MATRIX::num_of_rows)
      .def_readwrite("num_of_elems",&MATRIX::num_of_elems)
      .def(self+self)
      .def(self-self)
      .def(self*self)
      .def(self/double())
      .def(self*double())
      .def(double()*self)

      .def("Init",&MATRIX::Init)
      .def("InitSquareMatrix",&MATRIX::InitSquareMatrix)
      .def("Init_Unit_Matrix",&MATRIX::Init_Unit_Matrix)
      .def("Load_Matrix_From_File",&MATRIX::Load_Matrix_From_File)

      .def("Rotation",&MATRIX::Rotation)
      .def("Rx",&MATRIX::Rx)
      .def("Ry",&MATRIX::Ry)
      .def("Rz",&MATRIX::Rz)

      .def("show_num_of_rows",&MATRIX::show_num_of_rows)
      .def("show_num_of_cols",&MATRIX::show_num_of_cols)
      .def("show_num_of_elems",&MATRIX::show_num_of_elems)

      .def("Add_To_Element",&MATRIX::Add_To_Element)
      .def("FindMaxNondiagonalElement",&MATRIX::FindMaxNondiagonalElement)

      .def("RightRotation",&MATRIX::RightRotation)
      .def("LeftRotation",&MATRIX::LeftRotation)
      .def("Ortogonalization",&MATRIX::Ortogonalization)
      .def("Transpose",&MATRIX::Transpose)
      .def("T",&MATRIX::T)
      .def("Inverse",expt_Inverse_v2)    
      .def("tensor_product",&MATRIX::tensor_product)
      .def("dot_product",expt_dot_product_v2)
      .def("col", &MATRIX::col)     // return given column of the matrix 
      .def("row", &MATRIX::row)     // return given row of the matrix 

      .def("show_matrix", expt_show_matrix_v1)
      .def("show_matrix", expt_show_matrix_v2)

      .def("get_vectors", &MATRIX::get_vectors)
      .def("skew",&MATRIX::skew)
      .def("skew1",&MATRIX::skew1)
      .def("exp", expt_exp_v1)
//      .def("exp", expt_exp_v2)
      .def("JACOBY_EIGEN", expt_JACOBY_EIGEN_v1)
      .def("JACOBY_EIGEN", expt_JACOBY_EIGEN_v2)

      .def_readwrite("MATRIX_PRECISION",&MATRIX::MATRIX_PRECISION)
      .def_readwrite("MATRIX_WIDTH",&MATRIX::MATRIX_WIDTH)

      .def("Determinant",&MATRIX::Determinant)
      .def("tr",&MATRIX::tr)
      .def("max_elt",&MATRIX::max_elt)

      .def("bin_dump",&MATRIX::bin_dump)
      .def("bin_load",&MATRIX::bin_load)

  ;



  class_< MATRIXList >("MATRIXList")
      .def(vector_indexing_suite< MATRIXList >())
  ;

  class_< MATRIXMap >("MATRIXMap")
      .def(vector_indexing_suite< MATRIXMap >())
  ;



  class_<MATRIX3x3>("MATRIX3x3",init<>())      
      .def(init<const VECTOR&, const VECTOR&, const VECTOR&>())
      .def(init<const MATRIX3x3&>())
      .def("__copy__", &generic__copy__<MATRIX3x3>)
      .def("__deepcopy__", &generic__deepcopy__<MATRIX3x3>)


      .def("get",get1)
      .def("get",get2)
      .def("set",set1)
      .def("set",set2)

      .def_readwrite("xx",&MATRIX3x3::xx)
      .def_readwrite("xy",&MATRIX3x3::xy)
      .def_readwrite("xz",&MATRIX3x3::xz)
      .def_readwrite("yx",&MATRIX3x3::yx)
      .def_readwrite("yy",&MATRIX3x3::yy)
      .def_readwrite("yz",&MATRIX3x3::yz)
      .def_readwrite("zx",&MATRIX3x3::zx)
      .def_readwrite("zy",&MATRIX3x3::zy)
      .def_readwrite("zz",&MATRIX3x3::zz)

      .def(self+self)
      .def(self-self)
      .def(self*self)
      .def(self/double())
      .def(self*double())
      .def(double()*self)

  ;

  class_< MATRIX3x3List >("MATRIX3x3List")
      .def(vector_indexing_suite< MATRIX3x3List >())
  ;

  class_< MATRIX3x3Map >("MATRIX3x3Map")
      .def(vector_indexing_suite< MATRIX3x3Map >())
  ;



  class_<CMATRIX>("CMATRIX",init<>())      
      .def(init<int,int>())
      .def(init<MATRIX&>())
      .def(init<MATRIX&,MATRIX&>())
      .def(init<const CMATRIX&>())
      .def("__copy__", &generic__copy__<CMATRIX>)
      .def("__deepcopy__", &generic__deepcopy__<CMATRIX>)

      .def("get",get3)
      .def("get",get4)
      .def("set",set3)
      .def("set",set4)
      .def("set",set5)
      .def("set",set6)


      .def_readwrite("num_of_cols",&CMATRIX::n_cols)
      .def_readwrite("num_of_rows",&CMATRIX::n_rows)
      .def_readwrite("num_of_elems",&CMATRIX::n_elts)
      .def(self+self)
      .def(self-self)
      .def(self*self)
      .def(self/double())
      .def(self*double())
      .def(double()*self)

      .def("show_matrix", &CMATRIX::show_matrix)
      .def("conj", &CMATRIX::conj)   // return complex conjugate matrix
      .def("T", &CMATRIX::T)         // return transposed matrix
      .def("H", &CMATRIX::H)         // return Hermitian-conjugate matrix
      .def("load_identity", &CMATRIX::load_identity)
      .def("dot", &CMATRIX::dot)     // dot product of two matrices
      .def("col", &CMATRIX::col)     // return given column of the matrix 
      .def("row", &CMATRIX::row)     // return given row of the matrix 

      .def("QR", &CMATRIX::QR)       // QR for general (Hermitian or symmetric) matrices
      .def("QR1",&CMATRIX::QR1)      // QR for tridiagonal matrices

      .def("tridiagonalize", tridiagonalize1)  // for Hermitian or symmetric CMATRIX - only resulting tridiagonal CMATRIX
      .def("tridiagonalize", tridiagonalize2)  // ---//---  also keep track of Householder transformation matrices

      .def("eigen" , &CMATRIX::eigen )   // interface
      .def("eigen0", &CMATRIX::eigen0)   // Schur decomposition or Jacobi rotation
      .def("eigen1", &CMATRIX::eigen1)   // only eigenvalues - fast
      .def("eigen2", &CMATRIX::eigen2)   // also eigenvectors, slower - keep track of transformation matrixes
      .def("eigen3", &CMATRIX::eigen3)   // also eigenvectors - solve for each eigenvector independently

      .def("inverse",        &CMATRIX::inverse)        // interface
      .def("direct_inverse", &CMATRIX::direct_inverse)

      .def("bin_dump",&CMATRIX::bin_dump)
      .def("bin_load",&CMATRIX::bin_load)


/*
      .def("Init",&CMATRIX::Init)
      .def("InitSquareMatrix",&CMATRIX::InitSquareMatrix)
      .def("Init_Unit_Matrix",&CMATRIX::Init_Unit_Matrix)
      .def("Load_Matrix_From_File",&CMATRIX::Load_Matrix_From_File)
      .def("show_num_of_rows",&CMATRIX::show_num_of_rows)
      .def("show_num_of_cols",&CMATRIX::show_num_of_cols)
      .def("show_num_of_elems",&CMATRIX::show_num_of_elems)
      .def("RightRotation",&CMATRIX::RightRotation)
      .def("LeftRotation",&CMATRIX::LeftRotation)
      .def("Ortogonalization",&CMATRIX::Ortogonalization)
//      .def("Inverse",Inverse1)  
//      .def("Inverse",Inverse2)    
*/
  ;

  class_< CMATRIXList >("CMATRIXList")
      .def(vector_indexing_suite< CMATRIXList >())
  ;

  class_< CMATRIXMap >("CMATRIXMap")
      .def(vector_indexing_suite< CMATRIXMap >())
  ;





  class_<QUATERNION>("QUATERNION",init<>())
      .def(init<double,double,double,double>())
      .def("__copy__", &generic__copy__<QUATERNION>)
      .def("__deepcopy__", &generic__deepcopy__<QUATERNION>)

      .def_readwrite("Lt",&QUATERNION::Lt)
      .def_readwrite("Lx",&QUATERNION::Lx)
      .def_readwrite("Ly",&QUATERNION::Ly)
      .def_readwrite("Lz",&QUATERNION::Lz)
      .def(self+self)
      .def(self-self)
      .def(self*self)
      .def(self*double())
      .def(double()*self)
      .def(self+=self)
      .def(self-=self)

  ;

  class_< QUATERNIONList >("QUATERNIONList")
      .def(vector_indexing_suite< QUATERNIONList >())
  ;

  class_< QUATERNIONMap >("QUATERNIONMap")
      .def(vector_indexing_suite< QUATERNIONMap >())
  ;



}// export_linalg_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyglinalg){
#else
BOOST_PYTHON_MODULE(liblinalg){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_linalg_objects();

}





}// namespace liblinalg
}// namespace libmmath




