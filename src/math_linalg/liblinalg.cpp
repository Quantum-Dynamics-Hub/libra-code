/*********************************************************************************
* Copyright (C) 2017-2020 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file liblinalg.cpp
  \brief The file implements Python export function and data types
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
#endif 

#include "liblinalg.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;


/// liblinalg namespace
namespace liblinalg{



template <typename T1>
void export_base_matrix(){
  
  void (base_matrix<T1>::*expt_set_v1)(int i, T1 val) = &base_matrix<T1>::set;
  void (base_matrix<T1>::*expt_set_v2)(int i, int j, T1 val) = &base_matrix<T1>::set;

  T1 (base_matrix<T1>::*expt_get_v1)(int i) = &base_matrix<T1>::get;
  T1 (base_matrix<T1>::*expt_get_v2)(int i, int j) = &base_matrix<T1>::get;

  void (base_matrix<T1>::*expt_diag_v1)(int dim, T1 x) = &base_matrix<T1>::diag;
  void (base_matrix<T1>::*expt_diag_v2)(T1 x) = &base_matrix<T1>::diag;
  void (base_matrix<T1>::*expt_identity_v1)() = &base_matrix<T1>::identity;

  void (base_matrix<T1>::*expt_add_v1)(int row,int col, T1 x) = &base_matrix<T1>::add;
//  void (base_matrix<T1>::*expt_add_v1)(int row,int col, int x) = &base_matrix<T1>::add;
//  void (base_matrix<T1>::*expt_add_v2)(int row,int col, double x) = &base_matrix<T1>::add;
  void (base_matrix<T1>::*expt_scale_v1)(int row,int col, T1 x) = &base_matrix<T1>::scale;
//  void (base_matrix<T1>::*expt_scale_v1)(int row,int col, int x) = &base_matrix<T1>::scale;
//  void (base_matrix<T1>::*expt_scale_v2)(int row,int col, double x) = &base_matrix<T1>::scale;


  T1 (base_matrix<T1>::*expt_sum_col_v1)(int icol) = &base_matrix<T1>::sum_col;
  T1 (base_matrix<T1>::*expt_sum_col_v2)(int icol, int power) = &base_matrix<T1>::sum_col;
  T1 (base_matrix<T1>::*expt_sum_row_v1)(int irow) = &base_matrix<T1>::sum_row;
  T1 (base_matrix<T1>::*expt_sum_row_v2)(int irow, int power) = &base_matrix<T1>::sum_row;

  T1 (base_matrix<T1>::*expt_prod_col_v1)(int icol) = &base_matrix<T1>::prod_col;
  T1 (base_matrix<T1>::*expt_prod_row_v1)(int irow) = &base_matrix<T1>::prod_row;



  void (base_matrix<T1>::*expt_show_matrix_v1)() = &base_matrix<T1>::show_matrix;
  void (base_matrix<T1>::*expt_show_matrix_v2)(char * Output_File) = &base_matrix<T1>::show_matrix;


  // Also expose some functions
  void (*expt_pop_submatrix_v1)(base_matrix<T1>& X, base_matrix<T1>& x, vector<int>& subset) = &pop_submatrix;
  void (*expt_pop_submatrix_v2)(base_matrix<T1>& X, base_matrix<T1>& x, boost::python::list subset) = &pop_submatrix;
  void (*expt_pop_submatrix_v3)(base_matrix<T1>& X, base_matrix<T1>& x, vector<int>& subset,vector<int>& subset2) = &pop_submatrix;
  void (*expt_pop_submatrix_v4)(base_matrix<T1>& X, base_matrix<T1>& x, boost::python::list subset,boost::python::list subset2) = &pop_submatrix;

  void (*expt_push_submatrix_v1)(base_matrix<T1>& X, base_matrix<T1>& x, vector<int>& subset) = &push_submatrix;
  void (*expt_push_submatrix_v2)(base_matrix<T1>& X, base_matrix<T1>& x, boost::python::list subset) = &push_submatrix;
  void (*expt_push_submatrix_v3)(base_matrix<T1>& X, base_matrix<T1>& x, vector<int>& subset,vector<int>& subset2) = &push_submatrix;
  void (*expt_push_submatrix_v4)(base_matrix<T1>& X, base_matrix<T1>& x, boost::python::list subset,boost::python::list subset2) = &push_submatrix;

  void (*expt_add_submatrix_v1)(base_matrix<T1>& X, base_matrix<T1>& x, vector<int>& subset, T1 alpha) = &add_submatrix;
  void (*expt_add_submatrix_v2)(base_matrix<T1>& X, base_matrix<T1>& x, boost::python::list subset, T1 alpha) = &add_submatrix;
  void (*expt_add_submatrix_v3)(base_matrix<T1>& X, base_matrix<T1>& x, vector<int>& subset,vector<int>& subset2, T1 alpha) = &add_submatrix;
  void (*expt_add_submatrix_v4)(base_matrix<T1>& X, base_matrix<T1>& x, boost::python::list subset,boost::python::list subset2, T1 alpha) = &add_submatrix;


  class_<   base_matrix<T1> >("base_matrix_general",init<>())
      .def(init<int,int>())
      .def(init<const base_matrix<T1>& >())
      .def("__copy__", &generic__copy__<base_matrix<T1> >) 
      .def("__deepcopy__", &generic__deepcopy__<base_matrix<T1> >)

      /// Direct access to matrix basic properties
      .def_readwrite("num_of_cols",&base_matrix<T1>::n_cols)
      .def_readwrite("num_of_rows",&base_matrix<T1>::n_rows)
      .def_readwrite("num_of_elems",&base_matrix<T1>::n_elts)

      /// Getters and setters    
      .def("set", expt_set_v1)
      .def("set", expt_set_v2)
      .def("get", expt_get_v1)
      .def("get", expt_get_v2)

//      .def("col", &base_matrix<T1>::col )
//      .def("row", &base_matrix<T1>::row )

      /// Generic initializaitons 
      .def("diag", expt_diag_v1)
      .def("diag", expt_diag_v2)
      .def("identity", expt_identity_v1)
      .def("Init", &base_matrix<T1>::Init )
      .def("InitSquareMatrix", &base_matrix<T1>::InitSquareMatrix )
      .def("Init_Unit_Matrix", &base_matrix<T1>::Init_Unit_Matrix )

      /// Generic operations
      .def("add", expt_add_v1)
//      .def("add", expt_add_v2)
      .def("scale", expt_scale_v1)
//      .def("scale", expt_scale_v2)
      .def("product", &base_matrix<T1>::product )
      .def("dot_product", &base_matrix<T1>::dot_product )

      /// Generic matrix modifiers
      .def("Transpose", &base_matrix<T1>::Transpose )
//      .def("T", &base_matrix<T1>::T )
      .def("swap_cols", &base_matrix<T1>::swap_cols )
      .def("swap_rows", &base_matrix<T1>::swap_rows )
      .def("permute_cols", &base_matrix<T1>::permute_cols )
      .def("permute_rows", &base_matrix<T1>::permute_rows )
      .def("RightRotation", &base_matrix<T1>::RightRotation )
      .def("LeftRotation", &base_matrix<T1>::LeftRotation )

      /// Generic properties
      .def("tr", &base_matrix<T1>::tr )
      .def("sum", &base_matrix<T1>::sum )

      .def("sum_col", expt_sum_col_v1)
      .def("sum_col", expt_sum_col_v2)
      .def("sum_row", expt_sum_row_v1)
      .def("sum_row", expt_sum_row_v2)
      .def("prod_col", expt_prod_col_v1)
      .def("prod_row", expt_prod_row_v1)


      /// Generic IO operations
      .def("bin_dump", &base_matrix<T1>::bin_dump )
      .def("bin_load", &base_matrix<T1>::bin_load )
      .def("show_matrix", expt_show_matrix_v1)
      .def("show_matrix", expt_show_matrix_v2)
      .def("show_matrix_address", &base_matrix<T1>::show_matrix_address)
      .def("Load_Matrix_From_File", &base_matrix<T1>::Load_Matrix_From_File)

      /// Generic operator overloads
      .def(self+=self)
      .def(self+=int())
      .def(self+=double())

      .def(self-=self)
      .def(self-=int())
      .def(self-=double())

      .def(self*=int())
      .def(self*=double())


      .def(self/=int())
      .def(self/=double())



  ;

  def("pop_submatrix", expt_pop_submatrix_v1);
  def("pop_submatrix", expt_pop_submatrix_v2);
  def("pop_submatrix", expt_pop_submatrix_v3);
  def("pop_submatrix", expt_pop_submatrix_v4);

  def("push_submatrix", expt_push_submatrix_v1);
  def("push_submatrix", expt_push_submatrix_v2);
  def("push_submatrix", expt_push_submatrix_v3);
  def("push_submatrix", expt_push_submatrix_v4);

  def("add_submatrix", expt_add_submatrix_v1);
  def("add_submatrix", expt_add_submatrix_v2);
  def("add_submatrix", expt_add_submatrix_v3);
  def("add_submatrix", expt_add_submatrix_v4);


}



void export_IMATRIX(){

  export_base_matrix< int >();

  void (IMATRIX::*expt_max_col_elt_v1)(int, int&, int&) = &IMATRIX::max_col_elt;
  boost::python::list (IMATRIX::*expt_max_col_elt_v2)(int) = &IMATRIX::max_col_elt;

  void (IMATRIX::*expt_min_col_elt_v1)(int, int&, int&) = &IMATRIX::min_col_elt;
  boost::python::list (IMATRIX::*expt_min_col_elt_v2)(int) = &IMATRIX::min_col_elt;

  void (IMATRIX::*expt_max_row_elt_v1)(int, int&, int&) = &IMATRIX::max_row_elt;
  boost::python::list (IMATRIX::*expt_max_row_elt_v2)(int) = &IMATRIX::max_row_elt;

  void (IMATRIX::*expt_min_row_elt_v1)(int, int&, int&) = &IMATRIX::min_row_elt;
  boost::python::list (IMATRIX::*expt_min_row_elt_v2)(int) = &IMATRIX::min_row_elt;

//  void (IMATRIX::*expt_set_v1)(int, int) = &IMATRIX::set;
//  void (IMATRIX::*expt_set_v2)(int, int, int) = &IMATRIX::set;
/*
  void (IMATRIX::*expt_set_v2)(int, double) = &IMATRIX::set;
  void (IMATRIX::*expt_set_v4)(int, int, double) = &IMATRIX::set;
*/

  class_< IMATRIX, boost::python::bases<base_matrix< int > > >("IMATRIX",init<>())
      .def(init<int,int>())
      .def(init<const IMATRIX&>())
//      .def(init<const base_matrix<complex<double> >&>())
      .def("__copy__", &generic__copy__<IMATRIX>) 
      .def("__deepcopy__", &generic__deepcopy__<IMATRIX>)

      /// Setters and getters
//      .def("set", expt_set_v1)
//      .def("set", expt_set_v2)
/*
      .def("set", expt_set_v3)
      .def("set", expt_set_v4)
*/



      /// Initializations
//      .def("diag", &IMATRIX::diag)
//      .def("identity", &IMATRIX::identity)
      .def("load_identity", &IMATRIX::load_identity)

      .def("Init",&IMATRIX::Init)
      .def("InitSquareMatrix",&IMATRIX::InitSquareMatrix)
      .def("Init_Unit_Matrix",&IMATRIX::Init_Unit_Matrix)
      .def("Load_Matrix_From_File",&IMATRIX::Load_Matrix_From_File)


      /// Returning matrix derivatives
      .def("T", &IMATRIX::T)         // return a transposed matrix
      .def("col", &IMATRIX::col)     // return given column of the matrix 
      .def("row", &IMATRIX::row)     // return given row of the matrix 

      /// Properties of the matrix
      .def("max_elt", &IMATRIX::max_elt)
      .def("FindMaxNondiagonalElement", &IMATRIX::FindMaxNondiagonalElement)
      .def("max_nondiagonal", &IMATRIX::max_nondiagonal)
      .def("max_col_elt", expt_max_col_elt_v1)
      .def("max_col_elt", expt_max_col_elt_v2)
      .def("min_col_elt", expt_min_col_elt_v1)
      .def("min_col_elt", expt_min_col_elt_v2)
      .def("max_row_elt", expt_max_row_elt_v1)
      .def("max_row_elt", expt_max_row_elt_v2)
      .def("min_row_elt", expt_min_row_elt_v1)
      .def("min_row_elt", expt_min_row_elt_v2)


       /// Type-specific (returning this datatype) operator overloads 
      .def(self+self)
      .def(self+int())
      .def(self-self)
      .def(self-int())
      .def(self*self)
      .def(self*int())
      .def(int()*self)
      .def(self/int())

  ;

  class_< IMATRIXList >("IMATRIXList")
      .def(vector_indexing_suite< IMATRIXList >())
  ;

  class_< IMATRIXMap >("IMATRIXMap")
      .def(vector_indexing_suite< IMATRIXMap >())
  ;


}




struct MATRIX_pickle_suite : boost::python::pickle_suite
{
/** 
  Figured out it thanks to these sources:

  https://stackoverflow.com/questions/37180878/append-element-to-boostpythontuple

  https://valelab4.ucsf.edu/svn/3rdpartypublic/boost/libs/python/doc/v2/pickle.html


*/

  static boost::python::tuple getinitargs(const MATRIX& w)
  {
    boost::python::tuple res;
    res = boost::python::make_tuple(w.n_rows, w.n_cols);
      
    return res;
  }


  static boost::python::tuple getstate(const MATRIX& w)
  {
    boost::python::tuple res;

    for(int i=0; i<w.n_elts; i++){
      res = boost::python::extract< boost::python::tuple >( res +  boost::python::make_tuple( w.get(i) ) );
    }// for i

    return res;
  }

  static void setstate(MATRIX& w, boost::python::tuple state)
  {

    int sz = len(state);
    for(int i=0; i<sz; i++){  
      w.set(i, extract<double>(state[i]));
    }
  }

};





void export_MATRIX(){

  // It is important to wrap the base class - or the derived class wrapping won't work
  // unfortunately, it seems like we need to do this for every data type
  export_base_matrix<double>();


  void (MATRIX::*expt_init_v1)(const VECTOR& u1, const VECTOR& u2, const VECTOR& u3) = &MATRIX::init;

  void (MATRIX::*expt_max_col_elt_v1)(int I, double& val, int& max_elt_indx) = &MATRIX::max_col_elt;
  boost::python::list (MATRIX::*expt_max_col_elt_v2)(int I) = &MATRIX::max_col_elt;
  void (MATRIX::*expt_min_col_elt_v1)(int I, double& val, int& max_elt_indx) = &MATRIX::min_col_elt;
  boost::python::list (MATRIX::*expt_min_col_elt_v2)(int I) = &MATRIX::min_col_elt;

  void (MATRIX::*expt_max_row_elt_v1)(int I, double& val, int& max_elt_indx) = &MATRIX::max_row_elt;
  boost::python::list (MATRIX::*expt_max_row_elt_v2)(int I) = &MATRIX::max_row_elt;
  void (MATRIX::*expt_min_row_elt_v1)(int I, double& val, int& max_elt_indx) = &MATRIX::min_row_elt;
  boost::python::list (MATRIX::*expt_min_row_elt_v2)(int I) = &MATRIX::min_row_elt;


  class_< MATRIX, boost::python::bases<base_matrix<double> > >("MATRIX",init<>())
      .def(init<int,int>())
//      .def(init<const base_matrix<double>&>())
      .def(init<const MATRIX&>())
//      .def(init<boost::python::tuple>())
      .def(init<VECTOR,VECTOR,VECTOR>())
      .def("__copy__", &generic__copy__<MATRIX>) 
      .def("__deepcopy__", &generic__deepcopy__<MATRIX>)

      /// Initializations
      .def("Init",&MATRIX::Init)
      .def("InitSquareMatrix",&MATRIX::InitSquareMatrix)
      .def("Init_Unit_Matrix",&MATRIX::Init_Unit_Matrix)
      .def("Load_Matrix_From_File",&MATRIX::Load_Matrix_From_File)
      .def("init", expt_init_v1)

      /// Returning matrix derivatives
      .def("T", &MATRIX::T )
      .def("col", &MATRIX::col)
      .def("row", &MATRIX::row)

      /// Properties of the matrix
      .def("NonOrtogonality_Measure", &MATRIX::NonOrtogonality_Measure)
      .def("max_elt", &MATRIX::max_elt)
      .def("FindMaxNondiagonalElement", &MATRIX::FindMaxNondiagonalElement)
      .def("max_col_elt", expt_max_col_elt_v1)
      .def("max_col_elt", expt_max_col_elt_v2)
      .def("min_col_elt", expt_min_col_elt_v1)
      .def("min_col_elt", expt_min_col_elt_v2)
      .def("max_row_elt", expt_max_row_elt_v1)
      .def("max_row_elt", expt_max_row_elt_v2)
      .def("min_row_elt", expt_min_row_elt_v1)
      .def("min_row_elt", expt_min_row_elt_v2)


       /// Type-specific (returning this datatype) operator overloads 
      .def(self+self)
      .def(self+int())
      .def(self+double())

      .def(self-self)
      .def(self-int())
      .def(self-double())

      .def(self*self)
      .def(self*int())
      .def(self*double())
      .def(int()*self)
      .def(double()*self)

      .def(self/int())
      .def(self/double())


      /// Misc methods
      .def("tensor_product", &MATRIX::tensor_product)
      .def("get_vectors", &MATRIX::get_vectors)
      .def("skew", &MATRIX::skew)
      .def("skew1", &MATRIX::skew1)
      .def("Rotation",&MATRIX::Rotation)
      .def("Rx",&MATRIX::Rx)
      .def("Ry",&MATRIX::Ry)
      .def("Rz",&MATRIX::Rz)


      .def_pickle(MATRIX_pickle_suite())
      //.enable_pickling()

  ;


struct MATRIXList_pickle_suite : boost::python::pickle_suite
{
/** 
  Figured out it thanks to these sources:

  https://stackoverflow.com/questions/37180878/append-element-to-boostpythontuple

  https://valelab4.ucsf.edu/svn/3rdpartypublic/boost/libs/python/doc/v2/pickle.html


*/
/*
  static boost::python::tuple getinitargs(const MATRIX& w)
  {
    boost::python::tuple res;
    res = boost::python::make_tuple(w.n_rows, w.n_cols);
      
    return res;
  }
*/

  static boost::python::tuple getstate(const MATRIXList& w)
  {
    boost::python::tuple res;

    int sz = w.size();

    for(int i=0; i<sz; i++){
      res = boost::python::extract< boost::python::tuple >( res +  boost::python::make_tuple( w[i] ) );
    }// for i

    return res;
  }

  static void setstate(MATRIXList& w, boost::python::tuple state)
  {

    int sz = len(state);
    for(int i=0; i<sz; i++){  
      w.push_back( extract<MATRIX>(state[i]) );
    }
  }

};


  class_< MATRIXList >("MATRIXList")
      .def(vector_indexing_suite< MATRIXList >())
      .def_pickle(MATRIXList_pickle_suite())
//      .enable_pickling()
  ;

  class_< MATRIXMap >("MATRIXMap")
      .def(vector_indexing_suite< MATRIXMap >())
  ;


}



struct CMATRIX_pickle_suite : boost::python::pickle_suite
{
/** 
  Figured out it thanks to these sources:

  https://stackoverflow.com/questions/37180878/append-element-to-boostpythontuple

  https://valelab4.ucsf.edu/svn/3rdpartypublic/boost/libs/python/doc/v2/pickle.html


*/

  static boost::python::tuple getinitargs(const CMATRIX& w)
  {
    boost::python::tuple res;
    res = boost::python::make_tuple(w.n_rows, w.n_cols);
      
    return res;
  }


  static boost::python::tuple getstate(const CMATRIX& w)
  {
    boost::python::tuple res;

    for(int i=0; i<w.n_elts; i++){
      res = boost::python::extract< boost::python::tuple >( res +  boost::python::make_tuple( w.get(i) ) );
    }// for i

    return res;
  }

  static void setstate(CMATRIX& w, boost::python::tuple state)
  {

    int sz = len(state);
    for(int i=0; i<sz; i++){  
      w.set(i, extract< complex<double> >(state[i]));
    }
  }

};




void export_CMATRIX(){

  export_base_matrix< complex<double> >();

  void (CMATRIX::*expt_cmax_col_elt_v1)(int, complex<double>&, int&) = &CMATRIX::max_col_elt;
  boost::python::list (CMATRIX::*expt_cmax_col_elt_v2)(int) = &CMATRIX::max_col_elt;

  void (CMATRIX::*expt_cmin_col_elt_v1)(int, complex<double>&, int&) = &CMATRIX::min_col_elt;
  boost::python::list (CMATRIX::*expt_cmin_col_elt_v2)(int) = &CMATRIX::min_col_elt;

  void (CMATRIX::*expt_cmax_row_elt_v1)(int, complex<double>&, int&) = &CMATRIX::max_row_elt;
  boost::python::list (CMATRIX::*expt_cmax_row_elt_v2)(int) = &CMATRIX::max_row_elt;

  void (CMATRIX::*expt_cmin_row_elt_v1)(int, complex<double>&, int&) = &CMATRIX::min_row_elt;
  boost::python::list (CMATRIX::*expt_cmin_row_elt_v2)(int) = &CMATRIX::min_row_elt;

  void (CMATRIX::*expt_cset_v1)(int, double, double) = &CMATRIX::set;
  void (CMATRIX::*expt_cset_v2)(int, int, double, double) = &CMATRIX::set;
  void (CMATRIX::*expt_cset_v3)(int, complex<double>) = &CMATRIX::set;
  void (CMATRIX::*expt_cset_v4)(int, int, complex<double>) = &CMATRIX::set;



  class_< CMATRIX, boost::python::bases<base_matrix<complex<double> > > >("CMATRIX",init<>())
      .def(init<int,int>())
      .def(init<MATRIX&>())
      .def(init<MATRIX&,MATRIX&>())
      .def(init<const CMATRIX&>())
//      .def(init<const base_matrix<complex<double> >&>())
      .def("__copy__", &generic__copy__<CMATRIX>) 
      .def("__deepcopy__", &generic__deepcopy__<CMATRIX>)

      /// Setters and getters
      .def("set", expt_cset_v1)
      .def("set", expt_cset_v2)
      .def("set", expt_cset_v3)
      .def("set", expt_cset_v4)


      /// Initializations
      .def("load_identity", &CMATRIX::load_identity)

      .def("Init",&CMATRIX::Init)
      .def("InitSquareMatrix",&CMATRIX::InitSquareMatrix)
      .def("Init_Unit_Matrix",&CMATRIX::Init_Unit_Matrix)
      .def("Load_Matrix_From_File",&CMATRIX::Load_Matrix_From_File)


      /// Returning matrix derivatives
      .def("T", &CMATRIX::T)         // return a transposed matrix
      .def("H", &CMATRIX::H)         // return a Hermitian-conjugate matrix
      .def("conj", &CMATRIX::conj)   // return a complex conjugate matrix
      .def("real", &CMATRIX::real)
      .def("imag", &CMATRIX::imag)
      .def("get_components", &CMATRIX::get_components)
      .def("col", &CMATRIX::col)     // return given column of the matrix 
      .def("row", &CMATRIX::row)     // return given row of the matrix 

      /// Properties of the matrix
      .def("NonOrtogonality_Measure", &CMATRIX::NonOrtogonality_Measure)
      .def("max_elt", &CMATRIX::max_elt)
      .def("FindMaxNondiagonalElement", &CMATRIX::FindMaxNondiagonalElement)
      .def("max_nondiagonal", &CMATRIX::max_nondiagonal)
      .def("max_col_elt", expt_cmax_col_elt_v1)
      .def("max_col_elt", expt_cmax_col_elt_v2)
      .def("min_col_elt", expt_cmin_col_elt_v1)
      .def("min_col_elt", expt_cmin_col_elt_v2)
      .def("max_row_elt", expt_cmax_row_elt_v1)
      .def("max_row_elt", expt_cmax_row_elt_v2)
      .def("min_row_elt", expt_cmin_row_elt_v1)
      .def("min_row_elt", expt_cmin_row_elt_v2)


       /// Type-specific (returning this datatype) operator overloads 
      .def(self+=complex<double>())
      .def(self-=complex<double>())

      .def(self+self)
      .def(self+int())
      .def(self+double())
      .def(self+complex<double>())

      .def(self-self)
      .def(self-int())
      .def(self-double())
      .def(self-complex<double>())

      .def(self*=complex<double>())

      .def(self*self)
      .def(self*int())
      .def(self*double())
      .def(self*complex<double>())    
/*
      .def(MATRIX()*self)
      .def(self*MATRIX())
      .def(IMATRIX()*self)
      .def(self*IMATRIX())
*/

      .def(int()*self)
      .def(double()*self)
      .def(complex<double>()*self)

      .def(self/=complex<double>())

      .def(self/int())
      .def(self/double())
      .def(self/complex<double>())

      .def_pickle(CMATRIX_pickle_suite())
  ;



struct CMATRIXList_pickle_suite : boost::python::pickle_suite
{
/** 
  Figured out it thanks to these sources:

  https://stackoverflow.com/questions/37180878/append-element-to-boostpythontuple

  https://valelab4.ucsf.edu/svn/3rdpartypublic/boost/libs/python/doc/v2/pickle.html


*/
/*
  static boost::python::tuple getinitargs(const MATRIX& w)
  {
    boost::python::tuple res;
    res = boost::python::make_tuple(w.n_rows, w.n_cols);
      
    return res;
  }
*/

  static boost::python::tuple getstate(const CMATRIXList& w)
  {
    boost::python::tuple res;

    int sz = w.size();

    for(int i=0; i<sz; i++){
      res = boost::python::extract< boost::python::tuple >( res +  boost::python::make_tuple( w[i] ) );
    }// for i

    return res;
  }

  static void setstate(CMATRIXList& w, boost::python::tuple state)
  {

    int sz = len(state);
    for(int i=0; i<sz; i++){  
      w.push_back( extract<CMATRIX>(state[i]) );
    }
  }

};



  class_< CMATRIXList >("CMATRIXList")
      .def(vector_indexing_suite< CMATRIXList >())
      .def_pickle(CMATRIXList_pickle_suite())
  ;

  class_< CMATRIXMap >("CMATRIXMap")
      .def(vector_indexing_suite< CMATRIXMap >())
  ;


}


void export_VECTOR(){

  void (VECTOR::*expt_cross_v1)(const VECTOR &v1, const VECTOR &v2) = &VECTOR::cross;  

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
      .def(self+=double())
      .def(self+=self)

      .def(self-self)
      .def(self-double())
      .def(self-=double())
      .def(self-=self)

      .def(self*self)
      .def(self*double())
      .def(double()*self)
      .def(self*=double())
      .def(MATRIX()*self)

      .def(self/double())
      .def(self/=double())

      .def(self==other<VECTOR>())
      .def(self!=other<VECTOR>())

      .def("length",&VECTOR::length)
      .def("length2",&VECTOR::length2)
      .def("unit",&VECTOR::unit)
      .def("normalize",&VECTOR::normalize)
      .def("cross",expt_cross_v1)
      //.def("cross",cross2)
            
  ;

  VECTOR (*expt_cross_v2)(const double k, const VECTOR &v1, const VECTOR &v2) = &cross;

  def("cross", expt_cross_v2);


  class_< VECTORList >("VECTORList")
      .def(vector_indexing_suite< VECTORList >())
  ;

  class_< VECTORMap >("VECTORMap")
      .def(vector_indexing_suite< VECTORMap >())
  ;


}

void export_MATRIX3x3(){

  void (MATRIX3x3::*expt_eigen_v1)(MATRIX3x3&, MATRIX3x3&) = &MATRIX3x3::eigen;  
  void (MATRIX3x3::*expt_eigen_v2)(MATRIX3x3&, MATRIX3x3&, double) = &MATRIX3x3::eigen;  

  void (MATRIX3x3::*expt_diag_v1)(double x) = &MATRIX3x3::diag;  
  void (MATRIX3x3::*expt_diag_v2)(double x, double y, double z) = &MATRIX3x3::diag;  

  class_<MATRIX3x3>("MATRIX3x3",init<>())      
      .def(init<const VECTOR&, const VECTOR&, const VECTOR&>())
      .def(init<const MATRIX3x3&>())
      .def("__copy__", &generic__copy__<MATRIX3x3>)
      .def("__deepcopy__", &generic__deepcopy__<MATRIX3x3>)

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
      .def(self+=self)
      .def(self-=self)
      .def(self*self)
      .def(self/double())
      .def(self*double())
      .def(self/=double())
      .def(self*=double())
      .def(double()*self)
      .def(self*VECTOR())

      .def("tensor_product",&MATRIX3x3::tensor_product)
      .def("Determinant",&MATRIX3x3::Determinant)
      .def("tr",&MATRIX3x3::tr)
      .def("get_vectors",&MATRIX3x3::get_vectors)
      .def("inverse",&MATRIX3x3::inverse)
      .def("transpose",&MATRIX3x3::transpose)
      .def("T",&MATRIX3x3::T)
      .def("eigen",expt_eigen_v1)
      .def("eigen",expt_eigen_v2)
      .def("identity",&MATRIX3x3::identity)
      .def("diag",expt_diag_v1)
      .def("diag",expt_diag_v2)
      .def("skew",&MATRIX3x3::skew)
      .def("Rx",&MATRIX3x3::Rx)
      .def("Ry",&MATRIX3x3::Ry)
      .def("Rz",&MATRIX3x3::Rz)
      .def("Rotation",&MATRIX3x3::Rotation)

  ;

  class_< MATRIX3x3List >("MATRIX3x3List")
      .def(vector_indexing_suite< MATRIX3x3List >())
  ;

  class_< MATRIX3x3Map >("MATRIX3x3Map")
      .def(vector_indexing_suite< MATRIX3x3Map >())
  ;


}

void export_QUATERNION(){

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
      .def(self*VECTOR())
      .def(VECTOR()*self)
      .def(MATRIX()*self)

      .def("inverse", &QUATERNION::inverse)
      .def("conj", &QUATERNION::conj)
      .def("sqal", &QUATERNION::sqal)
      .def("vect", &QUATERNION::vect)
      .def("norm", &QUATERNION::norm)
      .def("mod", &QUATERNION::mod)
      .def("normalize", &QUATERNION::normalize)

  ;

  class_< QUATERNIONList >("QUATERNIONList")
      .def(vector_indexing_suite< QUATERNIONList >())
  ;

  class_< QUATERNIONMap >("QUATERNIONMap")
      .def(vector_indexing_suite< QUATERNIONMap >())
  ;


}


void export_FT(){

  void (*expt_dft_v1)(CMATRIX& in,CMATRIX& out) = &dft;
  void (*expt_inv_dft_v1)(CMATRIX& in,CMATRIX& out) = &inv_dft;

  def("dft", expt_dft_v1);
  def("inv_dft", expt_inv_dft_v1);


  void (*expt_cft_v1)(CMATRIX& in,CMATRIX& out,double xmin,double dx) = &cft;
  void (*expt_inv_cft_v1)(CMATRIX& in,CMATRIX& out,double xmin,double dx) = &inv_cft;

  void (*expt_cft1_v1)(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx) = &cft1;
  void (*expt_inv_cft1_v1)(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx) = &inv_cft1;

//  The functions below are commented because the inv_CFT[ CFT ] is not the identity, perhaps because dk 
//  would need to satisfy some additional conditions
//  void (*expt_cft2_v1)(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx,double dk) = &cft2;
//  void (*expt_inv_cft2_v1)(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx,double dk) = &inv_cft2;


  def("cft", expt_cft_v1);
  def("cft", expt_cft1_v1);
  def("inv_cft", expt_inv_cft_v1);
  def("inv_cft", expt_inv_cft1_v1);

//  The function below are commented because the inv_CFT[ CFT ] is not the identity, perhaps because dk 
//  would need to satisfy some additional conditions
//  def("cft", expt_cft2_v1);
//  def("inv_cft", expt_inv_cft2_v1);


  void (*expt_cft1_2D_v1)(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy) = &cft1_2D;
  void (*expt_inv_cft1_2D_v1)(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy) = &inv_cft1_2D;

  def("cft_2D", expt_cft1_2D_v1);
  def("inv_cft_2D", expt_inv_cft1_2D_v1);


  void (*expt_convolve_v1)(CMATRIX& f,CMATRIX& g, CMATRIX& conv,double dx) = &convolve;
  void (*expt_convolve_2D_v1)(CMATRIX& f,CMATRIX& g, CMATRIX& conv,double dx,double dy) = &convolve_2D;

  def("convolve", expt_convolve_v1);
  def("convolve_2D", expt_convolve_2D_v1);


  //-------- Fast Fourier Transforms -------------
  void (*expt_cfft1_v1)(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx) = &cfft1;  
  void (*expt_inv_cfft1_v1)(CMATRIX& in,CMATRIX& out,double xmin,double kmin,double dx) = &inv_cfft1;
  def("cfft", expt_cfft1_v1);
  def("inv_cfft", expt_inv_cfft1_v1);


  void (*expt_cfft1_2D_v1)(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy) = &cfft1_2D;
  void (*expt_inv_cfft1_2D_v1)(CMATRIX& in, CMATRIX& out,double xmin,double ymin, double kxmin, double kymin, double dx, double dy) = &inv_cfft1_2D;
  def("cfft_2D", expt_cfft1_2D_v1);
  def("inv_cfft_2D", expt_inv_cfft1_2D_v1);


}

void export_permutations(){


  vector<int> (*expt_id_permutation_v1)(int sz) = &id_permutation;
  def("id_permutation", expt_id_permutation_v1);

  vector<int> (*expt_inverse_permutation_v1)(vector<int>& perm) = &inverse_permutation;
  def("inverse_permutation", expt_inverse_permutation_v1);

  vector<int> (*expt_composite_permutation_v1)(vector<int>& perm_t, vector<int>& perm_cum) = &composite_permutation;
  def("composite_permutation", expt_composite_permutation_v1);


  void (*expt_update_permutation_v1)(vector<int>& perm_t, vector<int>& perm_cum) = &update_permutation;
  def("update_permutation", expt_update_permutation_v1);


  void (*expt_check_permutation_v1)(vector<int>& perm, int n) = &check_permutation;
  def("check_permutation", expt_check_permutation_v1);


}



struct intList_pickle_suite : boost::python::pickle_suite
{
/** 
  Figured out it thanks to these sources:

  https://stackoverflow.com/questions/37180878/append-element-to-boostpythontuple

  https://valelab4.ucsf.edu/svn/3rdpartypublic/boost/libs/python/doc/v2/pickle.html


*/

  static boost::python::tuple getstate(const intList& w)
  {
    boost::python::tuple res; 

    int sz = w.size();
    for(int i=0; i<sz; i++){
        res = boost::python::extract< boost::python::tuple >( res +  boost::python::make_tuple(w[i]) );
    }
      
    return res;

  }

  static void setstate(intList& w, boost::python::tuple state)
  {

    int sz = len(state);

    w = vector<int>(sz);
    for(int i=0; i<sz; i++){   
      w[i] = extract<int>(state[i]);  
    }
  }

};



void export_linalg_objects(){
/** 
  \brief Exporter of the liblinalg classes and functions

*/

  ///============== 1D data ======================

  class_< intList >("intList")
      .def(vector_indexing_suite< intList >())
      .def_pickle(intList_pickle_suite())
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



  ///============== 2D data ======================

  class_< intList2 >("intList2")
      .def(vector_indexing_suite< intList2 >())
  ;

  class_< floatList2 >("floatList2")
      .def(vector_indexing_suite< floatList2 >())
  ;

  class_< doubleList2 >("doubleList2")
      .def(vector_indexing_suite< doubleList2 >())
  ;

  class_< complexList2 >("complexList2")
      .def(vector_indexing_suite< complexList2 >())
  ;


  ///============== 3D data ======================

  class_< intList3 >("intList3")
      .def(vector_indexing_suite< intList3 >())
  ;

  class_< floatList3 >("floatList3")
      .def(vector_indexing_suite< floatList3 >())
  ;

  class_< doubleList3 >("doubleList3")
      .def(vector_indexing_suite< doubleList3 >())
  ;

  class_< complexList3 >("complexList3")
      .def(vector_indexing_suite< complexList3 >())
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


  export_IMATRIX();
  export_MATRIX();
  export_CMATRIX();

  export_VECTOR();
  export_MATRIX3x3();
  export_QUATERNION();

  export_FT();
  export_permutations();



  void (*expt_MATRIX_TO_QUATERNION_v1)(MATRIX&,QUATERNION&) = &MATRIX_TO_QUATERNION;
  void (*expt_MATRIX_TO_QUATERNION_v2)(MATRIX3x3&,QUATERNION&) = &MATRIX_TO_QUATERNION;

  void (*expt_QUATERNION_TO_MATRIX_v1)(QUATERNION&,MATRIX&) = &QUATERNION_TO_MATRIX;
  void (*expt_QUATERNION_TO_MATRIX_v2)(QUATERNION&,MATRIX3x3&) = &QUATERNION_TO_MATRIX;

  def("MATRIX_TO_QUATERNION", expt_MATRIX_TO_QUATERNION_v1);
  def("MATRIX_TO_QUATERNION", expt_MATRIX_TO_QUATERNION_v2);
  def("QUATERNION_TO_MATRIX", expt_QUATERNION_TO_MATRIX_v1);
  def("QUATERNION_TO_MATRIX", expt_QUATERNION_TO_MATRIX_v2);







}// export_linalg_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyglinalg){
#else
BOOST_PYTHON_MODULE(liblinalg){
#endif

  export_linalg_objects();

}




}// namespace liblinalg
}// namespace liblibra




