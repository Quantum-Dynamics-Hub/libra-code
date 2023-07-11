/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef DATA_H
#define DATA_H

#if defined(USING_PCH)
#include "../pch.h"
#else

#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#endif 

/// liblibra namespace
namespace liblibra{

using namespace boost::python;
using namespace std;


/// libdata namespace
namespace libdata{

//========================= Class DATA ====================================
// This class implements data analysis methods and in fact a collection of
// statistical analysis tools
// A tool may also be convenient for operations with 1D arrays of numerical data

class DATA{

  public:

  std::vector<double> Data;

  // Data properties are:
  // Estimators
  double ave;             int is_ave;        // average
  double var;             int is_var;        // variance
  double sd;              int is_sd;         // sample standard deviation (unbiased estimator for variance)
  double se;              int is_se;         // srandard deviation of mean
  double mse;             int is_mse;        // mean square error
  double mae;             int is_mae;        // mean absolute error
  double rmse;            int is_rmse;       // root mean square error

  // Minimax
  double min_val;         int is_min_val;    // minimal value
  int min_indx;           int is_min_indx;   // index of the minimal entry
  double max_val;         int is_max_val;        // maximal value
  int max_indx;           int is_max_indx;   // index of the maximal entry

  // Operations
  // Data scaling
  double scale_factor;    int is_scale_factor; // if we scale initial data
  double shift_amount;    int is_shift_amount; // if we shift initial data

  // Constructors
  DATA(){
     is_ave = 0;
     is_var = 0;
     is_sd  = 0;
     is_se  = 0;
     is_mse = 0;
     is_mae = 0;
     is_rmse= 0;

     is_min_val = 0;
     is_min_indx = 0;
     is_max_val = 0;
     is_max_indx = 0;

     scale_factor = 1.0;  is_scale_factor = 1;
     shift_amount = 0.0;  is_shift_amount = 1;

  }
  DATA(vector<double>);
  DATA(int,double*);
  DATA(boost::python::list);

  // Copy constructor
  DATA(const DATA& d);

  // Destructors
  ~DATA();


  // Overloaded operators
  DATA& operator=(const DATA&);

  friend int operator == (const DATA& d1, const DATA& d2);
  friend int operator != (const DATA& d1, const DATA& d2);

  // Data Manipulation
  int LinearTransformData(double,double);
  int invLinearTransformData();
  int ScaleData(double);
  int ScaleData(double,double);
  int ShiftData(double);
  int NormalizeData();

  // For data interface
  int PutData(vector<double>&); // from DATA to vector<double>
  int GetData(vector<double>&); // from vector<double> to DATA - similar to constructor

  // Descriptive Statistics
  int Calculate_Estimators(double&,double&,double&,double&,double&,double&,double&);
  int Calculate_Estimators();
  int Calculate_MiniMax(double&,int&,double&,int&);
  int Calculate_MiniMax();
  int Calculate_Distribution(vector<double>&,vector<double>&,vector<double>&);
  boost::python::list Calculate_Distribution(boost::python::list Interval);

  // Regression
  //int Lin_Regression(int, double& ,double& ,double&, double&);


};

  typedef std::vector<DATA> DATAList;

}// namespace libdata
}// namespace liblibra

#endif // DATA_H
