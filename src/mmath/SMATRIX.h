#ifndef SMATRIX_H
#define SMATRIX_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <vector>

using namespace std;

//============================================================
// Forward declared dependencies
class VECTOR;
class SMATRIX;

//===============================================================================

class SMATRIX{
        
public:
  // Data; 
  int num_of_rows;
  int num_of_cols;
  int num_of_elems; // !=  num_of_rows x num_of_cols

  vector<int> indx1; // first (row) index
  vector<int> indx2; // second (column) index
  vector<double> M;  // element

    
  // Constructors;
  SMATRIX(){ num_of_elems = 0; }// default constructor
  SMATRIX(const SMATRIX& ob);   // Copy constructor;
 ~SMATRIX();


  // Growers
  void reserve(int n){
    indx1.reserve(n); indx2.reserve(n); M.reserve(n);
  }

  void init(int n){
    indx1 = vector<int>(n,-1);
    indx2 = vector<int>(n,-1);
    M = vector<double>(n,0.0);
  }
  
  void add_element(int i,int j, double x){ 
  // Complexity: constant
    indx1.push_back(i);
    indx2.push_back(j);
    M.push_back(x);

    num_of_elems++;
  }

  int find_element(int i,int j){   
  // Complexity: linear in size of the matrix
    int res = -1;
    for(int k=0;k<num_of_elems;k++){
      if(indx1[k]==i && indx2[k]==j){ res = k; }
    }
    return res;
  }



  SMATRIX operator-();  
  SMATRIX operator*(const SMATRIX& ob);
  SMATRIX operator+(SMATRIX ob);
  SMATRIX operator-(SMATRIX ob);
  void operator+=(SMATRIX ob);
  void operator-=(SMATRIX ob);
  SMATRIX operator/(double num);
  SMATRIX operator=(SMATRIX ob);
  SMATRIX operator=(double num);



};


#endif // SMATRIX_H


