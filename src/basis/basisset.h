/*********************************************************************************
* Copyright (C) 2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
* This file is originally a part of the ErgoSCF code (see the info below). 
* Now, it is meant to use in the Libra code.
* The original file has been significantly modified: 
  1) the "struct" data types are now classes.
  2) the global preprocessor definitions are deprecated - now use
     the class constructor parameters to provide it. 
  3) we use the vector<> container instead of the pointers to the arrays. 
  4) the comparison operator are implemented in each class - this is needed
     to expose the datatypes to Python via indexing suite
  5) the do_output functions is temporarity deprecated - in the future, we are 
     likely to use the corresponding streams. Also, various verifications will 
     go silently, not producing the output message
  6) the free function for reading the basis set info is now a member-function of 
     the uppermost class
  7) ... this function now takes fewer arguments - the corresponding "magic" with
     the filenames constructions is commented, cause it's not needed anymore
*
*********************************************************************************/


/* Ergo, version 3.3, a program for linear scaling electronic structure
 * calculations.
 * Copyright (C) 2013 Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Primary academic reference:
 * Kohn−Sham Density Functional Theory Electronic Structure Calculations 
 * with Linearly Scaling Computational Time and Memory Usage,
 * Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek,
 * J. Chem. Theory Comput. 7, 340 (2011),
 * <http://dx.doi.org/10.1021/ct100611z>
 * 
 * For further information about Ergo, see <http://www.ergoscf.org>.
 */

/* -*-mode:c; c-style:k&r; c-basic-offset:4; indent-tabs-mode: nil -*- */
#ifndef BASISSET_HEADER
#define BASISSET_HEADER


#if defined(USING_PCH)
#include "../pch.h"
#else
#include <vector>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#endif
#include "../realtype.h"




/// We incorporate all the definitions into the liblibra namespace
/// liblibra namespace
namespace liblibra{

/// Also, to keep everything organized, lets also create a libergoescf namespace
/// to keep all the original developments in there
/// libergoescf namespace 
namespace libergoscf{


/* Here we specify which kinds of basis functions are allowed.
   0: s
   1: s, p
   2: s, p, d
   3: s, p, d, f
   4: s, p, d, f, g
   etc.
*/
#define BASIS_FUNC_POLY_MAX_DEGREE 5




using namespace std;

class basisset_shell_struct{

public:
/** Represents a single basis function as a contraction of GTOs */
  int type;
  int contrCount;
  int shell_ID;
  std::vector<ergo_real> exponentList; 
  std::vector<ergo_real> coeffList; 

  basisset_shell_struct(int max_no_of_contr){
    exponentList = std::vector<ergo_real>(max_no_of_contr, 0.0);
    coeffList = std::vector<ergo_real>(max_no_of_contr, 0.0);
  }
  basisset_shell_struct(const basisset_shell_struct& ob){
    type = ob.type; 
    contrCount = ob.contrCount;
    shell_ID = ob.shell_ID;
    exponentList = ob.exponentList;
    coeffList = ob.coeffList;
  }
  ~basisset_shell_struct(){
    exponentList.clear();
    coeffList.clear();
  }

  ///< Comparison: true if the ID sare different
  friend int operator !=(const basisset_shell_struct& m1, const basisset_shell_struct& m2){
    if(m1.type != m2.type){ return 1; } 
    if(m1.contrCount != m2.contrCount){ return 1; } 
    if(m1.shell_ID != m2.shell_ID){ return 1; } 
    if(m1.exponentList.size() != m2.exponentList.size()){ return 1; } 
    if(m1.coeffList.size() != m2.coeffList.size()){ return 1; } 
    
    int sz = m1.exponentList.size();
    for(int i=0;i<sz;i++){
      if(m1.exponentList[i] != m2.exponentList[i]) { return 1; }
      if(m1.coeffList[i] != m2.coeffList[i]) { return 1; }
    }
    sz = m1.coeffList.size();
    for(int i=0;i<sz;i++){
      if(m1.coeffList[i] != m2.coeffList[i]) { return 1; }
    }
    return 0;
  }

  ///< Comparison: 
  friend int operator ==(const basisset_shell_struct& m1, const basisset_shell_struct& m2){
    return !(m1!=m2);
  }



}; 


typedef std::vector<basisset_shell_struct> basisset_shell_structList;  
typedef std::vector< vector<basisset_shell_struct> > basisset_shell_structMap;  



class basisset_atom_struct{
/** Represents all atomic orbitals on a given atom */

public:

  int noOfShells;
  std::vector<basisset_shell_struct> shells;

  basisset_atom_struct(int max_no_of_shells_per_atom, int max_no_of_contr){
    noOfShells = 0;
    shells = std::vector<basisset_shell_struct>(max_no_of_shells_per_atom, basisset_shell_struct(max_no_of_contr) );
  }
  basisset_atom_struct(const basisset_atom_struct& ob){
    noOfShells = ob.noOfShells;
    shells = ob.shells;
  }
  ~basisset_atom_struct(){   shells.clear();  }

  ///< Comparison: 
  friend int operator !=(const basisset_atom_struct& m1, const basisset_atom_struct& m2){

    if(m1.noOfShells != m2.noOfShells){ return 1; } 
    
    int sz = m1.shells.size();
    for(int i=0;i<sz;i++){
      if(m1.shells[i] != m2.shells[i]) { return 1; }
    }

    return 0;
  }
  ///< Comparison: 
  friend int operator ==(const basisset_atom_struct& m1, const basisset_atom_struct& m2){
    return !(m1!=m2);
  }
  


};

typedef std::vector<basisset_atom_struct> basisset_atom_structList;  
typedef std::vector< std::vector<basisset_atom_struct> > basisset_atom_structMap;  



class basisset_struct{

  int MAX_NO_OF_CONTR;            ///< Max contraction size
  int MAX_NO_OF_SHELLS_PER_ATOM;  ///< Max number of AO orbitals per atom

  /// Max number of atomic types for which the basis is
  /// defined in the same file
  int MAX_NO_OF_ATOM_TYPES;

public:
  std::vector<basisset_atom_struct> atoms;



  basisset_struct(int max_no_of_atom_types, int max_no_of_shells_per_atom, int max_no_of_contr){ 
    
    MAX_NO_OF_CONTR = max_no_of_contr;    
    MAX_NO_OF_SHELLS_PER_ATOM = max_no_of_shells_per_atom;
    MAX_NO_OF_ATOM_TYPES = max_no_of_atom_types;

    atoms = std::vector<basisset_atom_struct>(max_no_of_atom_types, basisset_atom_struct(max_no_of_shells_per_atom, max_no_of_contr) );
  }

  int read_basisset_file(std::string fileName,int print_raw);

} ;

typedef std::vector<basisset_struct> basisset_structList;  
typedef std::vector< vector<basisset_struct> > basisset_structMap;  



/*
int read_basisset_file(basisset_struct* result, 
		       const char* fileName,
		       int dirc, 
		       const char *dirv[],
		       int print_raw);
*/

}// libergoscf
}// liblibra




#endif // BASISSET_HEADER
