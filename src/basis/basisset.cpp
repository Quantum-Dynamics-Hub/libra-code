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

/** @file basisset.cc

    \brief Code for parsing a text file specifying a basisset.

    @author: Elias Rudberg <em>responsible</em>. 
*/
/* -*-mode:c; indent-tabs-mode: nil -*- */
/* basisset.c: provides read_basisset_file() which creates a
   basisset_struct from a data contained in a file */

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#endif

#include "basisset.h"
//#include "output.h"
//#include "memorymanag.h"


/// We incorporate all the definitions into the liblibra namespace
/// liblibra namespace
namespace liblibra{

/// Also, to keep everything organized, lets also create a libergoescf namespace
/// to keep all the original developments in there
/// libergoescf namespace 
namespace libergoscf{


static void
remove_zeros(basisset_atom_struct* currAtom,
             int shellBaseIndex, int noOfShellsCurrBatch) {
  /**  Remove zero coefficients. */

  for(int i = 0; i < noOfShellsCurrBatch; i++) {
    int cnt = 0;

   /**Basically, we pack all non-zero contraction coefficients such that there are 
   no zeroes in a congruent region
   */

    for(int j = 0; j < currAtom->shells[shellBaseIndex+i].contrCount; j++) {

      ergo_real currCoeff    = currAtom->shells[shellBaseIndex+i].coeffList[j];
      ergo_real currExponent = currAtom->shells[shellBaseIndex+i].exponentList[j];

      if(currCoeff != 0) {
	currAtom->shells[shellBaseIndex+i].coeffList[cnt] = currCoeff;
	currAtom->shells[shellBaseIndex+i].exponentList[cnt] = currExponent;
	cnt++;

      }

    } /*  END FOR j */
    currAtom->shells[shellBaseIndex+i].contrCount = cnt;

  } /*  END FOR i         */

}


/** read_basisset_file: reads a basis set from fileName. The basis set
   exponents and contraction coefficients are placed in result.  The
   reading procedure is bit convoluted because the basis set file
   follows the Fortran syntax, with wrapping and skipping empty
   elements. We basically need to emulate fortran `format'
   statement. There is one "but" though: AhlrichsDenFit does not
   follow the format syntax so it will/may be misread by eg. dalton. What
   a mess...
   
   The parser is implemented as a state machine. It still cannot read
   ANO-type basis sets...
*/

/* The original signature
int 
read_basisset_file(basisset_struct* result, const char* fileName,
                   int dirc, const char *dirv[],
                   int print_raw)
*/ 

int basisset_struct::read_basisset_file(std::string fileName,int print_raw){

  enum { END_PARSING, ATOM_EXPECTED, SHELL_EXPECTED,
         SHELL_OR_ATOM_EXPECTED, CONTRACTION_BLOCK_EXPECTED } state;
  int uncontracted = 0;
  char line[256];
  basisset_atom_struct* currAtom = NULL;
  int atomType = -1;
  int spdf = -1;
  int shellBaseIndex = -1;
  int expNo = -1;
  FILE* f = NULL;


  f = fopen(fileName.c_str(), "rt");

/*
  if(!fileName) {
//    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: fileName == NULL.");
    return -1;
  }
  if(fileName[0] == '/')
    f = fopen(fileName, "rt");
  else {
    for(int i = 0; i < dirc; i++) {
      char ffname[256];
      int len = strlen(dirv[i]);
      strncpy(ffname, dirv[i], sizeof(ffname));
      if(len>0 && ffname[len-1] != '/')
        strncat(ffname, "/", sizeof(ffname)-1-len);
      strncat(ffname, fileName, sizeof(ffname)-2-len);
//      do_output(LOG_CAT_WARNING, LOG_AREA_INTEGRALS, 
//                "Trying basis set file '%s'...", ffname);      

      if( (f=fopen(ffname, "rt")) != NULL)
        break;
    }
  }
*/
      
  if(f == NULL)
    {
//      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error opening file '%s'", fileName);
      return -1;
    }
  
  /*  now create basis set by reading buf2 */

//  memset(result, 0, sizeof(basisset_struct));
  int noOfAtomTypes = 0;
  state = ATOM_EXPECTED;
  int lineNo = 0;
  int lineConsumed = 1;
  ergo_real currExponent = 0;

  /* start global parsing loop */
  do {
    int dummy;
    if(lineConsumed) {
      if(fgets(line, sizeof(line), f) == NULL) {
        state = END_PARSING;
        break;
      }
      lineConsumed = 0;
      lineNo++;
    }


    for(int cc = strlen(line)-1; cc>=0 && isspace(line[cc]); cc--)
      line[cc] = '\0';
        
    if(line[0] == '$' || line[0] == '!' || line[0] == '#'||
       line[0] == '*' || line[0] == '\0'|| line[0] == '\n') {
      lineConsumed = 1; /* skip the comment and move on */
      continue;
    }

    switch(state) {
    case ATOM_EXPECTED:
      if(line[0] == 'a' || line[0] == 'A') {
        noOfAtomTypes++;
        atomType = atoi(line+1);
//        int atomType = atoi(line+1);
//        if(print_raw) 
//          do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "Basis set for atom of Z=%d", atomType);

        state = SHELL_EXPECTED;
        if(atomType <= 0){
//            do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: (atomType <= 0) "
//                      " in line %d %s\n", lineNo, fileName);
          return -1;
        }

        if(atomType >= MAX_NO_OF_ATOM_TYPES){
//            do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: "
//                      "(atomType >= MAX_NO_OF_ATOM_TYPES) in line %d %s\n",
//                      lineNo, fileName);
          return -1;
        }

        currAtom = &atoms[atomType];

        //cout<<"currAtom->noOfShells = "<<currAtom->noOfShells<<endl;
        //cout<<result->atoms[atomType].noOfShells<<endl;
        spdf = 0;
        shellBaseIndex = 0;
      }
      lineConsumed = 1;
      break;

    case SHELL_OR_ATOM_EXPECTED:
      if(line[0] == 'a' || line[0] == 'A') {
        state = ATOM_EXPECTED;
        /* fininalize current atom data */
        if(currAtom == NULL || shellBaseIndex < 0)
          return -1;
        currAtom->noOfShells = shellBaseIndex;
        break;
      } /* else fall through */

    case SHELL_EXPECTED:
      if(shellBaseIndex < 0)
        return -1;
      int noOfExponents, noOfShellsCurrBatch;
      if(sscanf(line, "%d %d %d",
                &noOfExponents, &noOfShellsCurrBatch, &dummy ) != 3)
        {
//          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: "
//                    "Shell data expected in line %d:\n%s", lineNo, line);
          return -1;
        }
      if(noOfExponents <= 0)
        {
//          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: "
//                    "(noOfExponents <= 0) in line %d %s\n", lineNo, fileName);
          return -1;
        }
      if(noOfExponents >= MAX_NO_OF_CONTR)
        {
//          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: "
//                    "(noOfExponents >= MAX_NO_OF_CONTR) in line %d\n",
//                    lineNo);
          return -1;
        }
      if(noOfShellsCurrBatch < 0)
        {
//          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: "
//                    "(noOfShellsCurrBatch < 0) in line %d", lineNo);
          return -1;
        }
      if(noOfShellsCurrBatch == 0) {
        /*  special case: uncontracted, expect only one column */
        noOfShellsCurrBatch = noOfExponents;
        uncontracted = 1;
      } else uncontracted = 0;

      if(shellBaseIndex + noOfShellsCurrBatch >= MAX_NO_OF_SHELLS_PER_ATOM)
        {
//          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: too many shells.");
          return -1;
        }
      /* initialize shell data. Set the contraction count to its upper
         limit. remove_zeros() will later check for a better, lower
         value. */
      for(int i = 0; i < noOfShellsCurrBatch; i++) {
        cout<<"in for i\n";
	if(currAtom == NULL || spdf < 0){
	  return -1;
        }
        atoms[atomType].shells[shellBaseIndex+i].type = spdf;
        atoms[atomType].shells[shellBaseIndex+i].contrCount = noOfExponents;
      }
      expNo = 0;
      state = CONTRACTION_BLOCK_EXPECTED;
//      if(print_raw)
//        do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
//                  "Block for L=%d primitives: %d contracted: %d",
//                  spdf, noOfExponents, noOfShellsCurrBatch);

      lineConsumed = 1;
      break;

    case CONTRACTION_BLOCK_EXPECTED:
      currExponent = atof(line);
      if(currExponent <= 0) {
//	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (currExponent <= 0) in line %d", lineNo);
	return -1;
      }
      if(currAtom == NULL || shellBaseIndex < 0 || expNo < 0)
        return -1;
      if(uncontracted) {
        for(int i = 0; i < noOfShellsCurrBatch; i++) {
          currAtom->shells[shellBaseIndex+i].exponentList[expNo] =
            currExponent;
          currAtom->shells[shellBaseIndex+i].coeffList[expNo] =
            i == expNo ? 1.0 : 0.0;
        }
      } else {
        int idx = 0;
        /* skip exponent */
        while(line[idx] && isspace(line[idx]))  idx++;
        for(int i = 0; i < noOfShellsCurrBatch; i++) {
          currAtom->shells[shellBaseIndex+i].exponentList[expNo] =
            currExponent;
          while(line[idx] && !isspace(line[idx])) idx++;
          while(line[idx] && isspace(line[idx]))  idx++;
          if( !line[idx] ) {
	    /* Second line begins when we are about to read 7th
	       contraction coefficient (i=6), third line for 14th
	       (i=13), fourth for i=20 etc. If this pattern does not
	       match, warn the user. */
//	    if( i != 6 && (i+1) % 7 != 0 )
//	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "WARN: line %d has trailing data: '%s'"
//			"non-conformant basis set file for i=%d.",
//			lineNo, line + idx, i);
	    if(fgets(line, sizeof(line), f) == NULL) {
//	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "reading error when continuing shell data.");
	      return -1;
	    }
	    lineNo++;
	    idx = 0;
	    while(line[idx] && isspace(line[idx]))  idx++;
	  }
          currAtom->shells[shellBaseIndex+i].coeffList[expNo] =
            atof(line + idx);
        }  /*  END FOR i */
      }
      if(print_raw) {
        char line[256], eee[20];
        line[0] = '\0';
        for(int i = 0; i<noOfShellsCurrBatch; i++) {
          sprintf(eee, "%10.5f",
                  (double)currAtom->shells[shellBaseIndex+i].coeffList[expNo]);
          strcat(line, eee);
        }
//        do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
//                  "%d %12.6f: %s", expNo, (double)currExponent, line);
      }
      if(++expNo == noOfExponents) {
        remove_zeros(currAtom, shellBaseIndex, noOfShellsCurrBatch);
        shellBaseIndex += noOfShellsCurrBatch;
        spdf++;
        state = SHELL_OR_ATOM_EXPECTED;
      }
      lineConsumed = 1;
      break;
    case END_PARSING:
      /*  do nothing */
      break;
    }
  } while(state != END_PARSING);
  fclose(f);

  /* fininalize current atom data */
  if(currAtom == NULL || shellBaseIndex < 0){
    return -1;
  }
  currAtom->noOfShells = shellBaseIndex;


  /*  Postprocessing... */
  /*  set shell ID for each shell in basis set */
  int currShellID = 0;
  for(int i = 0; i < MAX_NO_OF_ATOM_TYPES; i++) {
    int noOfShells = atoms[i].noOfShells;

    for(int j = 0; j < noOfShells; j++) {
      currShellID++;
      atoms[i].shells[j].shell_ID = currShellID;

    } /*  END FOR j */
  } /*  END FOR i */
//  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "total number of shells in basis set: %i", currShellID);
//  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "Basis set file '%s' processed OK, noOfAtomTypes = %i",
//            fileName, noOfAtomTypes);

  cout<<"1.5\n";

  return 0;
}




}// libergoscf
}// liblibra

