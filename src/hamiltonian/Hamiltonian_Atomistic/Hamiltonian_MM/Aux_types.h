#ifndef Aux_types_H
#define Aux_types_H


/********************************************************
 This header will contain auxiliary data types for this 
 package. The data types are quite specific - for special
 purposes, so they are not part of more generic headers.
*********************************************************/

struct triple{
  int is_central;
  int n1,n2,n3;
};

struct quartet{
  int is_central;
  int j;        // index of the atom with which another one is interacting
  int n1,n2,n3; // translation vector of atom j
};

struct excl_scale{
  int at_indx1, at_indx2;
  double scale;
};


#endif // Aux_types_H
