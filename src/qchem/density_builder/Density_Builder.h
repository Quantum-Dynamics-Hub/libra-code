#ifndef DENSITY_BUILDER_H
#define DENSITY_BUILDER_H


#include "Mathematics.h"
using namespace std;

void pop_submatrix(MATRIX*,MATRIX*,vector<int>&);
void push_submatrix(MATRIX*,MATRIX*,vector<int>&);

void handle_subset(MATRIX*, MATRIX*, MATRIX*, 
                   MATRIX*, MATRIX*, MATRIX*, MATRIX*, MATRIX*,
                   vector< pair<int,double> >&, vector< pair<int,double> >&,                  
                   vector<int>&, int, int, int, std::string);

void handle_n_subsets(MATRIX*, MATRIX*, MATRIX*, 
                      int, MATRIX**, MATRIX**, MATRIX**, MATRIX**, MATRIX**,
                      vector< vector< pair<int,double> > >&, vector< vector< pair<int,double> > >&, 
                      vector< vector<int> >&, vector<int>&, vector<int>&, int, std::string, int);



#endif // DENSITY_BUILDER_H