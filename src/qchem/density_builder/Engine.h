#ifndef ENGINE_H
#define ENGINE_H

#include "Mathematics.h"
#include "MOAO.h"
#include "ERI_integrals.h"
using namespace std;

#include "AO.h"
#include "Control_Parameters.h"
#include "Timer.h"



// Defined in Engine.cpp
void solve_eigen(int, MATRIX*, MATRIX*, MATRIX*, MATRIX*);
void solve_eigen(int, MATRIX*, MATRIX*, MATRIX*);

int merge_sort(vector< pair<int,double> >&, vector< pair<int,double> >&);
void order_bands(int, MATRIX*, vector< pair<int,double> >&);

double fermi_integral(vector< pair<int,double> >&, double, double);
double fermi_energy(vector< pair<int,double> >& bnds, double, double);
double population(double, double, double);
void populate(int, int, int, double, int,vector< pair<int,double> >&,vector< pair<int,double> >&);


void compute_density_matrix(vector< pair<int,double> >&, MATRIX*, MATRIX*);

void annihilate(int, int, MATRIX*, MATRIX*, MATRIX*, MATRIX*);
void annihilate(int, int, MATRIX*, MATRIX*);


void Fock_to_P(int,int, int, double, std::string, int,MATRIX*, MATRIX*, MATRIX*, MATRIX*,
               vector< pair<int,double> >&, vector< pair<int,double> >&, MATRIX*, vector<Timer>& );

void update_Mull_orb_pop(MATRIX*, MATRIX*, vector<double>&, vector<double>&);

void update_Mull_charges(vector<int>&, vector<int>&, vector<vector<int> >&,vector<double>&,
                         vector<double>&, vector<double>&, vector<double>&, vector<double>&);

double energy_elec(int,MATRIX*,MATRIX*,MATRIX*);
double energy_elec(int Norb, 
                   MATRIX* P_alp, MATRIX* P_bet, 
                   MATRIX* Hao_alp, MATRIX* Hao_bet,
                   MATRIX* Fao_alp, MATRIX* Fao_bet,
                   MATRIX* dFao_alp_dP_alp, MATRIX* dFao_alp_dP_bet,
                   MATRIX* dFao_bet_dP_alp, MATRIX* dFao_bet_dP_bet,
                   MATRIX* temp
                  );

void show_bands(int, int, vector< pair<int,double> >&, vector< pair<int,double> >&);


#endif // ENGINE_H

