/*********************************************************************************
* Copyright (C) 2014 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
 \file Print.cpp
 \brief Implementation of the methods for printing properties

*/
/****************************************************************************
  This file contains following functions:

  void print_xyz(std::string filename,Nuclear& mol)
  void print_xyz(std::string filename,Nuclear& mol,double Eelec,double Enucl,double Etot)
  void print_xyzq(std::string filename,Nuclear& mol)
  void print_xyzq(std::string filename,Nuclear& mol,double Eelec,double Enucl,double Etot)
  void print_el_struct(std::string filename,Electronic* el,double Eelec)
  void print_Mulliken_charge(std::string filename,Nuclear& mol)
  void print_dipole(std::string filename,Electronic* el,double Eelec,MATRIX* mux, MATRIX* muy, MATRIX* muz)
  void print_excitations(std::string filename, Electronic* el, Control_Parameters& prms)

****************************************************************************/

#include <algorithm>

#include "Print.h"
#include "units.h"
#include "Engine.h"


using namespace std;


void print_xyz(std::string filename,Nuclear& mol){
  
  print_xyz(filename,mol,0.0,0.0,0.0);

}

void print_xyz(std::string filename,Nuclear& mol,double Eelec,double Enucl,double Etot){

  ofstream out(filename.c_str(),ios::out);

  out<<mol.Nnucl<<endl;
  out<<"Molecule Eelec= "<<Eelec<<" Enucl= "<<Enucl<<" Etot= "<<Etot<<endl;
  for(int n=0;n<mol.Nnucl;n++){
    out<<mol.at_type[n]<<"  "<<mol.R[n]/Angst<<endl;
  }
  out.close();

}

void print_xyzq(std::string filename,Nuclear& mol){
  
  print_xyzq(filename,mol,0.0,0.0,0.0);

}

void print_xyzq(std::string filename,Nuclear& mol,double Eelec,double Enucl,double Etot){

  ofstream out(filename.c_str(),ios::out);

  double mu_mod;
  VECTOR mu; 

  mol.dipole_moment(mu,mu_mod);

  out<<mol.Nnucl<<endl;
  out<<"Molecule: Eelec= "<<Eelec<<" Enucl= "<<Enucl<<" Etot= "<<Etot<<" Dipole_moment= "<<mu<<" |Dipole_moment|= "<<mu_mod<<endl;

  out<<"Atom_type         X          Y       Z       Mulliken_charge       Nuclear_effective_charge \n";
  for(int n=0;n<mol.Nnucl;n++){
    out<<mol.at_type[n]<<"  "<<mol.R[n]/Angst<<"   "<<mol.Mull_charges_gross[n]<<"   "<<mol.Zeff[n]<<endl;
  }
  out.close();

}


void print_el_struct(std::string filename,Electronic* el,double Eelec,double tol){
// This function prints the electronic structure details - energies of orbitals, HOMO, LUMO, Fermi energy

  ofstream out(filename.c_str(),ios::out); 

  out<<"Eelec= "<<Eelec<<" a.u. = "<<Eelec/eV<<" eV"<<endl;

  double _E_f_alp = fermi_energy(el->bands_alp,1.0 * el->Nocc_alp,1.0);
  double _E_f_bet = fermi_energy(el->bands_bet,1.0 * el->Nocc_bet,1.0);

  cout<<"E_Fermi, alp = "<<_E_f_alp<<endl;

  int homo_alp = 0;
  int homo_bet = 0;

//  double tol = 0.2;
  
  double tot_popa = 0.0;
  double tot_popb = 0.0;
  int first_a = 0;
  int first_b = 0;

  cout<<"Printing smeared populations\n";
  for(int i=0;i<el->Norb;i++){

    double popa = population(el->bands_alp[i].second,_E_f_alp,1.0);
    double popb = population(el->bands_bet[i].second,_E_f_bet,1.0);

    tot_popa += popa; 
    tot_popb += popb;

    if(fabs(tot_popa-el->Nocc_alp)<tol && !first_a) { homo_alp = i; first_a = 1; }
    if(fabs(tot_popb-el->Nocc_bet)<tol && !first_b) { homo_bet = i; first_b = 1; }

    cout<<" i = "<<i<<"tot_popa = "<<tot_popa<<" |tot_popa-Nocc_alp| = "<<fabs(tot_popa-el->Nocc_alp)<<endl;

  }// for i

  cout<<"homo_alp = "<<homo_alp<<endl;


  if((homo_alp+1)<el->Norb){  out<<"LUMO(alp)= "<<el->bands_alp[homo_alp+1].second/eV<<" eV";  }
  else{   out<<"LUMO(alp)= "<<el->bands_alp[el->Norb-1].second/eV<<" eV";  }
  
  out<<"   ";
  if((homo_alp+1)<el->Norb){  out<<"LUMO(bet)= "<<el->bands_bet[homo_bet+1].second/eV<<" eV";  }
  else{   out<<"LUMO(bet)= "<<el->bands_bet[el->Norb-1].second/eV<<" eV";  }
  out<<endl;

  out<<"E_F(alp)= "<<_E_f_alp/eV<<" eV";  out<<"   ";
  out<<"E_F(bet)= "<<_E_f_bet/eV<<" eV";  out<<endl;

  out<<"HOMO(alp)= "<<el->bands_alp[homo_alp].second/eV<<" eV";  out<<"   ";
  out<<"HOMO(bet)= "<<el->bands_bet[homo_bet].second/eV<<" eV";  out<<endl;

  

  for(int i=0;i<el->Norb;i++){
    out<<"i= "<<i<<" orbital_index= "<<el->bands_alp[i].first<<" E= "<<el->bands_alp[i].second/eV<<" eV, occupation_number= "<<el->occ_alp[i].second;
    out<<"   ";
    out<<"i= "<<i<<" orbital_index= "<<el->bands_bet[i].first<<" E= "<<el->bands_bet[i].second/eV<<" eV, occupation_number= "<<el->occ_bet[i].second;
    out<<endl;
  }

  out.close();

}


void print_Mulliken_charge(std::string filename,Nuclear& mol){
// Print Mulliken charges to the file <filename>
// Also compute multipoles, but keep in mind that these multipoles are only approximations 
// the true quantum-mechanical multipoles are computed in void print_multipoles(std::string spin_method, std::string filename)

  ofstream out(filename.c_str(),ios::out);

  VECTOR mu; mu = 0.0;  // dipole moment

  out<<mol.Nnucl<<endl;
  for(int n=0;n<mol.Nnucl;n++){ 

    out<<mol.at_type[n]<<"  Coord(Angst)= "<<mol.R[n]/Angst
                       <<"  Gross_charge= "<<mol.Mull_charges_gross[n]
                       <<"  Net_charge=   "<<mol.Mull_charges_net[n]<<endl;
    mu += mol.R[n] * mol.Mull_charges_gross[n];
  }// for n

  out<<" Dipole moment(total), Debye= "<<mu/Debye<<" , "<<mu.length()/Debye<<"\n";

  out.close();

}


void print_dipole(std::string filename,Electronic* el,double Eelec,MATRIX* mux, MATRIX* muy, MATRIX* muz,VECTOR& mu_nucl){
// This function prints the dipole moments: transition dipole moments, magnitudes, etc.

  int vlvl = 0;

  int i, ii;
  double mxa,mya,mza,mxb,myb,mzb,mtota,mtotb,mtot;
  MATRIX* tmp; tmp = new MATRIX(el->Norb,el->Norb);
  MATRIX* pop; pop = new MATRIX(el->Norb,el->Norb);


  ofstream out(filename.c_str(),ios::out); 

  out<<"Dipole momenta -<xi_i|r-Rcom(frag)|xi_j> - translationally invariant, but not rotationally invariant\n";
  if(vlvl>0){ out<<"AO basis: Transition dipole moment x-component, a.u. \n"<<*mux<<endl; }
  if(vlvl>0){ out<<"AO basis: Transition dipole moment y-component, a.u. \n"<<*muy<<endl; }
  if(vlvl>0){ out<<"AO basis: Transition dipole moment z-component, a.u. \n"<<*muz<<endl; }

  *pop = 0.0; 
  for(i=0;i<el->Norb;i++){ 
    ii = i;//el->occ_alp[i].first;
    pop->M[ii*el->Norb+ii] = el->occ_alp[i].second; 
  }

  *tmp = (*el->C_alp).T() * *mux * *el->C_alp;
  if(vlvl>0){ out<<"MO basis: Transition dipole moment x-component, alp, a.u. \n"<<*tmp<<endl; }  mxa = (*pop * *tmp).tr();
  *tmp = (*el->C_alp).T() * *muy * *el->C_alp;
  if(vlvl>0){ out<<"MO basis: Transition dipole moment y-component, alp, a.u. \n"<<*tmp<<endl; }  mya = (*pop * *tmp).tr();
  *tmp = (*el->C_alp).T() * *muz * *el->C_alp;
  if(vlvl>0){ out<<"MO basis: Transition dipole moment z-component, alp, a.u. \n"<<*tmp<<endl; }  mza = (*pop * *tmp).tr();
  mtota = sqrt(mxa*mxa+mya*mya+mza*mza);
  out<<"Components of the electronic dipole moment for alpha channel, Debye = ( "<<mxa/Debye<<" , "<<mya/Debye<<" , "<<mza/Debye<<" ). Absolute value, Debye = "<<mtota/Debye<<endl;


  *pop = 0.0;
  for(i=0;i<el->Norb;i++){ 
    ii = i;// el->occ_bet[i].first;
    pop->M[ii*el->Norb+ii] = el->occ_bet[i].second; 
  }

  *tmp = (*el->C_bet).T() * *mux * *el->C_bet;
  if(vlvl>0){ out<<"MO basis: Transition dipole moment x-component, bet, a.u. \n"<<*tmp<<endl; }  mxb = (*pop * *tmp).tr(); 
  *tmp = (*el->C_bet).T() * *muy * *el->C_bet;                                                         
  if(vlvl>0){ out<<"MO basis: Transition dipole moment y-component, bet, a.u. \n"<<*tmp<<endl; }  myb = (*pop * *tmp).tr(); 
  *tmp = (*el->C_bet).T() * *muz * *el->C_bet;                                                         
  if(vlvl>0){ out<<"MO basis: Transition dipole moment z-component, bet, a.u. \n"<<*tmp<<endl; }  mzb = (*pop * *tmp).tr(); 
  mtotb = sqrt(mxb*mxb+myb*myb+mzb*mzb);
  out<<"Components of the electronic dipole moment for beta channel, Debye = ( "<<mxb/Debye<<" , "<<myb/Debye<<" , "<<mzb/Debye<<" ). Absolute value, Debye = "<<mtotb/Debye<<endl;

  mtot = sqrt((mxa+mxb)*(mxa+mxb) + (mya+myb)*(mya+myb) + (mza+mzb)*(mza+mzb));
  out<<"Total electronic dipole moment, Debye = ( "<<(mxa+mxb)/Debye<<" , "<<(mya+myb)/Debye<<" , "<<(mza+mzb)/Debye<<" ). Absolute value, Debye = "<<mtot/Debye<<endl;
  out<<"Total nuclear dipole moment, Debye = ( "<<mu_nucl.x/Debye<<" , "<<mu_nucl.y/Debye<<" , "<<mu_nucl.z/Debye<<" ). Absolute value, Debye = "<<mu_nucl.length()/Debye<<endl;
  mtot = sqrt((mu_nucl.x+mxa+mxb)*(mu_nucl.x+mxa+mxb) + (mu_nucl.y+mya+myb)*(mu_nucl.y+mya+myb) + (mu_nucl.z+mza+mzb)*(mu_nucl.z+mza+mzb));
  out<<"Total dipole moment, Debye = ( "<<(mu_nucl.x+mxa+mxb)/Debye<<" , "<<(mu_nucl.y+mya+myb)/Debye<<" , "<<(mu_nucl.z+mza+mzb)/Debye<<" ). Absolute value, Debye = "<<mtot/Debye<<endl;
 
  out<<endl<<endl;


  *tmp = 0.5*( (*el->P_alp) * *mux + *mux * (*el->P_alp) );
  if(vlvl>0){ out<<"Density matrix(symmetrized): Transition dipole moment x-component, alp, a.u. \n"<<*tmp<<endl; }  mxa = tmp->tr();
  *tmp = 0.5*( (*el->P_alp) * *muy + *muy * (*el->P_alp) );
  if(vlvl>0){ out<<"Density matrix(symmetrized): Transition dipole moment y-component, alp, a.u. \n"<<*tmp<<endl; }  mya = tmp->tr();
  *tmp = 0.5*( (*el->P_alp) * *muz + *muz * (*el->P_alp) );
  if(vlvl>0){ out<<"Density matrix(symmetrized): Transition dipole moment z-component, alp, a.u. \n"<<*tmp<<endl; }  mza = tmp->tr();
  mtota = sqrt(mxa*mxa+mya*mya+mza*mza);
  out<<"Components of the electronic dipole moment for alpha channel, Debye = ( "<<mxa/Debye<<" , "<<mya/Debye<<" , "<<mza/Debye<<" ). Absolute value, Debye = "<<mtota/Debye<<endl;


  *tmp = 0.5*( (*el->P_bet) * *mux + *mux * (*el->P_bet) );
  if(vlvl>0){ out<<"Density matrix(symmetrized): Transition dipole moment x-component, bet, a.u. \n"<<*tmp<<endl; }  mxb = tmp->tr();
  *tmp = 0.5*( (*el->P_bet) * *muy + *muy * (*el->P_bet) );
  if(vlvl>0){ out<<"Density matrix(symmetrized): Transition dipole moment y-component, bet, a.u. \n"<<*tmp<<endl; }  myb = tmp->tr();
  *tmp = 0.5*( (*el->P_bet) * *muz + *muz * (*el->P_bet) );
  if(vlvl>0){ out<<"Density matrix(symmetrized): Transition dipole moment z-component, bet, a.u. \n"<<*tmp<<endl; }  mzb = tmp->tr();
  mtotb = sqrt(mxb*mxb+myb*myb+mzb*mzb);
  out<<"Components of the electronic dipole moment for beta channel, Debye = ( "<<mxb/Debye<<" , "<<myb/Debye<<" , "<<mzb/Debye<<" ). Absolute value, Debye = "<<mtotb/Debye<<endl;

  mtot = sqrt((mxa+mxb)*(mxa+mxb) + (mya+myb)*(mya+myb) + (mza+mzb)*(mza+mzb));
  out<<"Total electronic dipole moment, Debye = ( "<<(mxa+mxb)/Debye<<" , "<<(mya+myb)/Debye<<" , "<<(mza+mzb)/Debye<<" ). Absolute value, Debye = "<<mtot/Debye<<endl;
  out<<"Total nuclear dipole moment, Debye = ( "<<mu_nucl.x/Debye<<" , "<<mu_nucl.y/Debye<<" , "<<mu_nucl.z/Debye<<" ). Absolute value, Debye = "<<mu_nucl.length()/Debye<<endl;
  mtot = sqrt((mu_nucl.x+mxa+mxb)*(mu_nucl.x+mxa+mxb) + (mu_nucl.y+mya+myb)*(mu_nucl.y+mya+myb) + (mu_nucl.z+mza+mzb)*(mu_nucl.z+mza+mzb));
  out<<"Total dipole moment, Debye = ( "<<(mu_nucl.x+mxa+mxb)/Debye<<" , "<<(mu_nucl.y+mya+myb)/Debye<<" , "<<(mu_nucl.z+mza+mzb)/Debye<<" ). Absolute value, Debye = "<<mtot/Debye<<endl;




  out.close();

  delete tmp;
  delete pop;

}


void print_excitations(std::string filename, Electronic* el, Control_Parameters& prms){


  ofstream out(filename.c_str(),ios::out);

  vector<pair<double,double> > spectr;

  for(int e=0;e<prms.num_excitations;e++){

    int ii1 = el->Nocc_alp + prms.excitations[e].from_orbit[0] - 1; // index of the source MO
    int ii2 = el->Nocc_alp + prms.excitations[e].to_orbit[0] - 1;   // index of the target MO


    out<<"Excitation "<<e<<" : energy(a.u.)= "<<prms.excitations[e].Energy<<"  energy(eV)= "<<prms.excitations[e].Energy/eV<<"  "
       <<" orb_diff(eV)= "<<(el->bands_alp[ii2].second - el->bands_alp[ii1].second)/eV
       <<" oscillator_strength= "<<prms.excitations[e].f<<"  "
       <<"HOMO"<<((prms.excitations[e].from_orbit[0]>=0)?"+":"")<<prms.excitations[e].from_orbit[0]
       <<"( "<<ii1<<" ) "<<((prms.excitations[e].from_spin[0]==1)?"A":"B")<<" -> "
       <<"HOMO"<<((prms.excitations[e].to_orbit[0]>=0)?"+":"")<<prms.excitations[e].to_orbit[0]
       <<"( "<<ii2<<" ) "<<((prms.excitations[e].to_spin[0]==1)?"A":"B")<<"\n";


    spectr.push_back(pair<double,double>(prms.excitations[e].Energy/eV, prms.excitations[e].f));

  }// for e

  out.close();

  

  // Sort the spectr by the first value
  std:sort(spectr.begin(),spectr.end());

  ofstream out1((filename+"_sorted").c_str(),ios::out);

  double min_en = spectr[0].first;
  double max_en = spectr[prms.num_excitations-1].first;
  double de = 0.05; // in eV
  int npts = int((max_en - min_en + 5.0)/de) + 1;

  vector<double> intensity(npts,0.0);  // pure points
  vector<double> intensity1(npts,0.0); // convolved with Gaussian


  for(int i=0;i<spectr.size();i++){
  
    int indx = int(spectr[i].first / de);
    intensity[indx] += spectr[i].second;
 
  }


  // Convolution with Gaussian
  double var = prms.spectral_width;
  double prefac = 1.0/(var*sqrt(2.0*M_PI));
  double alp = 0.5/(var*var);

  for(int i=0;i<npts;i++){  
    double X = min_en + i*de;

    for(int j=0;j<npts;j++){
      double X0 = min_en + j*de;

      double w = prefac*exp(-alp*(X-X0)*(X-X0));  // normalized Gaussian
   
      intensity1[i] += w*intensity[j];

    }// for j 
  }// for i

    
  out1<<" min_en = "<<min_en<<"  max_en = "<<max_en<<endl;
  out1<<"  Energy(eV)         Oscillator strength(points)     Oscillator strength(convolved)\n";
  for(int i=0;i<npts;i++){
    out1<<(min_en + de * i)<<"  "<<intensity[i]<<"   "<<intensity1[i]<<endl;
  }
  

  out1.close();

//  exit(0);

}


