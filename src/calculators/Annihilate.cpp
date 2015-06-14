#include "Annihilate.h"

namespace libcalculators{

void annihilate(int Na, int Nb, MATRIX* Pa, MATRIX* Pb, MATRIX* Ra, MATRIX* Rb){
// Spin annihilation
// according to Eqs. 2.19 - 2.21 from Pople, Beveridge, Dobosh JCP 47, 2026 (1967)

  int DF = 0;

  if(DF){   cout<<"in annihilate...\n"; }

  int Norb = Pa->num_of_cols;
  MATRIX* Pab;   Pab   = new MATRIX(Norb,Norb);
  MATRIX* Pba;   Pba   = new MATRIX(Norb,Norb);
  MATRIX* Pabab; Pabab = new MATRIX(Norb,Norb);
  MATRIX* Pbaba; Pbaba = new MATRIX(Norb,Norb);

  *Pab = *Pa * *Pb; 
  *Pba = *Pb * *Pa; 
  *Pabab = *Pab * *Pab; 
  *Pbaba = *Pba * *Pba; 

  double Tr_ab = Pab->tr();
  double Tr_abab = Pabab->tr();

  if(DF){
    cout<<"Tr_ab = "<<Tr_ab<<endl;
    cout<<"Tr_abab = "<<Tr_abab<<endl;
  }

  double p = Na; // number of alpha electrons
  double q = Nb; // number of beta electrons
  double s = 0.5*(p-q); 
  double A = q - 2.0*(s+1);                                                                                              
  double M2 = A*A - 2.0*A*Tr_ab + p*q - (p+q)*Tr_ab + 2.0*Tr_ab + 2.0*Tr_ab*Tr_ab - 2.0*Tr_abab;

  if(DF){
    cout<<"p = "<<p<<endl;
    cout<<"q = "<<q<<endl;
    cout<<"s = "<<s<<endl;
    cout<<"A = "<<A<<endl;
    cout<<"M2 = "<<M2<<endl;
  }
  

  // Alpha   
  *Ra =  ( M2 + Tr_ab - q ) * (*Pa)
       + ( p - Tr_ab ) * (*Pb)
       + ( *Pb * *Pab ) 
       + ((p+q) - 4.0*Tr_ab - 3.0 + 2.0*A) * (*Pa * *Pba)
       + ( 2.0*Tr_ab - p + 1.0 - A)*(*Pab + *Pba)
       - 2.0 * (*Pabab + *Pbaba)
       + 4.0 * (*Pabab * *Pa);

  *Ra = (1.0/M2) * (*Ra);

  if(DF){
    cout<<"Original alpha density matrix =\n"<<*Pa<<endl;
    cout<<"Annihilated alpha density matrix =\n"<<*Ra<<endl;
  }


  // Beta
  *Rb =  ( M2 + Tr_ab - p ) * (*Pb)
       + ( q - Tr_ab ) * (*Pa)
       + ( *Pa * *Pba ) 
       + ((p+q) - 4.0*Tr_ab - 3.0 + 2.0*A) * (*Pb * *Pab)
       + ( 2.0*Tr_ab - q + 1.0 - A)*(*Pba + *Pab)
       - 2.0 * (*Pbaba + *Pabab)
       + 4.0 * (*Pbaba * *Pb);

  *Rb = (1.0/M2) * (*Rb);

  if(DF){
    cout<<"Original beta density matrix =\n"<<*Pb<<endl;
    cout<<"Annihilated beta density matrix =\n"<<*Rb<<endl;
  }

  // Clean temporary objects
  delete Pab;
  delete Pba;
  delete Pabab;
  delete Pbaba;

}

void annihilate(int Na, int Nb, MATRIX* Pa, MATRIX* Pb){

  int Norb = Pa->num_of_cols;
  MATRIX* Ra;   Ra   = new MATRIX(Norb,Norb);
  MATRIX* Rb;   Rb   = new MATRIX(Norb,Norb);

  annihilate(Na, Nb, Pa, Pb, Ra, Rb);

  *Pa = *Ra;
  *Pb = *Rb;

  delete Ra;
  delete Rb;

}

}// namespace libcalculators
