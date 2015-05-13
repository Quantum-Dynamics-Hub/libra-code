#include "SMATRIX.h"
#include "VECTOR.h"
#include <stdio.h>
#include <string.h>

namespace libmmath{
namespace liblinalg{


SMATRIX::SMATRIX(const SMATRIX& obj){

//  cout<<" In SMATRIX cc-tor\n";

  num_of_rows  = obj.num_of_rows;
  num_of_cols  = obj.num_of_cols;
  num_of_elems = obj.num_of_elems;

  indx1 = obj.indx1;
  indx2 = obj.indx2;
  M     = obj.M;

}
SMATRIX::~SMATRIX()
{

//  cout<<" In SMATRIX des-tor\n";

  indx1.clear();
  indx2.clear();
  M.clear();

  num_of_elems = 0;
}



SMATRIX SMATRIX::operator-(){

  SMATRIX tmp(*this);
  for(int i=0;i<tmp.num_of_elems;i++){
    tmp.M[i]=-tmp.M[i];
  }
  return tmp;

}

/*
  SMATRIX operator+(SMATRIX ob);
  SMATRIX operator-(SMATRIX ob);
  void operator+=(SMATRIX ob);
  void operator-=(SMATRIX ob);
  SMATRIX operator/(double num);
  SMATRIX operator=(SMATRIX ob);
  SMATRIX operator=(double num);
*/


SMATRIX SMATRIX::operator*(const SMATRIX& ob){

  SMATRIX tmp; //tmp.reserve(num_of_elems * ob.num_of_elems + 1);

  tmp.init(num_of_elems * ob.num_of_elems + 1);
  int kmax = 0;

  for(int n1=0;n1<num_of_elems;n1++){
    for(int n2=0;n2<ob.num_of_elems;n2++){
    
      if(indx2[n1]==ob.indx1[n2]){

        int k = tmp.find_element(indx1[n1],ob.indx2[n2]);        
        double x = M[n1]*ob.M[n2];
     

        if(fabs(x)>0.02){
        
          if(k==-1){ // New element
//            tmp.add_element(indx1[n1],ob.indx2[n2],x);

            tmp.indx1[kmax] = indx1[n1];
            tmp.indx2[kmax] = ob.indx2[n2];
            tmp.M[kmax] = x;
            kmax++;
            tmp.num_of_elems++;
          }
          else{  // add to existing element
            tmp.M[k] += x;
          }

        }// |x| > tol

      }// non-zero multiplication        

    }// for n2
  }// for n1

  return tmp;

}

/*
MATRIX MATRIX::operator +(MATRIX ob)
{   MATRIX Temp(num_of_cols,num_of_rows);
    for(int i=0;i<num_of_elems;i++) {Temp.M[i]=M[i]+ob.M[i];}
        return Temp;
}
MATRIX MATRIX::operator -(MATRIX ob)
{   MATRIX Temp(num_of_cols,num_of_rows);
    for(int i=0;i<num_of_elems;i++) {Temp.M[i]=M[i]-ob.M[i];}
        return Temp;
}

void MATRIX::operator+=(MATRIX ob)
{ for(int i=0;i<num_of_elems;i++) {M[i]+=ob.M[i];}
}
void MATRIX::operator-=(MATRIX ob)
{ for(int i=0;i<num_of_elems;i++) {M[i]-=ob.M[i];}
}

MATRIX MATRIX::operator /(double num)
{   MATRIX m(num_of_cols,num_of_rows);
    for(int i=0;i<num_of_elems;i++){  *(m.M+i)=*(M+i)/num;  }
        return m;
}
MATRIX MATRIX::operator =(MATRIX ob)
{  
     memcpy(M,ob.M,sizeof(double)*num_of_elems);

     return *this;
}
MATRIX MATRIX::operator=(double num){
    for(int i=0;i<num_of_elems;i++){
                M[i] = num;
        }
        return *this;
}

*/

}// namespace liblinalg
}// namespace libmmath

