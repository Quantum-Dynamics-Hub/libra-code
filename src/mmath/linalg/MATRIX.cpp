/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "MATRIX.h"
#include "VECTOR.h"
#include <stdio.h>
#include <string.h>
//using namespace std;

// ========================= Matrices ================================
// ------------------------- Constructors ----------------------------
using namespace libio;

namespace libmmath{
namespace liblinalg{


MATRIX::MATRIX(int m,int n)
{  MATRIX_PRECISION=8;
   MATRIX_WIDTH=15; 
   num_of_rows=m;
   num_of_cols=n;
  
   num_of_elems=num_of_rows*num_of_cols;
  
   M=new double[num_of_elems];
   for(int i=0;i<num_of_elems;i++)      {   *(M+i)=0.0; }
} 

MATRIX::MATRIX(const VECTOR& u1, const VECTOR& u2, const VECTOR& u3){
   MATRIX_PRECISION=8;
   MATRIX_WIDTH=15;
   num_of_rows=3;
   num_of_cols=3;
    
   num_of_elems=9;

   M=new double[num_of_elems];
   for(int i=0;i<num_of_elems;i++)      {   *(M+i)=0.0; }
  
   M[0] = u1.x;  M[1] = u2.x;  M[2] = u3.x;
   M[3] = u1.y;  M[4] = u2.y;  M[5] = u3.y;
   M[6] = u1.z;  M[7] = u2.z;  M[8] = u3.z;

} 

MATRIX::MATRIX(const MATRIX& obj){

  num_of_rows  = obj.num_of_rows;
  num_of_cols  = obj.num_of_cols;
  num_of_elems = obj.num_of_elems;

  M = new double[num_of_elems];

  for(int i=0;i<num_of_elems;i++){

          M[i] = obj.M[i];
  }

  MATRIX_PRECISION = obj.MATRIX_PRECISION;
  MATRIX_WIDTH     = obj.MATRIX_WIDTH;

  NonOrtMeasure  = obj.NonOrtMeasure;

}
MATRIX::~MATRIX()
{
//        cout<<"Matrix destructor! "<<this<<endl;
   delete [] M;
}

// --------------------------- Initialization -----------------------------


void MATRIX::Init(double x)
{   for(int i=0;i<num_of_elems;i++)     {   *(M+i)=x;   }
}
void MATRIX::InitSquareMatrix(int dim,double x)
{   num_of_rows = num_of_cols = dim;
        num_of_elems=num_of_rows*num_of_cols;

        M = new double[num_of_elems];
        for(int i=0;i<num_of_elems;i++) {   *(M+i)=x;   }
}

bool MATRIX::Init_Unit_Matrix(double x){
  bool res;

  if(num_of_rows!=num_of_cols){ res=false;   }
  else{
    res=true;
    int k=0;
    for(int i=0;i<num_of_cols;i++){
      for(int j=0;j<num_of_cols;j++){
        if(i==j) {*(M+k)=x;  }
        else     {*(M+k)=0.0;}
        k++;
      }
    }
  }
  return res;

}

bool MATRIX::Load_Matrix_From_File(char *FileName)
{  bool res;
   num_of_elems=num_of_rows*num_of_cols;
   ifstream ob;

   ob.open(FileName);
   if(ob.is_open())
   {
           for(int k=0;k<num_of_elems;k++)
           {   ob>>*(M+k);
           }

   ob.close();
   res=true;
   }
   else res=false;

   return res;
}


void MATRIX::init(VECTOR& u1, VECTOR& u2, VECTOR& u3){

   M[0] = u1.x;  M[1] = u2.x;  M[2] = u3.x;
   M[3] = u1.y;  M[4] = u2.y;  M[5] = u3.y;
   M[6] = u1.z;  M[7] = u2.z;  M[8] = u3.z;

}

void MATRIX::init(const VECTOR& u1, const VECTOR& u2, const VECTOR& u3){

   M[0] = u1.x;  M[1] = u2.x;  M[2] = u3.x;
   M[3] = u1.y;  M[4] = u2.y;  M[5] = u3.y;
   M[6] = u1.z;  M[7] = u2.z;  M[8] = u3.z;

}

void MATRIX::init(MATRIX& m){

  num_of_rows  = m.num_of_rows;
  num_of_cols  = m.num_of_cols;
  num_of_elems = m.num_of_elems;

  M = new double[num_of_elems];

  for(int i=0;i<num_of_elems;i++){

          M[i] = m.M[i];
  }

  MATRIX_PRECISION = m.MATRIX_PRECISION;
  MATRIX_WIDTH     = m.MATRIX_WIDTH;

  NonOrtMeasure  = m.NonOrtMeasure;


}

void MATRIX::init(const MATRIX& m){

  num_of_rows  = m.num_of_rows;
  num_of_cols  = m.num_of_cols;
  num_of_elems = m.num_of_elems;

  M = new double[num_of_elems];

  for(int i=0;i<num_of_elems;i++){

          M[i] = m.M[i];
  }

  MATRIX_PRECISION = m.MATRIX_PRECISION;
  MATRIX_WIDTH     = m.MATRIX_WIDTH;

  NonOrtMeasure  = m.NonOrtMeasure;


}

//------------------------ Tehniques for ortogonalization and decomposition --------------------

double MATRIX::NonOrtogonality_Measure()
{  double sum=0.0,aa=0.0;
    for(int i=0;i<num_of_cols-1;i++)
    {   for(int j=i+1;j<num_of_cols;j++)
        {   aa=0.0;
            for(int k=0;k<num_of_rows;k++)
            {   aa+=(M[k*num_of_cols+i]*M[k*num_of_cols+j]);
            }
            sum+=aa*aa;
        }
    }
    NonOrtMeasure=sum;
    return sum;
}
void MATRIX::Show_Orthogonalization_Matrix()
{  int k=0;
   cout<<endl;
   for(int row=0;row<num_of_rows;row++)
   {   for(int col=0;col<num_of_rows;col++)
       {   cout.width(MATRIX_WIDTH);
           cout.precision(MATRIX_PRECISION);
           cout<<*(Orthogonalization_Matrix+k);
           k++;
       }
       cout<<endl;
   }
}

void MATRIX::Show_Orthogonalization_Matrix(char * Output_File)
{  ofstream ob;
   ob.open(Output_File);
   if(ob.is_open())
   {
   int k=0;
   ob<<endl;
   for(int row=0;row<num_of_rows;row++)
   {   for(int col=0;col<num_of_rows;col++)
       {   ob.width(MATRIX_WIDTH);
           ob.precision(MATRIX_PRECISION);
           ob<<*(Orthogonalization_Matrix+k);
           k++;
       }
       ob<<endl;
   }
   }
   else;
}
void MATRIX::Delete_Orthogonalization_Matrix()
{   delete [] Orthogonalization_Matrix;
}
//------------------------------------------------------------------------------------------------------------
/*
     Óìíîæåíèå ñïðàâà íà ìàòðèöó ïîâîðîòà â ïëîñêîñòè (i,j), â êîòîðîé
          i-ûé è j-ûé ñòîëáåö è ñòðîêà èìåþò âèä
           
                        i                 | cos\f       sin\f |
            j         |-sin\f   cos\f |
                        
                                                  i       j

*/
void MATRIX::RightRotation(int i,int j,double sine,double cosine)
{   double a_row_i,a_row_j, A_row_i,A_row_j;
        int k_i,k_j;
    for(int row=0;row<num_of_rows;row++)
        {   k_i=row*num_of_cols+i;
            k_j=row*num_of_cols+j;
                a_row_i=*(M+k_i);
            a_row_j=*(M+k_j);

                A_row_i= a_row_i*cosine + a_row_j*sine;
                A_row_j=-a_row_i*sine   + a_row_j*cosine;

                *(M+k_i)=A_row_i;
                *(M+k_j)=A_row_j;
        }
}

/*
     Óìíîæåíèå ñëåâà íà ìàòðèöó ïîâîðîòà â ïëîñêîñòè (i,j), â êîòîðîé
          i-ûé è j-ûé ñòîëáåö è ñòðîêà èìåþò âèä
           
                        i                 | cos\f       -sin\f |
            j         | sin\f    cos\f |
                        
                                                  i       j

*/
void MATRIX::LeftRotation(int i,int j,double sine,double cosine)
{   double a_i_col,a_j_col, A_i_col,A_j_col;
        int k_i,k_j;
    for(int col=0;col<num_of_cols;col++)
        {   k_i=i*num_of_cols+col;
            k_j=j*num_of_cols+col;
                a_i_col=*(M+k_i);
            a_j_col=*(M+k_j);

                A_i_col= a_i_col*cosine + a_j_col*sine;
                A_j_col=-a_i_col*sine   + a_j_col*cosine;

                *(M+k_i)=A_i_col;
                *(M+k_j)=A_j_col;
        }

}

void MATRIX::Ortogonalization(double Eps)
{  double q,p,v,sine,cosine,Z=0.0,aa;
   Orthogonalization_Matrix=new double[num_of_cols*num_of_cols];
   MATRIX V(num_of_cols,num_of_cols), Vk(num_of_cols,num_of_cols);
          V.Init_Unit_Matrix(1.0);


   for(int i=0;i<num_of_cols-1;i++)
        {   for(int j=i+1;j<num_of_cols;j++)
            {      aa=0.0;
                for(int k=0;k<num_of_rows;k++)
                {   aa+=(M[k*num_of_cols+i]*M[k*num_of_cols+j]);
                }
                    Z+=aa*aa;
            }
        }

   while(Z>Eps)
   {

   for(int i=0;i<(num_of_cols-1);i++)
   {   for(int j=i+1;j<num_of_cols;j++)
       {      p=0.0;
              q=0.0;
           for(int k=0;k<num_of_rows;k++)
              {   p+=M[k*num_of_cols+i]*M[k*num_of_cols+j];
                  q+=(M[k*num_of_cols+i]*M[k*num_of_cols+i]-M[k*num_of_cols+j]*M[k*num_of_cols+j]);
              }

                   v=sqrt(4*p*p+q*q);

           if(q>=0)
           {  cosine=sqrt(0.5*((q/v)+1));
              sine=p/(v*cosine);
           }
           else
           {  if(p>=0) { sine= sqrt(0.5*(1.0-q/v)); }
              if(p<0)  { sine=-sqrt(0.5*(1.0-q/v)); }
                         cosine=p/(v*sine);
           }

           //-------------------------------------------
           Vk.Init_Unit_Matrix(1.0);

           Vk.set(i,i,cosine);
           Vk.set(j,j,cosine);
           Vk.set(i,j,-sine);
           Vk.set(j,i, sine);

           V=V*Vk;
           //-----------------------------------------
           // RightRotation-procedure.
            double a_row_i,a_row_j, A_row_i,A_row_j;
            int k_i,k_j;
        for(int row=0;row<num_of_rows;row++)
           {   k_i=row*num_of_cols+i;
               k_j=row*num_of_cols+j;
                a_row_i=*(M+k_i);
                a_row_j=*(M+k_j);

                A_row_i= a_row_i*cosine + a_row_j*sine;
                A_row_j=-a_row_i*sine   + a_row_j*cosine;

                *(M+k_i)=A_row_i;
                *(M+k_j)=A_row_j;
            }
           //------------------------------------------
           // NonOrthogonalityMeasure
            Z=0.0;
        for(int i1=0;i1<num_of_cols-1;i1++)
        {   for(int j1=i1+1;j1<num_of_cols;j1++)
            {   aa=0.0;
                for(int k=0;k<num_of_rows;k++)
                {   aa+=(M[k*num_of_cols+i1]*M[k*num_of_cols+j1]);
                }
                Z+=aa*aa;
            }
        }
           //------------------------------------------
           // cout<<Z<<endl;
       }
   }

   }
   int k=0;
   for(int row=0;row<num_of_rows;row++)
   {   for(int col=0;col<num_of_rows;col++)
       {   Orthogonalization_Matrix[k]=V.get(row,col);
           k++;
       }
   }
   V.Delete_Matrix();
   Vk.Delete_Matrix();
}

void MATRIX::RL_Decomposition(MATRIX *L,MATRIX *R)
{  double *R_time;   R_time=new double[num_of_rows*num_of_cols];
   double *L_time;   L_time=new double[num_of_rows*num_of_cols];

   // Âî âðåìåííûé ìàññèâ êîïèðóåì ñîäåðæèìîå äàííîé ìàòðèöû, ÷òîáû åå íå ïîâðåäèòü.
   for(int k=0;k<num_of_rows*num_of_cols;k++)
   {   R_time[k]=M[k];
   }
   // Ñîçäàåì åäèíè÷íöþ ìàòðèöó òîé æå ðàçìåðíîñòè, ÷òî è äàííàÿ ìàòðèöà.
   int k=0;
   for(int i=0;i<num_of_cols;i++)
   {   for(int j=0;j<num_of_cols;j++)
       {    if(i==j) {*(L_time+k)=1.0;  }
            else     {*(L_time+k)=0.0;}
            k++;
       }
   }

 // Îñíîâíàÿ ÷àñòü àëãîðèòìà ðàçëîæåíèÿ êâàäðàòíîé ìàòðèöû íà òðåóãîëüíûå ñîìíîæèòåëè.
   double alpha;
   for(int row1=0;row1<num_of_rows-1;row1++)
   {   if(R_time[row1*num_of_cols+row1]!=0)
       {
            for(int row2=row1+1;row2<num_of_rows;row2++)
            {   if(R_time[row2*num_of_cols+row1]!=0)
                {
                   alpha=-R_time[row2*num_of_cols+row1]/R_time[row1*num_of_cols+row1];

                    for(int col=0;col<num_of_cols;col++)
                    {   R_time[row2*num_of_cols+col]=R_time[row2*num_of_cols+col]+alpha*R_time[row1*num_of_cols+col];
                        L_time[row2*num_of_cols+col]=L_time[row2*num_of_cols+col]+alpha*L_time[row1*num_of_cols+col];
                    }
                }
                else continue;
            }
       }
   }
   k=0;
   for(int row=0;row<num_of_rows;row++)
   {   for(int col=0;col<num_of_cols;col++)
       {   L->set(row,col,L_time[k]);
           R->set(row,col,R_time[k]);
           k++;
       }
   }
   delete [] R_time;
   delete [] L_time;
}

// -------------------------------- Basic matrix operations ----------------------------------

void MATRIX::Transpose()
{  int old_num_of_rows=num_of_rows;
   int old_num_of_cols=num_of_cols;

   int new_num_of_rows=num_of_cols;
   int new_num_of_cols=num_of_rows;

   int N=num_of_rows*num_of_cols;
   double *Temp;
           Temp=new double[N];

   for(int row=0;row<old_num_of_rows;row++)
   {   for(int col=0;col<old_num_of_cols;col++)
       {   Temp[col*new_num_of_cols+row]=M[row*old_num_of_cols+col];
       }
   }
   for(int k=0;k<N;k++)
   {  M[k]=Temp[k];
   }
   delete [] Temp;

   num_of_rows=new_num_of_rows;
   num_of_cols=new_num_of_cols;
}
MATRIX MATRIX::T()
{

   MATRIX res(num_of_cols,num_of_rows);


   for(int row=0;row<num_of_rows;row++)
   {   for(int col=0;col<num_of_cols;col++)
       {   res.M[col*num_of_rows+row]=M[row*num_of_cols+col];
       }
   }

   return res;
}

void MATRIX::Inverse(MATRIX *INV)
{    double *R_time;   R_time=new double[num_of_rows*num_of_cols];
     double *L_time;   L_time=new double[num_of_rows*num_of_cols];

   // Âî âðåìåííûé ìàññèâ êîïèðóåì ñîäåðæèìîå äàííîé ìàòðèöû, ÷òîáû åå íå ïîâðåäèòü.
   for(int k=0;k<num_of_rows*num_of_cols;k++)
   {   R_time[k]=M[k];
   }
   // Ñîçäàåì åäèíè÷íöþ ìàòðèöó òîé æå ðàçìåðíîñòè, ÷òî è äàííàÿ ìàòðèöà.
   int k=0;
   for(int i=0;i<num_of_cols;i++)
   {   for(int j=0;j<num_of_cols;j++)
       {    if(i==j) {*(L_time+k)=1.0;  }
            else     {*(L_time+k)=0.0;}
            k++;
       }
   }

 // Ïðÿìîé õîä àëãîðèòìà.
   double alpha;
   for(int row1=0;row1<num_of_rows-1;row1++)
   {
       if(R_time[row1*num_of_cols+row1]==0)
       {   int row=row1+1;
           while(R_time[row*num_of_cols+row1]==0)
           {   row++;
           }
           double temp1,temp2;
           for(int col=0;col<num_of_cols;col++)
           {   temp1=R_time[row1*num_of_cols+col];
               R_time[row1*num_of_cols+col]=R_time[row*num_of_cols+col];
               R_time[row*num_of_cols+col]=temp1;


               temp2=L_time[row1*num_of_cols+col];
               L_time[row1*num_of_cols+col]=L_time[row*num_of_cols+col];
               L_time[row*num_of_cols+col]=temp2;

           }
       }


       if(R_time[row1*num_of_cols+row1]!=0)
       {
            for(int row2=row1+1;row2<num_of_rows;row2++)
            {   if(R_time[row2*num_of_cols+row1]!=0)
                {
                   alpha=-R_time[row2*num_of_cols+row1]/R_time[row1*num_of_cols+row1];

                    for(int col=0;col<num_of_cols;col++)
                    {   R_time[row2*num_of_cols+col]=R_time[row2*num_of_cols+col]+alpha*R_time[row1*num_of_cols+col];
                        L_time[row2*num_of_cols+col]=L_time[row2*num_of_cols+col]+alpha*L_time[row1*num_of_cols+col];
                    }
                }
                else continue;
            }

       }

   }
 // Îáðàòíûé õîä àëãîðèòìà.
     for(int row1=num_of_rows-1;row1>0;row1--)
     {
                  alpha=R_time[row1*num_of_cols+row1];
                        R_time[row1*num_of_cols+row1]=1.0;
         for(int col=(num_of_cols-1);col>=0;col--)
         {   L_time[row1*num_of_cols+col]=L_time[row1*num_of_cols+col]/alpha;
         }
         for(int row2=row1-1;row2>=0;row2--)
         {        alpha=-R_time[row2*num_of_cols+row1];
              for(int col=(num_of_cols-1);col>=0;col--)
              {   R_time[row2*num_of_cols+col]=R_time[row2*num_of_cols+col]+alpha*R_time[row1*num_of_cols+col];
                  L_time[row2*num_of_cols+col]=L_time[row2*num_of_cols+col]+alpha*L_time[row1*num_of_cols+col];
              }
         }
     }
                  alpha=R_time[0];
                        R_time[0]=1.0;
         for(int col=(num_of_cols-1);col>=0;col--)
         {   L_time[col]=L_time[col]/alpha;
         }

 //
   k=0;
   for(int row=0;row<num_of_rows;row++)
   {   for(int col=0;col<num_of_cols;col++)
       {   INV->set(row,col,L_time[k]);
           //R->Set_Element(row,col,R_time[k]);
           k++;
       }
   }
   delete [] R_time;
   delete [] L_time;
}

void MATRIX::Inverse(MATRIX& INV)
{    double *R_time;   R_time=new double[num_of_rows*num_of_cols];
     double *L_time;   L_time=new double[num_of_rows*num_of_cols];

   // Âî âðåìåííûé ìàññèâ êîïèðóåì ñîäåðæèìîå äàííîé ìàòðèöû, ÷òîáû åå íå ïîâðåäèòü.
   for(int k=0;k<num_of_rows*num_of_cols;k++)
   {   R_time[k]=M[k];
   }
   // Ñîçäàåì åäèíè÷íöþ ìàòðèöó òîé æå ðàçìåðíîñòè, ÷òî è äàííàÿ ìàòðèöà.
   int k=0;
   for(int i=0;i<num_of_cols;i++)
   {   for(int j=0;j<num_of_cols;j++)
       {    if(i==j) {*(L_time+k)=1.0;  }
            else     {*(L_time+k)=0.0;}
            k++;
       }
   }

 // Ïðÿìîé õîä àëãîðèòìà.
   double alpha;
   for(int row1=0;row1<num_of_rows-1;row1++)
   {
       if(R_time[row1*num_of_cols+row1]==0)
       {   int row=row1+1;
           while(R_time[row*num_of_cols+row1]==0)
           {   row++;
           }
           double temp1,temp2;
           for(int col=0;col<num_of_cols;col++)
           {   temp1=R_time[row1*num_of_cols+col];
               R_time[row1*num_of_cols+col]=R_time[row*num_of_cols+col];
               R_time[row*num_of_cols+col]=temp1;


               temp2=L_time[row1*num_of_cols+col];
               L_time[row1*num_of_cols+col]=L_time[row*num_of_cols+col];
               L_time[row*num_of_cols+col]=temp2;

           }
       }


       if(R_time[row1*num_of_cols+row1]!=0)
       {
            for(int row2=row1+1;row2<num_of_rows;row2++)
            {   if(R_time[row2*num_of_cols+row1]!=0)
                {
                   alpha=-R_time[row2*num_of_cols+row1]/R_time[row1*num_of_cols+row1];

                    for(int col=0;col<num_of_cols;col++)
                    {   R_time[row2*num_of_cols+col]=R_time[row2*num_of_cols+col]+alpha*R_time[row1*num_of_cols+col];
                        L_time[row2*num_of_cols+col]=L_time[row2*num_of_cols+col]+alpha*L_time[row1*num_of_cols+col];
                    }
                }
                else continue;
            }

       }

   }
 // Îáðàòíûé õîä àëãîðèòìà.
     for(int row1=num_of_rows-1;row1>0;row1--)
     {
                  alpha=R_time[row1*num_of_cols+row1];
                        R_time[row1*num_of_cols+row1]=1.0;
         for(int col=(num_of_cols-1);col>=0;col--)
         {   L_time[row1*num_of_cols+col]=L_time[row1*num_of_cols+col]/alpha;
         }
         for(int row2=row1-1;row2>=0;row2--)
         {        alpha=-R_time[row2*num_of_cols+row1];
              for(int col=(num_of_cols-1);col>=0;col--)
              {   R_time[row2*num_of_cols+col]=R_time[row2*num_of_cols+col]+alpha*R_time[row1*num_of_cols+col];
                  L_time[row2*num_of_cols+col]=L_time[row2*num_of_cols+col]+alpha*L_time[row1*num_of_cols+col];
              }
         }
     }
                  alpha=R_time[0];
                        R_time[0]=1.0;
         for(int col=(num_of_cols-1);col>=0;col--)
         {   L_time[col]=L_time[col]/alpha;
         }

 //
   k=0;
   for(int row=0;row<num_of_rows;row++)
   {   for(int col=0;col<num_of_cols;col++)
       {   INV.set(row,col,L_time[k]);
           //R->Set_Element(row,col,R_time[k]);
           k++;
       }
   }
   delete [] R_time;
   delete [] L_time;
}
double& MATRIX::operator[](const int indx){

    return M[indx];
}


MATRIX MATRIX::operator-(){
  MATRIX tmp(num_of_rows,num_of_cols);

  for(int i=0;i<num_of_elems;i++){
    tmp.M[i] = -M[i];
  }

  return tmp;
}
MATRIX MATRIX::operator*(const MATRIX& ob){
  
  int n=ob.num_of_cols; //show_num_of_cols();

  MATRIX Temp(num_of_rows,n);
  //  Âàæíî ÷òîáû num_of_cols==m=ob.show_num_of_rows(); !!! 
  //  Ýòî óñëîâèå çäåñü îäíàêî íå ïðîâåðÿåòñÿ 
  //      è ïîýòîìó äîëæíî ïðîâåðÿòüñÿ ïåðåä âûïîëíåíèåì ýòîé îïåðàöèè !!!

  for(int row=0;row<num_of_rows;row++){
    for(int col=0;col<n;col++){

      double d=0.0;
      for(int k=0;k<num_of_cols;k++){
        d+=M[k+row*num_of_cols]*ob.M[k*n+col];
      }// for k
      Temp.M[row*n+col] = d;

    }// for col
  }// for row

  return Temp;
}

MATRIX MATRIX::operator+(const MATRIX& ob){
  MATRIX Temp(num_of_cols,num_of_rows);
  for(int i=0;i<num_of_elems;i++){ 
    Temp.M[i]=M[i]+ob.M[i];
  }

  return Temp;
}

MATRIX MATRIX::operator-(const MATRIX& ob){
  MATRIX Temp(num_of_cols,num_of_rows);
  for(int i=0;i<num_of_elems;i++){
    Temp.M[i]=M[i]-ob.M[i];
  }

  return Temp;
}

MATRIX& MATRIX::operator+=(const MATRIX& ob){ 
  for(int i=0;i<num_of_elems;i++) { 
    M[i] += ob.M[i];
  }
  return *this; // return reference to allow chaining: (((A += 5) += 6) += 7)  etc..
}
                           
MATRIX& MATRIX::operator-=(const MATRIX& ob){
 
  for(int i=0;i<num_of_elems;i++){ 
    M[i] -= ob.M[i];
  }
  return *this; // return reference to allow chaining: (((A -= 5) -= 6) -= 7)  etc..
}

MATRIX& MATRIX::operator*=(double f){ 
  for(int i=0;i<num_of_elems;i++) { 
    M[i] *= f;
  }
  return *this; // return reference to allow chaining: (((A *= 5) *= 6) *= 7)  etc..
}


MATRIX MATRIX::operator /(double num)
{   MATRIX m(num_of_cols,num_of_rows);
    for(int i=0;i<num_of_elems;i++){  *(m.M+i)=*(M+i)/num;  }
        return m;
}

MATRIX MATRIX::operator=(const MATRIX& ob){  

  memcpy(M,ob.M,sizeof(double)*num_of_elems);

  return *this; // return reference to allow chaining: A = B = C =... !!! No: so crap doesn't happen in PYthon
}

MATRIX MATRIX::operator=(double num){
  for(int i=0;i<num_of_elems;i++){
    M[i] = num;
  }
  return *this;  // return reference to allow chaining: A = B = C = 7!!! No: so crap doesn't happen in PYthon
}


int operator ==(const MATRIX& m1,const MATRIX& m2){
        int res=1;
        for(int i=0;i<m1.num_of_elems;i++) {if(m1.M[i]!=m2.M[i]) {res=0;} else;}
        return res;
}
int operator !=(const MATRIX& m1,const MATRIX& m2){
        int res=0;
        for(int i=0;i<m1.num_of_elems;i++) {if(m1.M[i]!=m2.M[i]) {res=1;} else;}
        return res;
}
double operator %(MATRIX& m1,MATRIX& m2)
{    double res=0.0;
     for(int i=0;i<m1.num_of_elems;i++)  {   res+=*(m1.M+i)*(*(m2.M+i));          }
     return res;
}
MATRIX operator ^(VECTOR& v1,VECTOR& v2){
    MATRIX res(3,3);

    res.M[0] = v1.x*v2.x;   res.M[1] = v1.x*v2.y;  res.M[2] = v1.x*v2.z;
    res.M[3] = v1.y*v2.x;   res.M[4] = v1.y*v2.y;  res.M[5] = v1.y*v2.z;
    res.M[6] = v1.z*v2.x;   res.M[7] = v1.z*v2.y;  res.M[8] = v1.z*v2.z;

    return res;
}
MATRIX operator*(const double& f,  const MATRIX& m1){
        MATRIX m(m1.num_of_cols,m1.num_of_rows);
    for(int i=0;i<m1.num_of_elems;i++){  *(m.M+i)=*(m1.M+i)*f;  }
        return m;
}
MATRIX operator*(const MATRIX &m1, const double  &f){
        MATRIX m(m1.num_of_cols,m1.num_of_rows);
    for(int i=0;i<m1.num_of_elems;i++){  *(m.M+i)=*(m1.M+i)*f;  }
        return m;
}

// ------------------------------- Input & output functions ------------------------------

void MATRIX::show_matrix()
{  int k=0;
   cout<<endl;
   for(int row=0;row<num_of_rows;row++)
   {   for(int col=0;col<num_of_cols;col++)
       {   cout.width(MATRIX_WIDTH);
           cout.precision(MATRIX_PRECISION);
           cout<<*(M+k);
           k++;
       }
       cout<<endl;
   }
}

void MATRIX::show_matrix(char * Output_File)
{  ofstream ob;
   ob.open(Output_File);
   if(ob.is_open())
   {
   int k=0;
   ob<<endl;
   for(int row=0;row<num_of_rows;row++)
   {   for(int col=0;col<num_of_cols;col++)
       {   ob.width(MATRIX_WIDTH);
           ob.precision(MATRIX_PRECISION);
           ob<<*(M+k);
           k++;
       }
       ob<<endl;
   }
   }
   else;
}
ostream& operator<<(ostream &strm,MATRIX ob){
  strm.setf(ios::showpoint);
  for(int i=0;i<ob.num_of_rows;i++){
    for(int j=0;j<ob.num_of_cols;j++){

      strm.precision(ob.MATRIX_PRECISION);
      strm.width(ob.MATRIX_WIDTH);
      strm<<left;//right;
      strm<<ob.M[i*ob.num_of_cols+j]<<"  ";

    }// for j
    strm<<endl;
  }// for i
  return strm;

}

istream& operator>>(istream& strm,MATRIX &ob){
//     Do not defined for general case       !!!      
/*
       double val;
       string str;
       strm>>ob.x;//=val;
       strm>>ob.y;//=val;
       strm>>ob.z;//=val;
 */   return strm;
}

void MATRIX::Delete_Matrix(){  delete [] M;  }

void MATRIX::get_vectors(VECTOR& u1,VECTOR& u2,VECTOR& u3){

  u1.x = M[0];  u2.x = M[1];  u3.x = M[2];
  u1.y = M[3];  u2.y = M[4];  u3.y = M[5];
  u1.z = M[6];  u2.z = M[7];  u3.z = M[8];

}
void MATRIX::skew(VECTOR v) {

  M[0] = 0.0; M[1] =-v.z; M[2] = v.y;
  M[3] = v.z; M[4] = 0.0; M[5] =-v.x;
  M[6] =-v.y; M[7] = v.x; M[8] = 0.0;

}
void MATRIX::skew1(VECTOR v) {

  M[0] = 0.0; M[1] =-v.x; M[2] =-v.y; M[3] =-v.z;
  M[4] = v.x; M[5] = 0.0; M[6] = v.z; M[7] =-v.y;
  M[8] = v.y; M[9] =-v.z; M[10]= 0.0; M[11]= v.x;
  M[12]= v.z; M[13]= v.y; M[14]=-v.x; M[15]= 0.0;

}

void MATRIX::exp(MATRIX &m) {
    for(int i=0;i<num_of_elems;i++){
                M[i] = 0.0;
        }
        for(int i=0;i<num_of_rows;i++){ M[i * num_of_cols + i] = 1.0;}

    MATRIX &a=*this;
    MATRIX mul(num_of_rows,num_of_cols);
    mul = m;
        //cout<<(this)<<endl;
        //cout<<mul.M[0]<<"  "<<mul.M[1]<<"  "<<mul.M[2]<<endl;
        //cout<<mul.M[3]<<"  "<<mul.M[4]<<"  "<<mul.M[5]<<endl;
        //cout<<mul.M[6]<<"  "<<mul.M[7]<<"  "<<mul.M[8]<<endl;
    double sum;
    int j,i;
    for (i = 1; i < 40; i++) {
        mul = mul/double(i);
        a = a + mul;
        mul = mul * m;
        sum = 0.0;
        for (j = 0; j < num_of_elems; j++) {
            sum += fabs(mul.M[j]);
        }
        if (sum < 1e-16) break;
    }
    if (i > 35) {
        cout << "ERROR: exponent didn't converge" << endl;
                cout << "sum = "<<sum;
        exit(0);
    }
}
void MATRIX::exp(const MATRIX &m) {
    for(int i=0;i<num_of_elems;i++){
                M[i] = 0.0;
        }
        for(int i=0;i<num_of_rows;i++){ M[i * num_of_cols + i] = 1.0;}

    MATRIX &a=*this;
    MATRIX mul(num_of_rows,num_of_cols);
    mul = m;
        //cout<<(this)<<endl;
        //cout<<mul.M[0]<<"  "<<mul.M[1]<<"  "<<mul.M[2]<<endl;
        //cout<<mul.M[3]<<"  "<<mul.M[4]<<"  "<<mul.M[5]<<endl;
        //cout<<mul.M[6]<<"  "<<mul.M[7]<<"  "<<mul.M[8]<<endl;
    double sum;
    int i,j;
    for (i = 1; i < 40; i++) {
        mul = mul/double(i);
        a = a + mul;
        mul = mul * m;
        sum = 0.0;
        for (j = 0; j < num_of_elems; j++) {
            sum += fabs(mul.M[j]);
        }
        if (sum < 1e-16) break;
    }
    if (i > 35) {
        cout << "ERROR: exponent didn't converge" << endl;
                cout << "sum = "<<sum;
        exit(0);
    }
}
void MATRIX::FindMaxNondiagonalElement(int& row,int& col,double& value){
    int k=0;
    double elem, max_elem;
    value = M[1]; max_elem = fabs(value); row = 0 ; col = 1;

    for(int rw=0;rw<num_of_rows;rw++){
        for(int cl=rw+1;cl<num_of_cols;cl++){
            k = rw*num_of_cols + cl;
            elem = fabs(M[k]);
            if(elem>max_elem) {max_elem = elem; value = M[k]; col = cl; row = rw;}
        }
    }

}
void MATRIX::JACOBY_EIGEN(MATRIX& EVAL, MATRIX& EVECT) {

        // Description: only for symmetric matrices!
        // EVECT * EVAL  =  M * EVECT

        MATRIX V(num_of_rows,num_of_cols);
        MATRIX temp(num_of_rows,num_of_cols);

               temp.MATRIX_PRECISION = MATRIX_PRECISION;
                   for(int i=0;i<num_of_elems;i++){

                           temp.M[i] = M[i];
                   }

                   EVECT.Init_Unit_Matrix(1.0);


                   double EPS=1.0e-30;
                   double eps; eps=0.0;

                                  // Define convergence criteria.
                                  int k=0;
                                  for(int i=0;i<temp.num_of_rows;i++){
                                          for(int j=0;j<temp.num_of_cols;j++){
                                                  if(i!=j) {eps+=temp.M[k]*temp.M[k];} k++;
                                          }
                                  }

                                      int row,col;
                                          double  val;
                                          double phi;

                                  while(eps>EPS){

                                          temp.FindMaxNondiagonalElement(row,col,val);

                                          phi = 0.5*atan(2.0*val/(temp.M[col*(num_of_cols+1)]-temp.M[row*(num_of_cols+1)]));

                                          V.Init_Unit_Matrix(1.0);

                                          V.M[row*(num_of_cols+1)]   = cos(phi);   V.M[row*num_of_cols + col] = -sin(phi);
                                          V.M[col*num_of_cols + row] = sin(phi);   V.M[col*(num_of_cols+1)]   =  cos(phi);

                                          temp = V*temp;
                          EVECT = V*EVECT;
                                          V.Transpose();

                                          temp = temp*V;

                      eps = 0.0;
                                          int k=0;
                                               for(int i=0;i<temp.num_of_rows;i++){
                                                       for(int j=0;j<temp.num_of_cols;j++){
                                                               if(i!=j) {eps+=temp.M[k]*temp.M[k];} k++;
                                                                   }
                                                           }

                                  }
                                        EVAL = temp;
                                    EVECT.Transpose();
}
void MATRIX::JACOBY_EIGEN(MATRIX& EVAL, MATRIX &EVECT,double EPS) {

        // Description: only for symmetric matrices!
        // EVECT * EVAL  =  M * EVECT

// V = P^T
// EVAL = V_M V_{M-1} ... V_0 * M * V_0^T * V_1^T ... V_M^T = Q^T * M * Q
// EVECT = Q = V_0^T * V_1^T ... V_M^T

       MATRIX V(num_of_rows,num_of_cols);
       MATRIX temp(num_of_rows,num_of_cols);
       int row, col, i, j, k, num_iter;
       double val,phi,eps;


       temp.MATRIX_PRECISION = MATRIX_PRECISION;
       for(i=0;i<num_of_elems;i++){     temp.M[i] = M[i];   }

       EVECT = 0.0;
       EVECT.Init_Unit_Matrix(1.0);

       num_iter = 0;

       // Define convergence criteria.
       k=0; eps = 0.0;
       for(i=0;i<temp.num_of_rows;i++){
           for(j=0;j<temp.num_of_cols;j++){
               if(i!=j) {eps+=temp.M[k]*temp.M[k];}
               k++;
           }
       }

        while(eps>EPS){

            num_iter++;

            temp.FindMaxNondiagonalElement(row,col,val);

            phi = 0.5*atan(2.0*val/(temp.M[col*(num_of_cols+1)]-temp.M[row*(num_of_cols+1)]));

            V = 0.0;
            V.Init_Unit_Matrix(1.0);

            V.M[row*(num_of_cols+1)]   = cos(phi);   V.M[row*num_of_cols + col] = -sin(phi);
            V.M[col*num_of_cols + row] = sin(phi);   V.M[col*(num_of_cols+1)]   =  cos(phi);

            temp = V*temp*V.T();
            EVECT = EVECT*V.T();
            k = 0; eps = 0.0;
            for(i=0;i<temp.num_of_rows;i++){
                for(j=0;j<temp.num_of_cols;j++){
                    if(i!=j) {eps+=temp.M[k]*temp.M[k];}
                    k++;
                }// for j
            }// for i  

       }// while eps>EPS

       EVAL = temp;


}
double MATRIX::Determinant(){

  double res =0.0;
  if((num_of_cols==3)&&(num_of_rows==3)){

    res = M[0]*(M[4]*M[8]-M[5]*M[7])-
          M[3]*(M[1]*M[8]-M[2]*M[7])+
          M[6]*(M[1]*M[5]-M[2]*M[4]);
  }
  else{
    cout<<"Error in MATRIX::Determinant. The function is only implemented for 3x3 matrix\n";
    exit(0);
  }

  return res;
}

double MATRIX::dot_product(MATRIX& B){
// A->dot_product(B) = sum_of_elements_of (A^T * B)

  double res = 0.0;

  if((num_of_cols==num_of_rows) && (num_of_cols==B.num_of_cols) && (num_of_rows=B.num_of_rows)){  

    MATRIX* X; X = new MATRIX(num_of_rows,num_of_cols);

    *X = (*this).T() * B;

    for(int i=0;i<num_of_elems;i++){  res += X->M[i];   }

    delete X;


  }else{
    cout<<"Error in MATRIX::dot_product. Matrices must be square and of the same rank.\n";
    exit(0);
  }

  return res;
  
}

double MATRIX::dot_product(MATRIX* B){

  double res = 0.0;

  if((num_of_cols==num_of_rows) && (num_of_cols==B->num_of_cols) && (num_of_rows=B->num_of_rows)){

    MATRIX* X; X = new MATRIX(num_of_rows,num_of_cols);

    *X = (*this).T() * (*B);

    for(int i=0;i<num_of_elems;i++){  res += X->M[i];   }

    delete X;


  }else{
    cout<<"Error in MATRIX::dot_product. Matrices must be square and of the same rank.\n";
    exit(0);
  }

  return res;
  
}




double MATRIX::tr(){
// Trace
  double res =0.0;
  if(num_of_cols==num_of_rows){
    for(int i=0;i<num_of_cols;i++){
      res += M[i*num_of_cols+i];
    }
  }

  return res;

}
void MATRIX::Rotation(const VECTOR& u){
/*  This function initializes the matrix to a rotation matrix
    for rotation around the axis given by direction of the vector: (u/|u|)
    on amount given by norm of vector u: |u|.
    If the vector has a zero length - this will be an identity matrix
*/
    double umod,umod2;

    umod2 = u.x*u.x + u.y*u.y + u.z*u.z;

    if(umod2==0.0){
        M[0] = 1.0;  M[1] = 0.0;  M[2] = 0.0;
        M[3] = 0.0;  M[4] = 1.0;  M[5] = 0.0;
        M[6] = 0.0;  M[7] = 0.0;  M[8] = 1.0;
    }
    else{
        umod = sqrt(umod2);
        VECTOR n(u.x/umod, u.y/umod, u.z/umod);

        double cs,sn;
        cs = (1.0 - cos(umod));
        sn = sin(umod);

        /* 
        Here is efficient implementation of Rodrigues formula: 
        M = I + W * sin_psi + W*W*(1.0 - cos_psi)
        where I - identity matrix
                            | 0  -n.z  n.y |
              W = skew(n) = | n.z  0  -n.x |
                            |-n.y n.x  0   |

        */

        double x,y,z,xy,xz,yz,x2,y2,z2;
        x  = n.x * sn;
        y  = n.y * sn;
        z  = n.z * sn;
        x2 = (n.x * n.x - 1.0) * cs;
        y2 = (n.y * n.y - 1.0) * cs;
        z2 = (n.z * n.z - 1.0) * cs;
        xy = n.x * n.y * cs;
        xz = n.x * n.z * cs;
        yz = n.y * n.z * cs;

        M[0] = 1.0 + x2;  M[1] =-z   + xy;  M[2] = y   + xz;
        M[3] = z   + xy;  M[4] = 1.0 + y2;  M[5] =-x   + yz;
        M[6] =-y   + xz;  M[7] = x   + yz;  M[8] = 1.0 + z2;


    }

}
void MATRIX::Rx(double phi){
  double cs = cos(phi);
  double sn = sin(phi);
  M[0] = 1.0; M[1] = 0.0; M[2] = 0.0;
  M[3] = 0.0; M[4] = cs;  M[5] = -sn;
  M[6] = 0.0; M[7] = sn;  M[8] = cs;

}
void MATRIX::Ry(double phi){
  double cs = cos(phi);
  double sn = sin(phi);
  M[0] = cs;  M[1] = 0.0; M[2] = sn;
  M[3] = 0.0; M[4] = 1.0; M[5] = 0.0;
  M[6] = -sn; M[7] = 0.0; M[8] = cs;

}
void MATRIX::Rz(double phi){
  double cs = cos(phi);
  double sn = sin(phi);
  M[0] = cs;  M[1] = -sn; M[2] = 0.0;
  M[3] = sn;  M[4] = cs;  M[5] = 0.0;
  M[6] = 0.0; M[7] = 0.0; M[8] = 1.0;

}

void MATRIX::tensor_product(VECTOR v1,VECTOR v2){
  M[0] = v1.x*v2.x;   M[1] = v1.x*v2.y;  M[2] = v1.x*v2.z;
  M[3] = v1.y*v2.x;   M[4] = v1.y*v2.y;  M[5] = v1.y*v2.z;
  M[6] = v1.z*v2.x;   M[7] = v1.z*v2.y;  M[8] = v1.z*v2.z;
}

double MATRIX::max_elt(){
  double x = fabs(M[0]);
  double y;
  for(int i=0;i<num_of_elems;i++){  y = fabs(M[i]); if(y>=x){ x = y; } }
  return x;
}


void MATRIX::bin_dump(std::string filename){

  std::ofstream f(filename.c_str(), ios::out|ios::binary);

  if(f.is_open()){
    f.seekp(0);
    f.write((char*)M, sizeof(double)*num_of_elems);
    f.close();    
  }
  else{  cout<<"File "<<filename<<" cann't be open\n"; }
}
 
void MATRIX::bin_load(std::string filename){

  std::ifstream f(filename.c_str(), ios::in|ios::binary);

  if(f.is_open()){
    f.seekg(0);
    f.read((char *)M, sizeof(double)*num_of_elems);
    f.close();   
  }
  else{  cout<<"File "<<filename<<" cann't be open\n"; }


}



void pop_submatrix(MATRIX* X,MATRIX* x,vector<int>& subset){
// Extract the submatrix x from the matrix X according to indices given in <subset>
// Assume that memory for x is already allocated and its dimensions are consistent with the map dimensions:
// subset_size == x->num_of_cols = x->num_of_rows = subset.size()
// X->num_of_cols = X->num_of_rows >= subset_size

  int N = X->num_of_cols;
  int n = x->num_of_cols;

  int i,j,a,b;

  for(i=0;i<n;i++){
    a = subset[i];
    for(j=0;j<n;j++){      
      b = subset[j];
      x->M[i*n+j] = X->M[a*N+b];
    }// j
  }// i

}// void pop_submatrix(MATRIX* X,MATRIX* x,vector<int>& subset)

void push_submatrix(MATRIX* X,MATRIX* x,vector<int>& subset){
// Pushes the smaller submatrix x back to the bigger matrix X, according to indices given in <subset>
// Assume that memory for x is already allocated and its dimensions are consistent with the map dimensions:
// subset_size == x->num_of_cols = x->num_of_rows = subset.size()
// X->num_of_cols = X->num_of_rows >= subset_size

  int N = X->num_of_cols;
  int n = x->num_of_cols;

  int i,j,a,b;

  for(i=0;i<n;i++){
    a = subset[i];
    for(j=0;j<n;j++){      
      b = subset[j];
      X->M[a*N+b] = x->M[i*n+j];
    }// j
  }// i

}// void push_submatrix(MATRIX* X,MATRIX* x,vector<int>& subset)




void set_value(int& is_defined, MATRIX& value,boost::python::object obj, std::string attrName){

  int has_attr=0;
  has_attr = (int)hasattr(obj,attrName);
  if(has_attr){
      value = extract<MATRIX>(obj.attr(attrName.c_str()));
      is_defined = 1;
  }
}


// ----------- Save --------------
void save(boost::property_tree::ptree& pt,std::string path,MATRIX& vt){
  for(int i=0;i<vt.num_of_elems;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    pt.put(path+"."+rt,vt.M[i]);
  }
}

void save(boost::property_tree::ptree& pt,std::string path,vector<MATRIX>& vt){
  int sz = vt.size();
  for(int i=0;i<sz;i++){
    stringstream ss(stringstream::in | stringstream::out);
    std::string rt; ss<<i; ss>>rt;
    save(pt,path+"."+rt,vt[i]);
  }
}

// ----------- Load --------------
void load(boost::property_tree::ptree& pt,std::string path, MATRIX& vt, int& status){ std::cout<<"Sorry: load function for MATRIX object is not defined yet\n"; exit(0); }
void load(boost::property_tree::ptree& pt,std::string path,vector<MATRIX>& vt,int& status){ std::cout<<"Sorry: load function for vector<MATRIX> object is not defined yet\n"; exit(0); }




}// namespace liblinalg
}// namespace libmmath


