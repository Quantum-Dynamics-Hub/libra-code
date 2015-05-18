#include "Basis.h"
#include "Nuclear.h"

/****************************************************************************

  This file contains following functions:

  void update_overlap_matrix(int x_period,int y_period,int z_period,const VECTOR& t1, const VECTOR& t2, const VECTOR& t3,
                             vector<int>& basis_fo,vector<AO>& basis_ao,
                             MATRIX* Sao,vector<double*>& aux, int n_aux, Nuclear& mol)


****************************************************************************/

void update_overlap_matrix(int x_period,int y_period,int z_period,const VECTOR& t1, const VECTOR& t2, const VECTOR& t3,
                           vector<int>& basis_fo,vector<AO>& basis_ao,
                           MATRIX* Sao,vector<double*>& aux, int n_aux, Nuclear& mol){

  int i,j,n,I,J;
  VECTOR dIdA,dIdB,TV, Rij,Rij0;
  double dist, dist_min;  
  int opt_x,opt_y,opt_z;

  int Norb = basis_fo.size();

  *Sao = 0.0;
 

  for(i=0;i<Norb;i++){
    I = basis_fo[i];
    for(j=i;j<Norb;j++){      
      J = basis_fo[j];

      for(int nx=-x_period;nx<=x_period;nx++){
        for(int ny=-y_period;ny<=y_period;ny++){
          for(int nz=-z_period;nz<=z_period;nz++){

            // This summation corresponds to k = 0 (Gamma-point)
            TV = nx*t1 + ny*t2 + nz*t3;
            Sao->M[i*Norb+j] += OVERLAP_INTEGRAL(basis_ao[I],basis_ao[J],0,dIdA,dIdB,aux,n_aux,TV);

          }// for nz
        }// for ny
      }// for nx 

      Sao->M[j*Norb+i] = Sao->M[i*Norb+j];

    }// j
  }// i



/*   !!!!!!!!!!!!!! Don't need the below part untill we get to the gradients !!!!!!!!!!!

  else{

    for(n=0;n<Nnucl;n++){
      *dSao_dx[n] = 0.0;
      *dSao_dy[n] = 0.0;
      *dSao_dz[n] = 0.0;
    }

  
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){

        Sao->M[i*Norb+j] = OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);

        dSao_dx[basis_ao[i].at_indx]->M[i*Norb+j] += dIdA.x;
        dSao_dx[basis_ao[j].at_indx]->M[i*Norb+j] += dIdB.x;

        dSao_dy[basis_ao[i].at_indx]->M[i*Norb+j] += dIdA.y;
        dSao_dy[basis_ao[j].at_indx]->M[i*Norb+j] += dIdB.y;

        dSao_dz[basis_ao[i].at_indx]->M[i*Norb+j] += dIdA.z;
        dSao_dz[basis_ao[j].at_indx]->M[i*Norb+j] += dIdB.z;
*/
/*
      // Tests
      double dx = 0.0001;
      R[basis_ao[i].at_indx] += VECTOR(dx,0.0,0.0);
      double s2 = OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);
      R[basis_ao[i].at_indx] -= VECTOR(2.0*dx,0.0,0.0);
      double s1 = OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);
      R[basis_ao[i].at_indx] += VECTOR(dx,0.0,0.0);
      OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);

      cout<<"Overlap:  der(numer)= "<<0.5*(s2-s1)/dx<<" der(anal)= "<<dIdA.x<<endl;
*/
     

//    }// for j
//  }// for i
//  }// else : method!=eht0

//  for(int n=0;n<Nnucl;n++){
//    cout<<"dSao_dx["<<n<<"]= \n"<<*dSao_dx[n]<<endl;
//  }

//  cout<<"Overlap matrix in AO basis:\n";
//  cout<<*Sao<<endl;
  // So far use my inverter:
//  Sao->Inverse(Sao_inv);

}


void update_overlap_matrix_new(int x_period,int y_period,int z_period,const VECTOR& t1, const VECTOR& t2, const VECTOR& t3,
                               vector<int>& basis_fo,vector<AO>& basis_ao,
                               MATRIX* Sao,vector<double*>& aux, int n_aux, Nuclear& mol){

  int i,j,n,I,J;
  VECTOR dIdA,dIdB,TV, Rij,Rij0;
  double dist, dist_min;  
  int opt_x,opt_y,opt_z;

  int Norb = basis_fo.size();

  *Sao = 0.0;
 

  for(i=0;i<Norb;i++){
    I = basis_fo[i];
    for(j=i;j<Norb;j++){      
      J = basis_fo[j];

      for(int nx=-x_period;nx<=x_period;nx++){
        for(int ny=-y_period;ny<=y_period;ny++){
          for(int nz=-z_period;nz<=z_period;nz++){

            // This summation corresponds to k = 0 (Gamma-point)
            TV = nx*t1 + ny*t2 + nz*t3;

//            Sao->M[i*Norb+j] += OVERLAP_INTEGRAL_new(basis_ao[I],basis_ao[J],0,dIdA,dIdB,aux,n_aux,TV);

// Just test pseudopotential - this is too damn expensive:

            for(int a=0;a<mol.Nnucl;a++){

              Sao->M[i*Norb+j] += PSEUDO_02_INTEGRAL(0.0, 1.0,  1.0, mol.R[a], 
                                                     basis_ao[I],basis_ao[J],0,dIdA,dIdB,aux,n_aux,TV);
            }


//double PSEUDO_02_INTEGRAL(double C0, double C2, double alp, VECTOR& R,
//                          AO& aoa,AO& aob,int is_normalize,VECTOR& dIdA,VECTOR& dIdB,vector<double*>& aux,int n_aux,
//                          VECTOR& TV){



          }// for nz
        }// for ny
      }// for nx 

      Sao->M[j*Norb+i] = Sao->M[i*Norb+j];

    }// j
  }// i



/*   !!!!!!!!!!!!!! Don't need the below part untill we get to the gradients !!!!!!!!!!!

  else{

    for(n=0;n<Nnucl;n++){
      *dSao_dx[n] = 0.0;
      *dSao_dy[n] = 0.0;
      *dSao_dz[n] = 0.0;
    }

  
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){

        Sao->M[i*Norb+j] = OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);

        dSao_dx[basis_ao[i].at_indx]->M[i*Norb+j] += dIdA.x;
        dSao_dx[basis_ao[j].at_indx]->M[i*Norb+j] += dIdB.x;

        dSao_dy[basis_ao[i].at_indx]->M[i*Norb+j] += dIdA.y;
        dSao_dy[basis_ao[j].at_indx]->M[i*Norb+j] += dIdB.y;

        dSao_dz[basis_ao[i].at_indx]->M[i*Norb+j] += dIdA.z;
        dSao_dz[basis_ao[j].at_indx]->M[i*Norb+j] += dIdB.z;
*/
/*
      // Tests
      double dx = 0.0001;
      R[basis_ao[i].at_indx] += VECTOR(dx,0.0,0.0);
      double s2 = OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);
      R[basis_ao[i].at_indx] -= VECTOR(2.0*dx,0.0,0.0);
      double s1 = OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);
      R[basis_ao[i].at_indx] += VECTOR(dx,0.0,0.0);
      OVERLAP_INTEGRAL(basis_ao[i],basis_ao[j],0,dIdA,dIdB);

      cout<<"Overlap:  der(numer)= "<<0.5*(s2-s1)/dx<<" der(anal)= "<<dIdA.x<<endl;
*/
     

//    }// for j
//  }// for i
//  }// else : method!=eht0

//  for(int n=0;n<Nnucl;n++){
//    cout<<"dSao_dx["<<n<<"]= \n"<<*dSao_dx[n]<<endl;
//  }

//  cout<<"Overlap matrix in AO basis:\n";
//  cout<<*Sao<<endl;
  // So far use my inverter:
//  Sao->Inverse(Sao_inv);

}

