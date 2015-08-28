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

#include "DOS.h"
//#include "units.h"

#include <sstream>
using namespace std;


/*********************************************************************************
   This file contains the following functions:
   

*********************************************************************************/

void compute_dos(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                 vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                 Electronic* el, Memory* mem){

// Calculation of the total density of states based on the density matrix

  int iter = 0;
  int a,b,c,d,ii,jj,i,j,k,kk_a,kk_b,n,I,J,A;
  double Eelec_prev,dE,den_err;
  std::string eigen_method; 
  std::string core_method;

  double all_tot_s = 0.0;
  double all_tot_p = 0.0;
  double all_tot_d = 0.0;


  // Complexity: Norb x N_AO x Norb         


  for(int a=0;a<fragment.size();a++){  // loop over all atoms in the given fragment
 
    A = fragment[a];                   // global atom index for a-th atom in the fragment

    stringstream ss(stringstream::in | stringstream::out);
    std::string out;
    (ss << A);  ss >> out;

    
    // Do calculations for given atom A
    FILE* fpa; 
    std::string filename;  filename = prms.dos_prefix+"_alpha_wfc_atom" + out;    
    fpa = fopen(filename.c_str(),"w");

    FILE* fpb; 
    filename = prms.dos_prefix+"_beta_wfc_atom" + out;    
    fpb = fopen(filename.c_str(),"w");



    fprintf(fpa,"   Energy(alpha)   Total(alpha)    DOS(s,alpha)   DOS(p,alpha)   DOS(d,alpha)\n");
    fprintf(fpb,"   Energy(beta)    Total(beta)     DOS(s,beta)    DOS(p,beta)    DOS(d,beta)\n");


    // total s, p, d populations on this atom (integral over all energies)
    // alpha
    double at_tot_s_a = 0.0;
    double at_tot_p_a = 0.0;
    double at_tot_d_a = 0.0;
    // beta
    double at_tot_s_b = 0.0;
    double at_tot_p_b = 0.0;
    double at_tot_d_b = 0.0;




    for(k=0;k<el->Norb;k++){      // loop over all energies (ordered bands)

      kk_a = el->bands_alp[k].first;  // local index of the MO with this energy (usually it is just "k")
      kk_b = el->bands_bet[k].first;  // local index of the MO with this energy (usually it is just "k")


      // energy-specific s, p, d populations on this atom 
      //(projectoin of MO with given energy on AOs of corresponding type)
      // alpha
      double pops_a = 0.0;
      double popp_a = 0.0;
      double popd_a = 0.0;

      // beta
      double pops_b = 0.0;
      double popp_b = 0.0;
      double popd_b = 0.0;

      
      for(n=0;n<at_orbitals[A].size();n++){      // Loop over all AOs on atom A

          int I = at_orbitals[A][n];             // global index of n-th AO centered on atom A
                                             

          // Well, this is quick: find the local index (i) of the global I AO withing AO basis of given fragment          
          // so i is such that: basis_fo[i] = I - is the global index of the i-th AO in the given fragment 
          i = -1;
          for(jj=0;jj<basis_fo.size();jj++){    if(basis_fo[jj]==I){  i = jj; }             }



          // alpha
          double tmp = 0.0;
          ii = el->bands_alp[i].first;         

          for(j=0;j<el->Norb;j++){       
            jj = el->bands_alp[j].first;

            tmp += el->C_alp->M[ii*el->Norb+kk_a] * el->C_alp->M[jj*el->Norb+kk_a] * el->Sao->M[jj*el->Norb+ii]; 

          }// for j
          if(basis_ao[I].ao_shell_type=="s"){ pops_a += tmp;  }
          else if(basis_ao[I].ao_shell_type=="p"){ popp_a += tmp;  }
          else if(basis_ao[I].ao_shell_type=="d"){ popd_a += tmp;  }



          // beta
          tmp = 0.0;
          ii = el->bands_bet[i].first;         

          for(j=0;j<el->Norb;j++){       
            jj = el->bands_bet[j].first;

            tmp += el->C_bet->M[ii*el->Norb+kk_b] * el->C_bet->M[jj*el->Norb+kk_b] * el->Sao->M[jj*el->Norb+ii]; 

          }// for j
          if(basis_ao[I].ao_shell_type=="s"){ pops_b += tmp;  }
          else if(basis_ao[I].ao_shell_type=="p"){ popp_b += tmp;  }
          else if(basis_ao[I].ao_shell_type=="d"){ popd_b += tmp;  }

          
      }// for n - over all AOs on atom A

      // Sum rule for number of electrons
      /*
      at_tot_s_a += pops_a*el->occ_alp[k].second;
      at_tot_p_a += popp_a*el->occ_alp[k].second;
      at_tot_d_a += popd_a*el->occ_alp[k].second;

      at_tot_s_b += pops_b*el->occ_bet[k].second;
      at_tot_p_b += popp_b*el->occ_bet[k].second;
      at_tot_d_b += popd_b*el->occ_bet[k].second;
      */

      // Sum rule for number of states
      at_tot_s_a += pops_a;
      at_tot_p_a += popp_a;
      at_tot_d_a += popd_a;

      at_tot_s_b += pops_b;
      at_tot_p_b += popp_b;
      at_tot_d_b += popd_b;


     
      fprintf(fpa," %8.3f  %8.3f %8.3f %8.3f %8.3f   \n",el->bands_alp[k].second, pops_a+popp_a+popd_a, pops_a, popp_a, popd_a );
      fprintf(fpb," %8.3f  %8.3f %8.3f %8.3f %8.3f   \n",el->bands_bet[k].second, pops_b+popp_b+popd_b, pops_b, popp_b, popd_b );

    }// for k - all MOs (energy levels)


    fprintf(fpa,"atomic total(s,alpha) %8.3f  \n",at_tot_s_a );
    fprintf(fpa,"atomic total(p,alpha) %8.3f  \n",at_tot_p_a );
    fprintf(fpa,"atomic total(d,alpha) %8.3f  \n",at_tot_d_a );
    fprintf(fpa,"atomic total,alpha    %8.3f  \n",at_tot_s_a + at_tot_p_a + at_tot_d_a );

    all_tot_s += at_tot_s_a;
    all_tot_p += at_tot_p_a;
    all_tot_d += at_tot_d_a;

    fprintf(fpb,"atomic total(s,beta) %8.3f  \n",at_tot_s_b );
    fprintf(fpb,"atomic total(p,beta) %8.3f  \n",at_tot_p_b );
    fprintf(fpb,"atomic total(d,beta) %8.3f  \n",at_tot_d_b );
    fprintf(fpb,"atomic total,beta    %8.3f  \n",at_tot_s_b + at_tot_p_b + at_tot_d_b );

    all_tot_s += at_tot_s_b;
    all_tot_p += at_tot_p_b;
    all_tot_d += at_tot_d_b;


    fclose(fpa);
    fclose(fpb);

  }// for a - all nuclei


  cout<<"all_tot_s = "<<all_tot_s<<endl;
  cout<<"all_tot_p = "<<all_tot_p<<endl;
  cout<<"all_tot_d = "<<all_tot_d<<endl;
  cout<<"all_tot = "<<all_tot_s+all_tot_p+all_tot_d<<endl;

}

