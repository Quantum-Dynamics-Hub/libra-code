/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file DOS.cpp
  \brief This file implements function for computing Densities of States (atomically and orbitally-resolved)
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <sstream>
#endif 

#include "DOS.h"

using namespace std;

/// liblibra namespace
namespace liblibra{


/// libqchem_tools namespace
namespace libqchem_tools{



void compute_dos
( Electronic_Structure& el, vector<AO>& basis_ao, Control_Parameters& prms,
  vector<int>& fragment, vector< vector<int> >& atom_to_ao_map
){
/** Calculation of the total density of states based on the density matrix
  
  \param[in] el Variable containig the electronic variables
  \param[in] basis_ao AO basis for the entire system (the one used to produce electronic structure)
  \param[in] prms The parameters controlling DOS calculations
  \param[in] fragment The list of integers - the atom indices (not the IDs!) onto which we project wavefunction
  \param[in] atom_to_ao_map The mapping between global atomic and AO indices, so atom_to_ao_map[a][i] is the global index
  of i-th orbital of the atom with index a
*/

  double all_tot_s = 0.0;
  double all_tot_p = 0.0;
  double all_tot_d = 0.0;


  // Complexity: Norb x N_AO x Norb         
  for(int a=0;a<fragment.size();a++){  // loop over all atoms in the given fragment
 
    int A = fragment[a];                   // global atom index for a-th atom in the fragment

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




    for(int k=0;k<el.Norb;k++){      // loop over all energies (ordered bands)

      int kk_a = el.bands_alp[k].first;  // local index of the MO with this energy (usually it is just "k")
      int kk_b = el.bands_bet[k].first;  // local index of the MO with this energy (usually it is just "k")


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

      
      for(int n=0;n<atom_to_ao_map[A].size();n++){  // Loop over all AOs on atom A
          int I = atom_to_ao_map[A][n];         // global index of n-th AO centered on atom A
                                             
          // alpha
          double tmp = 0.0;     
          int ii = el.bands_alp[I].first;         

          for(int j=0;j<el.Norb;j++){       
            int jj = el.bands_alp[j].first;
            tmp += el.C_alp->M[ii*el.Norb+kk_a] * el.C_alp->M[jj*el.Norb+kk_a] * el.Sao->M[jj*el.Norb+ii]; 

          }// for j

          if(basis_ao[I].ao_shell_type=="s"){ pops_a += tmp;  }
          else if(basis_ao[I].ao_shell_type=="p"){ popp_a += tmp;  }
          else if(basis_ao[I].ao_shell_type=="d"){ popd_a += tmp;  }



          // beta
          tmp = 0.0;
          ii = el.bands_bet[I].first;         

          for(int j=0;j<el.Norb;j++){       
            int jj = el.bands_bet[j].first;
            tmp += el.C_bet->M[ii*el.Norb+kk_b] * el.C_bet->M[jj*el.Norb+kk_b] * el.Sao->M[jj*el.Norb+ii]; 

          }// for j

          if(basis_ao[I].ao_shell_type=="s"){ pops_b += tmp;  }
          else if(basis_ao[I].ao_shell_type=="p"){ popp_b += tmp;  }
          else if(basis_ao[I].ao_shell_type=="d"){ popd_b += tmp;  }

          
      }// for n - over all AOs on atom A

      // Sum rule for number of states
      at_tot_s_a += pops_a;
      at_tot_p_a += popp_a;
      at_tot_d_a += popd_a;

      at_tot_s_b += pops_b;
      at_tot_p_b += popp_b;
      at_tot_d_b += popd_b;


     
      fprintf(fpa," %8.3f  %8.3f %8.3f %8.3f %8.3f   \n",el.bands_alp[k].second, pops_a+popp_a+popd_a, pops_a, popp_a, popd_a );
      fprintf(fpb," %8.3f  %8.3f %8.3f %8.3f %8.3f   \n",el.bands_bet[k].second, pops_b+popp_b+popd_b, pops_b, popp_b, popd_b );

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


void compute_dos
( Electronic_Structure& el, vector<AO>& basis_ao, Control_Parameters& prms,
  boost::python::list fragment, vector< vector<int> >& atom_to_ao_map
){
/** Calculation of the total density of states based on the density matrix
  Python-friendly version.
  
  \param[in] el Variable containig the electronic variables
  \param[in] basis_ao AO basis for the entire system (the one used to produce electronic structure)
  \param[in] prms The parameters controlling DOS calculations
  \param[in] fragment The list of integers - the atom indices (not the IDs!) onto which we project wavefunction
  \param[in] atom_to_ao_map The mapping between global atomic and AO indices, so atom_to_ao_map[a][i] is the global index
  of i-th orbital of the atom with index a
*/


  int sz = len(fragment);
  vector<int> fragm(sz,0);
  for(int i=0;i<sz;i++){ fragm[i] = boost::python::extract<double>(fragment[i]);   }

  compute_dos(el, basis_ao, prms, fragm, atom_to_ao_map);

}


void compute_dos
( Electronic_Structure& el, System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, vector< vector<int> >& atom_to_ao_map
){
/** Calculation of the total density of states based on the density matrix
  In this case, we consider all atoms of the system
  
  \param[in] el Variable containig the electronic variables
  \param[in] syst Variable containig nuclear variables - the chemical structure
  \param[in] basis_ao AO basis for the entire system (the one used to produce electronic structure)
  \param[in] prms The parameters controlling DOS calculations
  \param[in] atom_to_ao_map The mapping between global atomic and AO indices, so atom_to_ao_map[a][i] is the global index
  of i-th orbital of the atom with index a
*/


  int sz = syst.Number_of_atoms;
  vector<int> fragm(sz,0);
  for(int i=0;i<sz;i++){ fragm[i] = i;  }

  compute_dos(el, basis_ao, prms, fragm, atom_to_ao_map);
}




}// namespace libqchem_tools
}// liblibra
