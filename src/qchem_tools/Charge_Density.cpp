/*********************************************************************************
* Copyright (C) 2016-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Charge_Density.cpp
  \brief This file implements functions for computing and printing charge density information

*/


#if defined(USING_PCH)
#include "../pch.h"
#else
#include <sstream>
#endif

#include "Charge_Density.h"
#include "../converters/libconverters.h"

/// liblibra namespace
namespace liblibra{

using namespace std;
using namespace libconverters;


/// libqchem_tools namespace
namespace libqchem_tools{



void charge_density( Electronic_Structure& el, System& syst, vector<AO>& basis_ao, Control_Parameters& prms){
/**
  \param[in] el The electronic structure of the system
  \param[in] syst The nuclear structure of the system
  \param[in] basis_ao The atomic basis used to compute electronic structure
  \param[in] prms Parameters controlling how to execute the calculations

  This function compute wfc of a given state (orbital) and
  prints it in CUBE format for visualization
  See about Gaussian CUBE format at: http://paulbourke.net/dataformats/cube/

  prms.charge_density_prefix - specifies the directory to which the files will be written
  prms.orbs - the indices of the MOs to print
*/

  cout<<"Printing parameters...\n";
  cout<<"Orbitals to handle: "; for(int i=0;i<prms.orbs.size();i++){ cout<<prms.orbs[i]<<" "; } cout<<endl;
  cout<<"prms.nx_grid = "<<prms.nx_grid<<endl;
  cout<<"prms.ny_grid = "<<prms.ny_grid<<endl;
  cout<<"prms.nz_grid = "<<prms.nz_grid<<endl;
  cout<<"prms.charge_density_prefix = "<<prms.charge_density_prefix<<endl;

  for(int i=0;i<prms.orbs.size();i++){

    int orb = prms.orbs[i];

    cout<<"Computing charge density for orbital "<<orb<<endl;

    stringstream ss(stringstream::in | stringstream::out);
    std::string out;
    (ss << orb);  ss >> out;
  
    FILE* fp;
    std::string filename;
    filename = prms.charge_density_prefix+"_orbital_" + out+".cube";    
    fp = fopen(filename.c_str(),"w");
    if(fp==NULL){
      cout<<"Error: Can not create/open file with prefix "<<prms.charge_density_prefix<<endl;
      cout<<"If this prefix is the directory name, please create the directory first\n";
      exit(0);
    }
  
  
    fprintf(fp,"EHT CHARGE DENSITY  \n");
    fprintf(fp,"Comment line  \n");
  
  
    // Determine box around
    VECTOR min_pos; min_pos = 0.0;
    VECTOR max_pos; max_pos = 0.0;
  
    for(int n=0;n<syst.Number_of_atoms;n++){
      double X = syst.Atoms[n].Atom_RB.rb_cm.x;
      double Y = syst.Atoms[n].Atom_RB.rb_cm.y;
      double Z = syst.Atoms[n].Atom_RB.rb_cm.z;

      if(X < min_pos.x) { min_pos.x = X; }
      if(Y < min_pos.y) { min_pos.y = Y; }
      if(Z < min_pos.z) { min_pos.z = Z; }
  
      if(X > max_pos.x) { max_pos.x = X; }
      if(Y > max_pos.y) { max_pos.y = Y; }
      if(Z > max_pos.z) { max_pos.z = Z; }
  
    }
    // Add padding
    min_pos -= 5.0;
    max_pos += 5.0;

    //cout<<"Cube dimension = "<<min_pos<<"  "<<max_pos<<endl;
  
  
    fprintf(fp,"%5i%12.6f%12.6f%12.6f\n",syst.Number_of_atoms, min_pos.x, min_pos.y, min_pos.z);
  
    // Size of voxels
    VECTOR dr;
    dr.x = (max_pos.x - min_pos.x)/float(prms.nx_grid); 
    dr.y = (max_pos.y - min_pos.y)/float(prms.ny_grid); 
    dr.z = (max_pos.z - min_pos.z)/float(prms.nz_grid); 
    fprintf(fp,"%5i%12.6f%12.6f%12.6f\n",prms.nx_grid, dr.x,0.00,0.00);
    fprintf(fp,"%5i%12.6f%12.6f%12.6f\n",prms.ny_grid, 0.00,dr.y,0.00);
    fprintf(fp,"%5i%12.6f%12.6f%12.6f\n",prms.nz_grid, 0.00,0.00,dr.z);
  
  
  
    for(int n=0;n<syst.Number_of_atoms;n++){
      fprintf(fp,"%5i%12.6f%12.6f%12.6f%12.6f\n",syst.Atoms[n].Atom_Z, 0.0, syst.Atoms[n].Atom_RB.rb_cm.x, syst.Atoms[n].Atom_RB.rb_cm.y, syst.Atoms[n].Atom_RB.rb_cm.z); 
    }
  
  
    for(int ix=0;ix<prms.nx_grid;ix++){
      for(int iy=0;iy<prms.ny_grid;iy++){
        for(int iz=0;iz<prms.nz_grid;iz++){
  
           VECTOR pos;
           pos.x = min_pos.x + ix * dr.x;
           pos.y = min_pos.y + iy * dr.y;
           pos.z = min_pos.z + iz * dr.z;
           
           
           double psi = 0.0;
           for(int a=0;a<el.Norb;a++){
             psi += el.C_alp->M[a*el.Norb+orb] * basis_ao[a].compute( pos );
           }

           //cout<<ix<<"  "<<iy<<"  "<<iz<<"  "<<psi<<endl;
  
           double charg_dens = psi; // * psi;
           
           
           fprintf(fp, "%g ",charg_dens);  // charg_dens = data[ix][iy][iz]
  
          if (iz % 6 == 5)   fprintf(fp,"\n");
        }// for iz
        fprintf(fp,"\n");
      }// for iy
    }// for ix

    fclose(fp);

  }// for i

}// charge_density




void charge_density(MATRIX& C, vector<listHamiltonian_QM>& ham, System& syst, vector<vector<int> >& active_orb, Control_Parameters& prms){
/**
  \brief To compute the charge densities for the combinations of fragment MOs

  \param[in] C The matrix containing the transformation coefficients = in which combintation to take the MOs of the sub-systems
  \param[in] ham Is the list of Hamiltonian objects (of type listHamiltonian), each containing both electronic structure (MO)
             information and the basis (AO) information about each sub-system  
  \param[in] syst The nuclear structure of the system (combined)
  \param[in] active_orb Holds indices of the fragment MOs each fragment contributes to the total pool of the fragment MOs             
  \param[in] prms Parameters controlling how to execute the calculations

  This function compute wfc of a given state (orbital) and
  prints it in CUBE format for visualization
  See about Gaussian CUBE format at: http://paulbourke.net/dataformats/cube/

  prms.charge_density_prefix - specifies the directory to which the files will be written
  prms.orbs - the indices of the MOs to print
*/

//vector<AO>& basis_ao, 
//Electronic_Structure& el
  int i;

  cout<<"Printing parameters...\n";
  cout<<"Orbitals to handle: "; for(int i=0;i<prms.orbs.size();i++){ cout<<prms.orbs[i]<<" "; } cout<<endl;
  cout<<"prms.nx_grid = "<<prms.nx_grid<<endl;
  cout<<"prms.ny_grid = "<<prms.ny_grid<<endl;
  cout<<"prms.nz_grid = "<<prms.nz_grid<<endl;
  cout<<"prms.charge_density_prefix = "<<prms.charge_density_prefix<<endl;


  int nfrags = active_orb.size();  // this is the total number of fragments
  int nfmo = C.n_rows;        // total number of FMOs

  int summ = 0;
  for(i=0;i<nfrags;i++){
    summ += active_orb[i].size();  // how many FMOs does the fragment i contribute to the total pool of the orbitals
  }
  if(summ!=nfmo){
    cout<<"Error in void charge_density(...) : The total number of FMOs due to all fragments is not equal to the number expected from the LC matrix\n";
    cout<<"the # of fragments = "<<nfrags<<endl;
    cout<<"the total # of FMOs from all fragments = "<<summ<<endl;
    cout<<"the # of expected FMOs = "<<nfmo<<endl;
    cout<<"Exiting now...\n"; exit(0);
  }


  for(i=0;i<prms.orbs.size();i++){ // Note that these "orbs" will now have a meaning of the superpositions of fragment states

    int orb = prms.orbs[i];

    if(orb>=nfmo){
      cout<<"Error in void charge_density(...) : can only compute orbitals with indices below the total # of FMOs\n";
      cout<<"Note: the indices denote the FMO states\n";
      cout<<"You tried to compute orbital = "<<orb<<endl;
      cout<<"...but the # of FMOs is = "<<nfmo<<endl;
      cout<<"Exiting now...\n"; exit(0);
    }

    cout<<"Computing charge density for orbital "<<orb<<endl;

    stringstream ss(stringstream::in | stringstream::out);
    std::string out;
    (ss << orb);  ss >> out;
  
    FILE* fp;
    std::string filename;
    filename = prms.charge_density_prefix+"_orbital_" + out+".cube";    
    fp = fopen(filename.c_str(),"w");
    if(fp==NULL){
      cout<<"Error: Can not create/open file with prefix "<<prms.charge_density_prefix<<endl;
      cout<<"If this prefix is the directory name, please create the directory first\n";
      exit(0);
    }
  
  
    fprintf(fp,"EHT CHARGE DENSITY  \n");
    fprintf(fp,"Comment line  \n");
  
  
    // Determine box around
    VECTOR min_pos; min_pos = 0.0;
    VECTOR max_pos; max_pos = 0.0;
  
    for(int n=0;n<syst.Number_of_atoms;n++){
      double X = syst.Atoms[n].Atom_RB.rb_cm.x;
      double Y = syst.Atoms[n].Atom_RB.rb_cm.y;
      double Z = syst.Atoms[n].Atom_RB.rb_cm.z;

      if(X < min_pos.x) { min_pos.x = X; }
      if(Y < min_pos.y) { min_pos.y = Y; }
      if(Z < min_pos.z) { min_pos.z = Z; }
  
      if(X > max_pos.x) { max_pos.x = X; }
      if(Y > max_pos.y) { max_pos.y = Y; }
      if(Z > max_pos.z) { max_pos.z = Z; }
  
    }
    // Add padding
    min_pos -= 5.0;
    max_pos += 5.0;

    //cout<<"Cube dimension = "<<min_pos<<"  "<<max_pos<<endl;
  
  
    fprintf(fp,"%5i%12.6f%12.6f%12.6f\n",syst.Number_of_atoms, min_pos.x, min_pos.y, min_pos.z);
  
    // Size of voxels
    VECTOR dr;
    dr.x = (max_pos.x - min_pos.x)/float(prms.nx_grid); 
    dr.y = (max_pos.y - min_pos.y)/float(prms.ny_grid); 
    dr.z = (max_pos.z - min_pos.z)/float(prms.nz_grid); 
    fprintf(fp,"%5i%12.6f%12.6f%12.6f\n",prms.nx_grid, dr.x,0.00,0.00);
    fprintf(fp,"%5i%12.6f%12.6f%12.6f\n",prms.ny_grid, 0.00,dr.y,0.00);
    fprintf(fp,"%5i%12.6f%12.6f%12.6f\n",prms.nz_grid, 0.00,0.00,dr.z);
  
  
  
    for(int n=0;n<syst.Number_of_atoms;n++){
      fprintf(fp,"%5i%12.6f%12.6f%12.6f%12.6f\n",syst.Atoms[n].Atom_Z, 0.0, syst.Atoms[n].Atom_RB.rb_cm.x, syst.Atoms[n].Atom_RB.rb_cm.y, syst.Atoms[n].Atom_RB.rb_cm.z); 
    }


  
    for(int ix=0;ix<prms.nx_grid;ix++){
      for(int iy=0;iy<prms.ny_grid;iy++){
        for(int iz=0;iz<prms.nz_grid;iz++){
  
           VECTOR pos;
           pos.x = min_pos.x + ix * dr.x;
           pos.y = min_pos.y + iy * dr.y;
           pos.z = min_pos.z + iz * dr.z;
           
           
           double psi = 0.0;
           // 
           //  | ADI_FMO_i> = sum C_Fi  | DIA_FMO_f >
           //                  f
           //
           // Here f = (fragment fr, i of fragment fr)

           int f = 0;

           for(int fr=0;fr<nfrags;fr++){ 
           for(int fr_i=0;fr_i<active_orb[fr].size();fr_i++){ 
           

             // 
             //  | DIA_FMO_f> = sum MO_af  | AO_a >
             //                  a
             //           

             double psi_f = 0.0;

             int Norb = ham[fr].el->Norb;
             for(int a=0;a<Norb;a++){

               psi_f += ham[fr].el->C_alp->M[a*Norb + active_orb[fr][fr_i]] * ham[fr].basis_ao[a].compute( pos );

             }// for 


             psi += C.M[f*nfmo + orb] * psi_f;

             f++;

           }// for fr_i - all orbitals in the fragment fr
           }// for fr - all fragments

           //cout<<ix<<"  "<<iy<<"  "<<iz<<"  "<<psi<<endl;
  
           double charg_dens = psi; // * psi;
           
           
           fprintf(fp, "%g ",charg_dens);  // charg_dens = data[ix][iy][iz]
  
          if (iz % 6 == 5)   fprintf(fp,"\n");
        }// for iz
        fprintf(fp,"\n");
      }// for iy
    }// for ix

    fclose(fp);

  }// for i

}// charge_density


void charge_density(MATRIX& C, boost::python::list ham, System& syst, boost::python::list active_orb, Control_Parameters& prms){

  vector<listHamiltonian_QM> ham_int;
  vector<vector<int> > active_orb_int;

  int sz = len(ham);
  for(int i=0;i<sz;i++){   ham_int.push_back( extract<listHamiltonian_QM>(ham[i]) );  }   // extract<int>(at.values()[i])

  sz = len(active_orb);  
  for(int fr=0;fr<sz;fr++){
    boost::python::list lst = extract<boost::python::list>(active_orb[fr]);

    vector<int> lst_i; 
    for(int i=0;i<len(lst);i++){  lst_i.push_back( extract<int>(lst[i]) ); }

    active_orb_int.push_back( lst_i );
  }

  charge_density(C, ham_int, syst, active_orb_int, prms);
  
}// charge_density


}// namespace libqchem_tools
}// liblibra
