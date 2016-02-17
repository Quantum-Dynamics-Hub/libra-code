/*********************************************************************************
* Copyright (C) 2016 Alexey V. Akimov
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

#include "Charge_Density.h"

#include <sstream>
using namespace std;

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



}// namespace libqchem_tools

