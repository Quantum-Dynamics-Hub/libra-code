/*********************************************************************************
* Copyright (C) 2021 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Wfcgrid2_io.cpp
  \brief The file implements the methods for printing out/saving data from the Wfcgrid2 class 
    
*/

#include "Wfcgrid2.h"
#include "../../math_meigen/libmeigen.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

using namespace libwfcgrid;

/// libwfcgrid namespace
namespace libwfcgrid2{

using namespace libmeigen;



void Wfcgrid2::print_wfc_1D(std::string prefix, int rep, vector<int>& states, int do_real, int do_imag, int do_dens){
/**
  \brief Print a 1D wavefunction into 3 files: real part, imaginary part, and the probability density (if requested)
  \param[in] prefix The prefix of the filenames to which the wfc will be printed out
  \param[in] rep The wavefunction representation: 0 - diabatic, 1 - adiabatic 
  \param[in] states Indices of the states whose properties is to be printed out
  \param[in] do_real Print or not the info for the real part of the wfc
  \param[in] do_imag Print or not the info for the imaginary part of the wfc
  \param[in] do_dens Print or not the info for the probability density of the wfc

  for 1D profile on XY plane
*/

  int nx; 
  int Nx = npts[0];
  int st;
  int nst = states.size();

  std::string filename, reps;
  stringstream ss1(stringstream::in | stringstream::out);

  ss1 << rep;  ss1 >> reps;

  //====================== Real part ====================
  if(do_real){ 
    filename = prefix+"_real_rep_"+reps;
    ofstream out(filename.c_str(),ios::out);

    if(rep==0){
      for(nx=0;nx<Nx;nx++){
        out<<rgrid[0]->get(nx);
        for(st=0; st<nst; st++){   out<<"  "<<PSI_dia[nx].get(states[st],0).real(); }
        out<<endl;
      }// for nx
    }// rep ==0
    else if(rep==1){
      for(nx=0;nx<Nx;nx++){
        out<<rgrid[0]->get(nx);
        for(st=0; st<nst; st++){   out<<"  "<<PSI_adi[nx].get(states[st],0).real(); }
        out<<endl;
      }// for nx
    }// rep ==1

    out.close();
  }

  //====================== Imaginary part ====================
  if(do_imag){ 
    filename = prefix+"_imag_rep_"+reps;
    ofstream out(filename.c_str(),ios::out);

    if(rep==0){
      for(nx=0;nx<Nx;nx++){
        out<<rgrid[0]->get(nx);
        for(st=0; st<nst; st++){   out<<"  "<<PSI_dia[nx].get(states[st],0).imag(); }
        out<<endl;
      }// for nx
    }// rep ==0
    else if(rep==1){
      for(nx=0;nx<Nx;nx++){
        out<<rgrid[0]->get(nx);
        for(st=0; st<nst; st++){   out<<"  "<<PSI_adi[nx].get(states[st],0).imag(); }
        out<<endl;
      }// for nx
    }// rep ==1

    out.close();
  }

  //====================== Probability density ====================
  if(do_dens){ 
    filename = prefix+"_dens_rep_"+reps;
    ofstream out(filename.c_str(),ios::out);

    if(rep==0){
      for(nx=0;nx<Nx;nx++){
        out<<rgrid[0]->get(nx);
        for(st=0; st<nst; st++){   
          double dens = real(std::conj(PSI_dia[nx].get(states[st],0)) * PSI_dia[nx].get(states[st],0));
          out<<"  "<<dens;
        }
        out<<endl;
      }// for nx
    }// rep ==0
    else if(rep==1){
      for(nx=0;nx<Nx;nx++){
        out<<rgrid[0]->get(nx);
        for(st=0; st<nst; st++){
          double dens = real(std::conj(PSI_adi[nx].get(states[st],0)) * PSI_adi[nx].get(states[st],0));
          out<<"  "<<dens; 
        }
        out<<endl;
      }// for nx
    }// rep ==1

    out.close();
  }


}// print_wfc_1D



void Wfcgrid2::print_reci_wfc_1D(std::string prefix, int rep, vector<int>& states, int do_real, int do_imag, int do_dens){
/**
  \brief Print reciprocal representation of a 1D wavefunction into 3 files: 
  real part, imaginary part, and the probability density (if requested)
  \param[in] prefix The prefix of the filenames to which the wfc will be printed out
  \param[in] rep The wavefunction representation: 0 - diabatic, 1 - adiabatic 
  \param[in] states Indices of the states whose properties is to be printed out
  \param[in] do_real Print or not the info for the real part of the wfc
  \param[in] do_imag Print or not the info for the imaginary part of the wfc
  \param[in] do_dens Print or not the info for the probability density of the wfc

  for 1D profile on XY plane
*/

  int nx; 
  int Nx = npts[0];
  int st;
  int nst = states.size();

  std::string filename, reps;
  stringstream ss1(stringstream::in | stringstream::out);

  ss1 << rep;  ss1 >> reps;

  //====================== Real part ====================
  if(do_real){ 
    filename = prefix+"_real_rep_"+reps;
    ofstream out(filename.c_str(),ios::out);

    if(rep==0){
      for(nx=0;nx<Nx;nx++){
        out<<2.0*M_PI*kgrid[0]->get(nx);
        for(st=0; st<nst; st++){   out<<"  "<<reciPSI_dia[nx].get(states[st],0).real(); }       
        out<<endl;
      }// for nx
    }// rep ==0
    else if(rep==1){
      for(nx=0;nx<Nx;nx++){
        out<<2.0*M_PI*kgrid[0]->get(nx);
        for(st=0; st<nst; st++){   out<<"  "<<reciPSI_adi[nx].get(states[st],0).real(); }       
        out<<endl;
      }// for nx
    }// rep ==1

    out.close();
  }

  //====================== Imaginary part ====================
  if(do_imag){ 
    filename = prefix+"_imag_rep_"+reps;
    ofstream out(filename.c_str(),ios::out);

    if(rep==0){
      for(nx=0;nx<Nx;nx++){
        out<<2.0*M_PI*kgrid[0]->get(nx);
        for(st=0; st<nst; st++){   out<<"  "<<reciPSI_dia[nx].get(states[st],0).imag(); }       
        out<<endl;
      }// for nx
    }// rep ==0
    else if(rep==1){
      for(nx=0;nx<Nx;nx++){
        out<<2.0*M_PI*kgrid[0]->get(nx);
        for(st=0; st<nst; st++){   out<<"  "<<reciPSI_adi[nx].get(states[st],0).imag(); }       
        out<<endl;
      }// for nx
    }// rep ==1

    out.close();
  }

  //====================== Probability density ====================
  if(do_dens){ 
    filename = prefix+"_dens_rep_"+reps;
    ofstream out(filename.c_str(),ios::out);

    if(rep==0){
      for(nx=0;nx<Nx;nx++){
        out<<2.0*M_PI*kgrid[0]->get(nx);
        for(st=0; st<nst; st++){   
          double dens = real(std::conj(reciPSI_dia[nx].get(states[st],0)) * reciPSI_dia[nx].get(states[st],0));
          out<<"  "<<dens; 
        }       
        out<<endl;
      }// for nx
    }// rep ==0
    else if(rep==1){
      for(nx=0;nx<Nx;nx++){
        out<<2.0*M_PI*kgrid[0]->get(nx);        
        for(st=0; st<nst; st++){   
          double dens = real(std::conj(reciPSI_adi[nx].get(states[st],0)) * reciPSI_adi[nx].get(states[st],0));
          out<<"  "<<dens; 
        }       
        out<<endl;
      }// for nx
    }// rep ==1

    out.close();
  }


}// print_reci_wfc_1D





void Wfcgrid2::print_wfc_2D(std::string prefix, int rep, int state, int do_real, int do_imag, int do_dens){
/**
  \brief Print a 2D wavefunction into 3 files: real part, imaginary part, and the probability density (if requested)
  \param[in] prefix The prefix of the filenames to which the wfc will be printed out
  \param[in] rep The wavefunction representation: 0 - diabatic, 1 - adiabatic 
  \param[in] state Index of the state whose properties is to be printed out
  \param[in] do_real Print or not the info for the real part of the wfc
  \param[in] do_imag Print or not the info for the imaginary part of the wfc
  \param[in] do_dens Print or not the info for the probability density of the wfc

  for 2D profile on XY plane
*/

  int nx, ny, indx; 
  int Nx = npts[0];
  int Ny = npts[1];
  vector<int> ipt(2, 0);

  std::string filename, reps;
  stringstream ss1(stringstream::in | stringstream::out);

  ss1 << rep;  ss1 >> reps;

  //====================== Real part ====================
  if(do_real){ 
    filename = prefix+"_real_rep_"+reps;
    ofstream out(filename.c_str(),ios::out);

    if(rep==0){
      for(nx=0;nx<Nx;nx++){
        for(ny=0;ny<Ny;ny++){
          ipt[0] = nx; ipt[1] = ny;
          indx = imap(ipt);
          out<<rgrid[0]->get(nx)<<"  "<<rgrid[1]->get(ny)<<"  "<<PSI_dia[indx].get(state,0).real()<<endl;
        }// for ny
        out<<endl;
      }// for nx
    }// rep ==0

    else if(rep==1){
      for(nx=0;nx<Nx;nx++){
        for(ny=0;ny<Ny;ny++){
          ipt[0] = nx; ipt[1] = ny;
          indx = imap(ipt);
          out<<rgrid[0]->get(nx)<<"  "<<rgrid[1]->get(ny)<<"  "<<PSI_adi[indx].get(state,0).real()<<endl;
        }// for ny
        out<<endl;
      }// for nx

    }// rep ==1

    out.close();
  }

  //====================== Imaginary part ====================
  if(do_imag){ 
    filename = prefix+"_imag_rep_"+reps;
    ofstream out(filename.c_str(),ios::out);

    if(rep==0){
      for(nx=0;nx<Nx;nx++){
        for(ny=0;ny<Ny;ny++){
          ipt[0] = nx; ipt[1] = ny;
          indx = imap(ipt);
          out<<rgrid[0]->get(nx)<<"  "<<rgrid[1]->get(ny)<<"  "<<PSI_dia[indx].get(state,0).imag()<<endl;
        }// for ny
        out<<endl;
      }// for nx
    }// rep ==0

    else if(rep==1){
      for(nx=0;nx<Nx;nx++){
        for(ny=0;ny<Ny;ny++){
          ipt[0] = nx; ipt[1] = ny;
          indx = imap(ipt);
          out<<rgrid[0]->get(nx)<<"  "<<rgrid[1]->get(ny)<<"  "<<PSI_adi[indx].get(state,0).imag()<<endl;
        }// for ny
        out<<endl;
      }// for nx

    }// rep ==1

    out.close();
  }

  //====================== Probability density ====================
  if(do_dens){ 
    filename = prefix+"_dens_rep_"+reps;
    ofstream out(filename.c_str(),ios::out);

    if(rep==0){
      for(nx=0;nx<Nx;nx++){
        for(ny=0;ny<Ny;ny++){
          ipt[0] = nx; ipt[1] = ny;
          indx = imap(ipt);
          double dens = (std::conj(PSI_dia[indx].get(state,0)) * PSI_dia[indx].get(state,0) ).real();
          out<<rgrid[0]->get(nx)<<"  "<<rgrid[1]->get(ny)<<"  "<<dens<<endl;
        }// for ny
        out<<endl;
      }// for nx
    }// rep ==0

    else if(rep==1){
      for(nx=0;nx<Nx;nx++){
        for(ny=0;ny<Ny;ny++){
          ipt[0] = nx; ipt[1] = ny;
          indx = imap(ipt);
          double dens = (std::conj(PSI_adi[indx].get(state,0)) * PSI_adi[indx].get(state,0) ).real();
          out<<rgrid[0]->get(nx)<<"  "<<rgrid[1]->get(ny)<<"  "<<dens<<endl;
        }// for ny
        out<<endl;
      }// for nx

    }// rep ==1

    out.close();
  }


}// print_wfc_2D




void Wfcgrid2::print_reci_wfc_2D(std::string prefix, int rep, int state, int do_real, int do_imag, int do_dens){
/**
  \brief Print a reciprocal-space 2D wavefunction into 3 files: real part, imaginary part, and the probability density (if requested)
  \param[in] prefix The prefix of the filenames to which the wfc will be printed out
  \param[in] rep The wavefunction representation: 0 - diabatic, 1 - adiabatic 
  \param[in] state Index of the state whose properties is to be printed out
  \param[in] do_real Print or not the info for the real part of the wfc
  \param[in] do_imag Print or not the info for the imaginary part of the wfc
  \param[in] do_dens Print or not the info for the probability density of the wfc

  for 2D profile on XY plane
*/

  int nx, ny, indx; 
  int Nx = npts[0];
  int Ny = npts[1];
  vector<int> ipt(2, 0);

  std::string filename, reps;
  stringstream ss1(stringstream::in | stringstream::out);

  ss1 << rep;  ss1 >> reps;

  //====================== Real part ====================
  if(do_real){ 
    filename = prefix+"_real_rep_"+reps;
    ofstream out(filename.c_str(),ios::out);

    if(rep==0){
      for(nx=0;nx<Nx;nx++){
        for(ny=0;ny<Ny;ny++){
          ipt[0] = nx; ipt[1] = ny;
          indx = imap(ipt);
          out<<2.0*M_PI*kgrid[0]->get(nx)<<"  "<<2.0*M_PI*kgrid[1]->get(ny)<<"  "<<reciPSI_dia[indx].get(state,0).real()<<endl;
        }// for ny
        out<<endl;
      }// for nx
    }// rep ==0

    else if(rep==1){
      for(nx=0;nx<Nx;nx++){
        for(ny=0;ny<Ny;ny++){
          ipt[0] = nx; ipt[1] = ny;
          indx = imap(ipt);
          out<<2.0*M_PI*kgrid[0]->get(nx)<<"  "<<2.0*M_PI*kgrid[1]->get(ny)<<"  "<<reciPSI_adi[indx].get(state,0).real()<<endl;
        }// for ny
        out<<endl;
      }// for nx

    }// rep ==1

    out.close();
  }

  //====================== Imaginary part ====================
  if(do_imag){ 
    filename = prefix+"_imag_rep_"+reps;
    ofstream out(filename.c_str(),ios::out);

    if(rep==0){
      for(nx=0;nx<Nx;nx++){
        for(ny=0;ny<Ny;ny++){
          ipt[0] = nx; ipt[1] = ny;
          indx = imap(ipt);
          out<<2.0*M_PI*kgrid[0]->get(nx)<<"  "<<2.0*M_PI*kgrid[1]->get(ny)<<"  "<<reciPSI_dia[indx].get(state,0).imag()<<endl;
        }// for ny
        out<<endl;
      }// for nx
    }// rep ==0

    else if(rep==1){
      for(nx=0;nx<Nx;nx++){
        for(ny=0;ny<Ny;ny++){
          ipt[0] = nx; ipt[1] = ny;
          indx = imap(ipt);
          out<<2.0*M_PI*kgrid[0]->get(nx)<<"  "<<2.0*M_PI*kgrid[1]->get(ny)<<"  "<<reciPSI_adi[indx].get(state,0).imag()<<endl;
        }// for ny
        out<<endl;
      }// for nx

    }// rep ==1

    out.close();
  }

  //====================== Probability density ====================
  if(do_dens){ 
    filename = prefix+"_dens_rep_"+reps;
    ofstream out(filename.c_str(),ios::out);

    if(rep==0){
      for(nx=0;nx<Nx;nx++){
        for(ny=0;ny<Ny;ny++){
          ipt[0] = nx; ipt[1] = ny;
          indx = imap(ipt);
          double dens = (std::conj(reciPSI_dia[indx].get(state,0)) * reciPSI_dia[indx].get(state,0) ).real();
          out<<2.0*M_PI*kgrid[0]->get(nx)<<"  "<<2.0*M_PI*kgrid[1]->get(ny)<<"  "<<dens<<endl;
        }// for ny
        out<<endl;
      }// for nx
    }// rep ==0

    else if(rep==1){
      for(nx=0;nx<Nx;nx++){
        for(ny=0;ny<Ny;ny++){
          ipt[0] = nx; ipt[1] = ny;
          indx = imap(ipt);
          double dens = (std::conj(reciPSI_adi[indx].get(state,0)) * reciPSI_adi[indx].get(state,0) ).real();
          out<<2.0*M_PI*kgrid[0]->get(nx)<<"  "<<2.0*M_PI*kgrid[1]->get(ny)<<"  "<<dens<<endl;
        }// for ny
        out<<endl;
      }// for nx

    }// rep ==1

    out.close();
  }


}// print_reci_wfc_2D





}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra

