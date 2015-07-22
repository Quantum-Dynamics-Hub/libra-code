#include "Hamiltonian.h"

/****************************************************************************

  This file contains following functions:

  void Hamiltonian_core(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                        vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                        MATRIX* Hao, MATRIX* Sao, Memory* mem)


  void Hamiltonian_Fock(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                           vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                           Electronic* el,Electronic* el0, Memory* mem)


****************************************************************************/

void Hamiltonian_core(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                      vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                      MATRIX* Hao, MATRIX* Sao, Memory* mem){

  if(prms.hamiltonian=="hf"){

    Hamiltonian_core_hf(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, Hao, Sao, mem);

  }
  else if(prms.hamiltonian=="eht"){

    Hamiltonian_core_eht(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, Hao, Sao, mem);

  }
  else if(prms.hamiltonian=="geht"){

    Hamiltonian_core_geht(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, Hao, Sao, mem);

  }
  else if(prms.hamiltonian=="geht1"){

    Hamiltonian_core_geht1(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, Hao, Sao, mem);
//    Hamiltonian_core_geht1(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, Hao, Sao, mem, mem->eri, mem->V_AB);

  }
  else if(prms.hamiltonian=="geht2"){

    Hamiltonian_core_geht2(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, Hao, Sao, mem);

  }
  else if(prms.hamiltonian=="indo"){

    Hamiltonian_core_indo(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, Hao, Sao, mem, mem->eri, mem->V_AB);

  }


}


void Hamiltonian_Fock(Control_Parameters& prms,Model_Parameters& modprms,Nuclear& mol,
                      vector<int>& fragment, vector<int>& basis_fo,vector<AO>& basis_ao,vector<vector<int> >& at_orbitals,
                      Electronic* el,Electronic* el0, Memory* mem){

  if(prms.hamiltonian=="hf"){

    Hamiltonian_Fock_hf(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, el, el0, mem);

  }
  else if(prms.hamiltonian=="eht"){

    Hamiltonian_Fock_eht(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, el, el0, mem);

  }
  else if(prms.hamiltonian=="geht"){

    Hamiltonian_Fock_geht(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, el, el0, mem);

  }
  else if(prms.hamiltonian=="geht1"){

    Hamiltonian_Fock_geht1_v1(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, el, el0, mem);
//    Hamiltonian_Fock_geht1(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, el, el0, mem);
//    Hamiltonian_Fock_geht1(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, el, mem, mem->eri, mem->V_AB);

  }
  else if(prms.hamiltonian=="geht2"){

    Hamiltonian_Fock_geht2(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, el, el0, mem);

  }

  else if(prms.hamiltonian=="indo"){

    Hamiltonian_Fock_indo(prms, modprms, mol, fragment, basis_fo, basis_ao, at_orbitals, el, mem, mem->eri, mem->V_AB);

  }



}

