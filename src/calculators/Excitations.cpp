#include "Excitations.h"

namespace libcalculators{


void excite(int I, int J, vector< pair<int,double> >& occ_ini, vector< pair<int,double> >& occ_fin){
  // Create I --> J ecxitation
  // basically, vector< pair<int,double> > datastructure represents Slater determinant (by selecting ordering of orbitals in the active space)

  int Norb = occ_ini.size();

  if(I<0){ std::cout<<"Error: source orbital index("<<I<<") can not be negative\n";  exit(0); }
  if(J<0){ std::cout<<"Error: target orbital index("<<J<<") can not be negative\n";  exit(0); }
  if(I>=Norb){ std::cout<<"Error: source orbital index("<<I<<") exceeds the number of orbitals in the active space("<<Norb<<")\n";  exit(0); }
  if(J>=Norb){ std::cout<<"Error: target orbital index("<<J<<") exceeds the number of orbitals in the active space("<<Norb<<")\n";  exit(0); }


  if(occ_fin.size()>0){  occ_fin.clear(); }                    // clean result

  for(int i=0;i<Norb;i++){  occ_fin.push_back(occ_ini[i]);  }  // copy initial to final

  // Swap populations I-th and J-th orbitals
  double pop = occ_fin[I].second; 
  occ_fin[I].second = occ_fin[J].second;
  occ_fin[J].second = pop;
   
}


}// namespace libcalculators
