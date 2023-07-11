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
  \file System_methods7.cpp
  \brief This file implements formatted printing and output of the System's properties
    
*/

#include <sstream>

#include "System.h"
#include "../../Units.h"

/// liblibra namespace
namespace liblibra{

/// libchemobjects namespace
namespace libchemobjects{

/// libchemsys namespace
namespace libchemsys{



void System::print_ent(std::string filename){
/** 
  \param[in] filename The name of the file where the info will be trinted out
   
  Print the state of the system into file in Brookhaven PDB format (ENT). This version will fold atomic coordinates into 3D box
*/

  print_ent(filename,1,"abc"); // by default - fold in 3D
}

void System::print_ent(std::string filename,int fold,std::string pbc_type){
/** 
  \param[in] filename The name of the file where the info will be trinted out
  \param[in] fold Controlls the folding of the coordinates into the unit-cell 
  if fold==1 - will output coordinates folded into simulation box
  \param[in] pbc_type The parameter controlling the periodicity (when and if folding) of the unit cell
  Can take values: "a", "b", "c", "ab", "ac", "bc", and "abc"

  Print the state of the system into file in Brookhaven PDB format (ENT)  
*/

  FILE* fp;
  fp = fopen(filename.c_str(),"w");

// Crystal structure information, CRYST1 record see: 
//http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
  if(is_Box){
    VECTOR tv1,tv2,tv3;
    Box.get_vectors(tv1,tv2,tv3);
    double a,b,c,alp,bet,gam;
    a = tv1.length();
    b = tv2.length();
    c = tv3.length();
    tv1 = tv1/a;
    tv2 = tv2/b;
    tv3 = tv3/c;
    alp = acos(tv2*tv3)*radians_to_degrees;
    bet = acos(tv3*tv1)*radians_to_degrees;
    gam = acos(tv1*tv2)*radians_to_degrees;
    fprintf(fp,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f    \n",a * bohr, b * bohr, c * bohr,alp,bet,gam); // convert from internal units
     // (Bohrs) to Angstroms
  }

  for(int i=0;i<Number_of_atoms;i++){
    int at_ser_num = 1;
    std::string at_name = "H";
    int grp_num;
    VECTOR r;
    double occup = 1.0;
    double b_fact = 0.0;

    if(Atoms[i].is_Atom_id){at_ser_num = Atoms[i].Atom_id;}
    if(Atoms[i].is_Atom_element){at_name = Atoms[i].Atom_element;}
    grp_num = Atoms[i].globGroup_Index;
    if(Atoms[i].is_Atom_RB){  if(Atoms[i].Atom_RB.is_rb_cm){   r = Atoms[i].Atom_RB.rb_cm;  }  }
    if(Atoms[i].is_Atom_charge){ b_fact = Atoms[i].Atom_charge; }

    if(fold && is_Box){
      MATRIX3x3 invBox; invBox = Box.inverse();
      VECTOR borig; borig = 0.0; // box origin

      r = invBox*(r - borig);

      if((pbc_type=="a")||(pbc_type=="ab")||(pbc_type=="ac")||(pbc_type=="abc")) {
        r.x = r.x - floor(r.x);
      }
      if((pbc_type=="b")||(pbc_type=="ab")||(pbc_type=="bc")||(pbc_type=="abc")) {
        r.y = r.y - floor(r.y);
      }
      if((pbc_type=="c")||(pbc_type=="ac")||(pbc_type=="bc")||(pbc_type=="abc")) {
        r.z = r.z - floor(r.z);
      }
      r = Box * r + borig;
    }// fold && is_Box

    r *= bohr; // convert from internal units (a.u.) to conventional units (Angstrom)

    fprintf(fp,"HETATM%5i %4s MOL C%4i    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",at_ser_num,at_name.c_str(),grp_num,r.x,r.y,r.z,occup,b_fact);
  }
  fclose(fp);
}

void System::print_ent(std::string filename,boost::python::list atoms_list){
/** 
  \param[in] filename The name of the file where the info will be trinted out
  \param[in] atoms_list The list of the atom IDs for the atoms to be printed out 

  Print the state of the system into file in Brookhaven PDB format (ENT). Write only specified atoms which IDs are in the atomlist
  This version will fold atomic coordinates into 3D box
*/

  print_ent(filename,atoms_list,1,"abc"); // by default - fold in 3D
}

void System::print_ent(std::string filename,boost::python::list atoms_list,int fold,std::string pbc_type){
/** 
  \param[in] filename The name of the file where the info will be trinted out
  \param[in] atoms_list The list of the atom IDs for the atoms to be printed out 
  \param[in] fold Controlls the folding of the coordinates into the unit-cell 
  if fold==1 - will output coordinates folded into simulation box
  \param[in] pbc_type The parameter controlling the periodicity (when and if folding) of the unit cell
  Can take values: "a", "b", "c", "ab", "ac", "bc", and "abc"

  Print the state of the system into file in Brookhaven PDB format (ENT). Write only specified atoms which IDs are in the atomlist
*/
  int i;
  FILE* fp;
  fp = fopen(filename.c_str(),"w");

// Crystal structure information, CRYST1 record see:
//http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
  if(is_Box){
    VECTOR tv1,tv2,tv3;
    Box.get_vectors(tv1,tv2,tv3);
    double a,b,c,alp,bet,gam;
    a = tv1.length();
    b = tv2.length();
    c = tv3.length();
    tv1 = tv1/a;
    tv2 = tv2/b;
    tv3 = tv3/c;
    alp = acos(tv2*tv3)*radians_to_degrees;
    bet = acos(tv3*tv1)*radians_to_degrees;
    gam = acos(tv1*tv2)*radians_to_degrees;
    fprintf(fp,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f    \n",a * bohr, b * bohr, c * bohr, alp,bet,gam); // convert from internal units
     // (Bohrs) to Angstroms

  }


  vector<int> at_indexes;
  int sz = len(atoms_list);
  for(i=0;i<sz;i++){
    int id = extract<int>(atoms_list[i]);
    int j = get_atom_index_by_atom_id(id);
    if(j!=-1){ at_indexes.push_back(j); }
  }

  sz = at_indexes.size();
  for(i=0;i<sz;i++){
    int indx = at_indexes[i];
    int at_ser_num = 1;
    std::string at_name = "H";
    int grp_num;
    VECTOR r;
    double occup = 1.0;
    double b_fact = 0.0;

    if(Atoms[indx].is_Atom_id){at_ser_num = Atoms[indx].Atom_id;}
    if(Atoms[indx].is_Atom_element){at_name = Atoms[indx].Atom_element;}
    grp_num = Atoms[indx].globGroup_Index;
    if(Atoms[indx].is_Atom_RB){  if(Atoms[indx].Atom_RB.is_rb_cm){   r = Atoms[indx].Atom_RB.rb_cm;  }  }
    if(Atoms[indx].is_Atom_charge){ b_fact = Atoms[indx].Atom_charge; }

    if(fold && is_Box){
      MATRIX3x3 invBox; invBox = Box.inverse();
      VECTOR borig; borig = 0.0; // box origin

      r = invBox*(r - borig);

      if((pbc_type=="a")||(pbc_type=="ab")||(pbc_type=="ac")||(pbc_type=="abc")) {
        r.x = r.x - floor(r.x);
      }
      if((pbc_type=="b")||(pbc_type=="ab")||(pbc_type=="bc")||(pbc_type=="abc")) {
        r.y = r.y - floor(r.y);
      }
      if((pbc_type=="c")||(pbc_type=="ac")||(pbc_type=="bc")||(pbc_type=="abc")) {
        r.z = r.z - floor(r.z);
      }
      r = Box * r + borig;
    }// fold && is_Box

    r *= bohr; // convert from internal units (a.u.) to conventional units (Angstrom)

    fprintf(fp,"HETATM%5i %4s MOL C%4i    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",at_ser_num,at_name.c_str(),grp_num,r.x,r.y,r.z,occup,b_fact);
  }
  fclose(fp);
}


//==================== XYZ format ==========================
void System::print_xyz(std::string filename,int frame){
/** 
  \param[in] filename The name of the file where the info will be trinted out
  \param[in] frame The integer index (in MD - the MD step, for instance) to be printed in the label of each record

  Print the state of the system into file in XYZ format. Note that this functions will append the output to the
  file fileanme, if the latter is not empty. This way, we can print out the whole MD trajectory in a single file. 
  If this is not what you want, you may need to delete the older version of the file
  This version will fold atomic coordinates into 3D box
*/

  print_xyz(filename,1,"abc",frame); // by default - fold in 3D
}

void System::print_xyz(std::string filename,int fold,std::string pbc_type,int frame){
/** 
  \param[in] filename The name of the file where the info will be trinted out
  \param[in] fold Controlls the folding of the coordinates into the unit-cell 
  if fold==1 - will output coordinates folded into simulation box
  \param[in] pbc_type The parameter controlling the periodicity (when and if folding) of the unit cell
  Can take values: "a", "b", "c", "ab", "ac", "bc", and "abc"
  \param[in] frame The integer index (in MD - the MD step, for instance) to be printed in the label of each record

  Print the state of the system into file in XYZ format. Note that this functions will append the output to the
  file fileanme, if the latter is not empty. This way, we can print out the whole MD trajectory in a single file. 
  If this is not what you want, you may need to delete the older version of the file
*/

  FILE* fp;
  fp = fopen(filename.c_str(),"a");

  fprintf(fp,"%6d  \n",Number_of_atoms);

// Crystal structure information, CRYST1 record see: 
//http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
  if(is_Box){
    VECTOR tv1,tv2,tv3;
    Box.get_vectors(tv1,tv2,tv3);
    double a,b,c,alp,bet,gam;
    a = tv1.length();
    b = tv2.length();
    c = tv3.length();
    tv1 = tv1/a;
    tv2 = tv2/b;
    tv3 = tv3/c;
    alp = acos(tv2*tv3)*radians_to_degrees;
    bet = acos(tv3*tv1)*radians_to_degrees;
    gam = acos(tv1*tv2)*radians_to_degrees;
    fprintf(fp,"Molecule frame= %6d %9.3f%9.3f%9.3f%7.2f%7.2f%7.2f    \n",frame, a * bohr,b * bohr,c * bohr,alp,bet,gam);
  }
  else{
    fprintf(fp,"Molecule frame= %6d \n", frame);
  }

  for(int i=0;i<Number_of_atoms;i++){
    int at_ser_num = 1;
    std::string at_name = "H";
    int grp_num;
    VECTOR r;
    double occup = 1.0;
    double b_fact = 0.0;

    if(Atoms[i].is_Atom_id){at_ser_num = Atoms[i].Atom_id;}
    if(Atoms[i].is_Atom_element){at_name = Atoms[i].Atom_element;}
    grp_num = Atoms[i].globGroup_Index;
    if(Atoms[i].is_Atom_RB){  if(Atoms[i].Atom_RB.is_rb_cm){   r = Atoms[i].Atom_RB.rb_cm;  }  }
    if(Atoms[i].is_Atom_charge){ b_fact = Atoms[i].Atom_charge; }

    if(fold && is_Box){
      MATRIX3x3 invBox; invBox = Box.inverse();
      VECTOR borig; borig = 0.0; // box origin

      r = invBox*(r - borig);

      if((pbc_type=="a")||(pbc_type=="ab")||(pbc_type=="ac")||(pbc_type=="abc")) {
        r.x = r.x - floor(r.x);
      }
      if((pbc_type=="b")||(pbc_type=="ab")||(pbc_type=="bc")||(pbc_type=="abc")) {
        r.y = r.y - floor(r.y);
      }
      if((pbc_type=="c")||(pbc_type=="ac")||(pbc_type=="bc")||(pbc_type=="abc")) {
        r.z = r.z - floor(r.z);
      }
      r = Box * r + borig;
    }// fold && is_Box

    r *= bohr; // convert from internal units (a.u.) to conventional units (Angstrom)

    fprintf(fp," %4s  %8.3f%8.3f%8.3f   \n",at_name.c_str(),r.x,r.y,r.z);
  }
  fclose(fp);
}


std::string System::get_xyz(int fold,std::string pbc_type,int frame){
/**Returns a string representing the system in an xyz format 
  
  \param[in] filename The name of the file where the info will be trinted out
  \param[in] fold Controlls the folding of the coordinates into the unit-cell 
  if fold==1 - will output coordinates folded into simulation box
  \param[in] pbc_type The parameter controlling the periodicity (when and if folding) of the unit cell
  Can take values: "a", "b", "c", "ab", "ac", "bc", and "abc"
  \param[in] frame The integer index (in MD - the MD step, for instance) to be printed in the label of each record

  Print the state of the system into file in XYZ format. Note that this functions will append the output to the
  file fileanme, if the latter is not empty. This way, we can print out the whole MD trajectory in a single file. 
  If this is not what you want, you may need to delete the older version of the file
*/

  stringstream ss; //(stringstream::in | stringstream::out);
  //std::string out;

  ss << Number_of_atoms << "\n";  // %6d

  // Crystal structure information, CRYST1 record see: 
  //http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
  if(is_Box){
    VECTOR tv1,tv2,tv3;
    Box.get_vectors(tv1,tv2,tv3);
    double a,b,c,alp,bet,gam;
    a = tv1.length();
    b = tv2.length();
    c = tv3.length();
    tv1 = tv1/a;
    tv2 = tv2/b;
    tv3 = tv3/c;
    alp = acos(tv2*tv3)*radians_to_degrees;
    bet = acos(tv3*tv1)*radians_to_degrees;
    gam = acos(tv1*tv2)*radians_to_degrees;

    ss << "Molecule frame= "<<frame<<" "<<a*bohr<<" "<<b*bohr<<" "<<c*bohr; 
    ss <<" "<<alp<<" "<<bet<<" "<<gam <<"\n";   // %6d %9.3f%9.3f%9.3f%7.2f%7.2f%7.2f
    //ss>>out;
  } 
  else{
    ss << "Molecule frame= "<<frame<<" \n"; // %6d 
    //ss>>out;
  }

  for(int i=0;i<Number_of_atoms;i++){
    int at_ser_num = 1;
    std::string at_name = "H";
    int grp_num;
    VECTOR r;
    double occup = 1.0;
    double b_fact = 0.0;

    if(Atoms[i].is_Atom_id){at_ser_num = Atoms[i].Atom_id;}
    if(Atoms[i].is_Atom_element){at_name = Atoms[i].Atom_element;}
    grp_num = Atoms[i].globGroup_Index;
    if(Atoms[i].is_Atom_RB){  if(Atoms[i].Atom_RB.is_rb_cm){   r = Atoms[i].Atom_RB.rb_cm;  }  }
    if(Atoms[i].is_Atom_charge){ b_fact = Atoms[i].Atom_charge; }

    if(fold && is_Box){
      MATRIX3x3 invBox; invBox = Box.inverse();
      VECTOR borig; borig = 0.0; // box origin

      r = invBox * (r - borig);

      if((pbc_type=="a")||(pbc_type=="ab")||(pbc_type=="ac")||(pbc_type=="abc")) {
        r.x = r.x - floor(r.x);
      }
      if((pbc_type=="b")||(pbc_type=="ab")||(pbc_type=="bc")||(pbc_type=="abc")) {
        r.y = r.y - floor(r.y);
      }
      if((pbc_type=="c")||(pbc_type=="ac")||(pbc_type=="bc")||(pbc_type=="abc")) {
        r.z = r.z - floor(r.z);
      }
      r = Box * r + borig;
    }// fold && is_Box

    r *= bohr; // convert from internal units (a.u.) to conventional units (Angstrom)

    ss << at_name<<" "<<r.x<<" "<<r.y<<" "<<r.z<<"\n";  // %4s  %8.3f%8.3f%8.3f
    //ss>>out;
  }

  return ss.str();
}



}// namespace libchemsys
}// namespace libchemobjects
}// liblibra



