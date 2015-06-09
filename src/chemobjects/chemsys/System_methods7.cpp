#include "System.h"


namespace libchemobjects{
namespace libchemsys{



void System::print_ent(std::string filename){
  print_ent(filename,1,"abc"); // by default - fold in 3D
}

void System::print_ent(std::string filename,int fold,std::string pbc_type){
// if fold==1 - will output coordinates folded into simulation box
// Print state of the system into file filename
// which has a Brookhaven PDB format (ENT)
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
    fprintf(fp,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f    \n",a,b,c,alp,bet,gam);
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

    fprintf(fp,"HETATM%5i %4s MOL C%4i    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",at_ser_num,at_name.c_str(),grp_num,r.x,r.y,r.z,occup,b_fact);
  }
  fclose(fp);
}

void System::print_ent(std::string filename,boost::python::list atoms_list){
  print_ent(filename,atoms_list,1,"abc"); // by default - fold in 3D
}

void System::print_ent(std::string filename,boost::python::list atoms_list,int fold,std::string pbc_type){
// Print state of the system into file filename
// which has a Brookhaven PDB format (ENT)
// Write only those atoms which IDs are in atomlist

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
    fprintf(fp,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f    \n",a,b,c,alp,bet,gam);
  }


  vector<int> at_indexes;
  int sz = len(atoms_list);
  for(int i=0;i<sz;i++){
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

    fprintf(fp,"HETATM%5i %4s MOL C%4i    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",at_ser_num,at_name.c_str(),grp_num,r.x,r.y,r.z,occup,b_fact);
  }
  fclose(fp);
}


//==================== XYZ format ==========================
void System::print_xyz(std::string filename,int frame){
  print_xyz(filename,1,"abc",frame); // by default - fold in 3D
}

void System::print_xyz(std::string filename,int fold,std::string pbc_type,int frame){
// if fold==1 - will output coordinates folded into simulation box
// Print state of the system into file filename
// which has a Brookhaven PDB format (ENT)
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
    fprintf(fp,"Molecule frame= %6d %9.3f%9.3f%9.3f%7.2f%7.2f%7.2f    \n",frame, a,b,c,alp,bet,gam);
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

    fprintf(fp," %4s  %8.3f%8.3f%8.3f   \n",at_name.c_str(),r.x,r.y,r.z);
  }
  fclose(fp);
}


}// namespace libchemsys
}// namespace libchemobjects




