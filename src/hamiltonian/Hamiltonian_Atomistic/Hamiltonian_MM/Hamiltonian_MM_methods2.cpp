#include "Hamiltonian_MM.h"


namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_mm{


int is_in_vector(int indx,vector<int> vect){
/*************************************************************************
   This functions is not to be exposed to user.
   It searches for index "indx" of int type in a vector of ints
   Return 1 if "indx" has been found, 0 - otherwise
**************************************************************************/
  int res = 0;
  int sz = vect.size();
  for(int i=0;i<sz;i++){
    if(vect[i]==indx) { res = 1; break;}
  }
  return res;
}



//========================================================
// Setting up interactions in listHamiltonian_MM class
//========================================================

int listHamiltonian_MM::is_new_interaction(Hamiltonian_MM& inter){
  int res = 1; // Assume it is new
/*  This verification is time-consuming (~N^4) so for big systems we do not do this!
  int sz = interactions.size();
  for(int i=0;i<sz;i++){
    if(inter==interactions[i]){ res = 0; }
  }
*/
  return res;
}

void listHamiltonian_MM::show_interactions_statistics(){
  int sz = interactions.size();
  cout<<"Total number of interactions is "<<sz<<endl;
  int b,a,d,o,v,e,mb,cg, mbe;
  int ab,aa,ad,ao,av,ae,amb,acg, ambe;
  b = a = d = o = v = e = mb = cg = mbe = 0;
  ab = aa = ad = ao = av = ae = amb = acg = ambe = 0;

  for(int i=0;i<sz;i++){
    int t = interactions[i].get_type();
    int is_active = interactions[i].get_status();
    if(t==0){ b++; ab += is_active;}
    else if(t==1){ a++; aa += is_active; }
    else if(t==2){ d++; ad += is_active; }
    else if(t==3){ o++; ao += is_active; }
    else if(t==4){ v++; av += is_active; }
    else if(t==5){ e++;ae += is_active; }
    else if(t==6){ mb++; amb += is_active; }
    else if(t==7){ cg++; acg += is_active; }
    else if(t==8){ mbe++; ambe += is_active; }
  }
  std::cout<<"Number of bond interactions = "<<b<<"\n";
  std::cout<<"-------------------- active = "<<ab<<"\n";
  std::cout<<"Number of angle interactions = "<<a<<"\n";
  std::cout<<"--------------------- active = "<<aa<<"\n";
  std::cout<<"Number of dihedral interactions = "<<d<<"\n";
  std::cout<<"------------------------ active = "<<ad<<"\n";
  std::cout<<"Number of oop interactions = "<<o<<"\n";
  std::cout<<"------------------- active = "<<ao<<"\n";
  std::cout<<"Number of vdw interactions = "<<v<<"\n";
  std::cout<<"------------------- active = "<<av<<"\n";
  std::cout<<"Number of elec interactions = "<<e<<"\n";
  std::cout<<"-------------------- active = "<<ae<<"\n";
  std::cout<<"Number of mb interactions = "<<mb<<"\n";
  std::cout<<"------------------ active = "<<amb<<"\n";
  std::cout<<"Number of cg interactions = "<<cg<<"\n";
  std::cout<<"------------------ active = "<<acg<<"\n";
  std::cout<<"Number of mb_excl interactions = "<<mbe<<"\n";
  std::cout<<"----------------------- active = "<<ambe<<"\n";

}

void listHamiltonian_MM::set_atom_types(System& syst, vector<int>& lst,ForceField& ff){
/********************************************************************
  This function checks if the atom type of the atoms in the list are
  valid for a given force field.
  If they are valid and defined - no action is required, otherwise the
  typization will be performed according to the rules of the given
  force field
*********************************************************************/
  cout<<"listHamiltonian_MM::set_atom_types\n";
  cout<<"Setting up atom types for the force field "<< ff.ForceField_Name<<endl;

  int sz = lst.size();
  for(int i=0;i<sz;i++){
      int at = lst[i];
      int _defined = 0;
      if(syst.Atoms[at].is_Atom_ff_type){ 
        _defined = ff.is_valid_atom_type(syst.Atoms[at].Atom_ff_type);  

        cout<<"Atom i="<<i<<" type is defined manually: "<<syst.Atoms[at].Atom_ff_type<<endl;
      }
      if(!_defined){
        string ff_type = ff.get_atom_type(syst.Atoms[at].Atom_element,
                                          syst.Atoms[at].globAtom_Adjacent_Atoms.size(),
                                          syst.Atoms[at].Atom_functional_group,
                                          syst.Atoms[at].Atom_min_ring_size,
                                          syst.Atoms[at].Atom_coordination);
        // Since we have calculated some atomic properties we must declare that they are defined
        syst.Atoms[at].is_Atom_coordination = 1;
        if(ff_type!="none"){
          syst.Atoms[at].Atom_ff_type = ff_type;
          syst.Atoms[at].is_Atom_ff_type = 1;
        }
        cout<<"Atom i="<<i<<" type is defined automatically: "<<syst.Atoms[at].Atom_ff_type<<endl;
      }// if !defined
  }// for i
}

void listHamiltonian_MM::set_fragment_types(System& syst, vector<int>& lst,ForceField& ff){
/********************************************************************
  This function checks if the atom type of the atoms in the list are
  valid for a given force field.
  If they are valid and defined - no action is required, otherwise the
  typization will be performed according to the rules of the given
  force field
*********************************************************************/
  int sz = lst.size();
  for(int i=0;i<sz;i++){

      int fr = lst[i];
      int _defined = 0;
      if(syst.Fragments[fr].is_Group_ff_type){ _defined = ff.is_valid_fragment_type(syst.Fragments[fr].Group_ff_type);   }
      if(!_defined){
/*
        string ff_type = ff.get_atom_type(Atoms[at].Atom_element,
                                          Atoms[at].globAtom_Adjacent_Atoms.size(),
                                          Atoms[at].Atom_functional_group,
                                          Atoms[at].Atom_min_ring_size,
                                          Atoms[at].Atom_coordination);
        // Since we have calculated some atomic properties we must declare that they are defined
        Atoms[at].is_Atom_coordination = 1;
        if(ff_type!="none"){
          Atoms[at].Atom_ff_type = ff_type;
          Atoms[at].is_Atom_ff_type = 1;
        }
*/
      }// if !defined

  }// for i
}


bool listHamiltonian_MM::is_active(Atom& at1,Atom& at2){
  return !(at1.globGroup_Index==at2.globGroup_Index);
}
bool listHamiltonian_MM::is_active(Atom& at1,Atom& at2,Atom& at3){
  return !((at1.globGroup_Index==at2.globGroup_Index) &&
           (at2.globGroup_Index==at3.globGroup_Index)
          );
}
bool listHamiltonian_MM::is_active(Atom& at1,Atom& at2,Atom& at3,Atom& at4){
  return !((at1.globGroup_Index==at2.globGroup_Index) &&
           (at2.globGroup_Index==at3.globGroup_Index) &&
           (at3.globGroup_Index==at4.globGroup_Index)
          );
}


void listHamiltonian_MM::set_atom_interactions_for_atoms
(System& syst, string int_type,vector<Atom>& top_elt,vector<int>& lst1,vector<int>& lst2,ForceField& ff,int verb){
/**************************************************************
  This function determines the parameters for many-body interactions
  corresponding to the force field and the functional form of the
  potential in use. It then creates the interaction objects
 (interactions) for all atoms defined in the system.
  If the interactions were already defined - they are not included.
  Also if not all necessary parameters obtained - the interaction
  is not created.
  <verb> determines verbosity level
***************************************************************/
  int sz = top_elt.size();
  int* at; at = new int[sz];
  int* id; id = new int[sz];
  bool x = true;
  for(int j=0;j<sz;j++){
    at[j] = top_elt[j].globAtom_Index;
    id[j] = top_elt[j].Atom_id;
    x = (x && (is_in_vector(at[j],lst1)||is_in_vector(at[j],lst2)));
  }
  // If now x is true - then each of the atoms of the topological element
  // belongs to one of the lists
  if(x){

  //--------------- Get the interaction parameters from the force field --------------
  map<string,double> prms;
  int res = 1;

  if(int_type=="mb"){  res = ff.get_mb_parameters(prms); }
  else if(int_type=="mb_excl") { res = ff.get_mb_excl_parameters(prms); }

  if(verb>0){
    cout<<"In set_atom_interactions_for_atoms, int_type = "<<int_type<<endl;
    cout<<"In set_atom_interactions_for_atoms, mb_functional = "<<ff.mb_functional<<endl;
    cout<<"In set_atom_interactions_for_atoms, res after get_mb_parameters or get_mb_excl_parameters = "<<res<<endl;
  }
  //------------------------------------------------------------------------------------
  // Create the interaction and set up the interaction if the parameters
  // exist, no interaction created otherwise
  if(res){
    Hamiltonian_MM inter;
    if(syst.is_Box){ inter.set_pbc(&syst.Box,0,0,0); }
    if(int_type=="mb"||int_type=="mb_excl"){
    
      VECTOR** r; r = new VECTOR*[sz];
      VECTOR** f; f = new VECTOR*[sz];
      VECTOR** g; g = new VECTOR*[sz];
      VECTOR** m; m = new VECTOR*[sz]; 
      double** q; q = new double*[sz];
      double** epsilon; epsilon = new double*[sz];
      double** sigma;  sigma = new double*[sz];
      double** displr_2; displr_2 = new double*[sz];
//      double*  displT_2; displT_2 = &dT_2;
      vector<std::string> types;

      for(int i=0; i<sz; i++){
        types.push_back(syst.Atoms[at[i]].Atom_ff_type);
        r[i] = &syst.Atoms[at[i]].Atom_RB.rb_cm;
        int g_indx = syst.Atoms[at[i]].globGroup_Index;
        int m_indx = syst.Atoms[at[i]].globMolecule_Index;
        g[i] = &syst.Fragments[g_indx].Group_RB.rb_cm;
        m[i] = &syst.Molecules[m_indx].Molecule_RB.rb_cm;
        f[i] = &syst.Atoms[at[i]].Atom_RB.rb_force;
        q[i] = &syst.Atoms[at[i]].Atom_charge;
        displr_2[i] = &syst.Atoms[at[i]].Atom_displ2;
        epsilon[i] = new double; *epsilon[i] = 0.0;
        sigma[i] = new double;  *sigma[i] = 3.0;
      }
      res = ff.set_ff_charges(sz, types, r, q);
      res = ff.get_vdw_parameters(sz,types,epsilon,sigma);
//      res = ff.set_ff_epsilon_and_sigma(sz, types, epsilon, sigma); // vdw parameters


      vector<int> vexcl1, vexcl2;
      vector<double> vscale;
      vector< vector<excl_scale> > excl_scales; // combined structure

      double scale12,scale13,scale14;
      scale12 = 0.0; scale13 = 0.0; scale14 = 1.0; // default values
      if(int_type=="mb"){
        if(ff.mb_functional=="Ewald_3D"){ scale12 = ff.elec_scale12; scale13 = ff.elec_scale13; scale14 = ff.elec_scale14; }
        else if(ff.mb_functional=="vdw_LJ"||ff.mb_functional=="vdw_LJ1"){ scale12 = ff.vdw_scale12; scale13 = ff.vdw_scale13; scale14 = ff.vdw_scale14; }

        else if(ff.mb_functional=="LJ_Coulomb"){ 
          //!!! For now assume only vdw-based scaling factors
          scale12 = ff.vdw_scale12; scale13 = ff.vdw_scale13; scale14 = ff.vdw_scale14;
        }
      }
      else if(int_type=="mb_excl"){
        if(ff.mb_excl_functional=="vdw_LJ1"){ scale12 = ff.vdw_scale12; scale13 = ff.vdw_scale13; scale14 = ff.vdw_scale14; }
      }

      if(verb>1){  cout<<"=============== Exclusion map ========================\n"; }

      int num_excl = 0;
      for(i=0;i<sz;i++){
        vector<int> ve1,ve2,vs;
        for(int j=i;j<sz;j++){ 
          if(syst.is_group_pair(at[i],at[j])){ vexcl1.push_back(at[i]); vexcl2.push_back(at[j]); vscale.push_back(0.0); 
                                               ve1.push_back(at[i]); ve2.push_back(at[j]); vs.push_back(0.0);
                                        }
          else{
            if(syst.is_12pair(at[i],at[j])) { vexcl1.push_back(at[i]); vexcl2.push_back(at[j]); vscale.push_back(scale12); 
                                              ve1.push_back(at[i]); ve2.push_back(at[j]); vs.push_back(scale12);
                                       }
            else{
              if(syst.is_13pair(at[i],at[j])) { vexcl1.push_back(at[i]); vexcl2.push_back(at[j]); vscale.push_back(scale13); 
                                                ve1.push_back(at[i]); ve2.push_back(at[j]); vs.push_back(scale13);
                                         }
              else{
                if(syst.is_14pair(at[i],at[j])) { vexcl1.push_back(at[i]); vexcl2.push_back(at[j]); vscale.push_back(scale14);
                                                  ve1.push_back(at[i]); ve2.push_back(at[j]); vs.push_back(scale14);
                                           }         
              }
            }
          }
        }//for j

        // Add new entry in excl_scales:
        vector<excl_scale> excli;
        for(j=0;j<ve1.size();j++){
          excl_scale exclij;          
          exclij.at_indx1 = ve1[j];
          exclij.at_indx2 = ve2[j];
          exclij.scale = vs[j];
      //    if(int_type=="mb"){ excli.push_back(exclij); num_excl++; }
      //    else if(int_type=="mb_excl"){  if(vs[j]>0.0) { excli.push_back(exclij); num_excl++; } }          
          excli.push_back(exclij); num_excl++;
        }
        if(!excli.empty()){ excl_scales.push_back(excli);  }
//        cout<<"i = "<<i<<" size = "<<excli.size()<<endl;
      }// for i

//      cout<<"=============== Exclusion map ========================\n";
//      for()
//      cout<<"======================================================\n";

      int nexcl = vscale.size();
      int* excl1; excl1 = new int[nexcl];
      int* excl2; excl2 = new int[nexcl];
      double* scale; scale = new double[nexcl];
      cout<<"###############################\n";
      if(verb>0){ cout<<" nexcl = "<<nexcl<<endl; }
      
      for(i=0;i<nexcl;i++){
        excl1[i] = vexcl1[i];
        excl2[i] = vexcl2[i];
        scale[i] = vscale[i];
        if(verb>1){ cout<<"i = "<<i<<" excl1 = "<<excl1[i]<<" excl2 = "<<excl2[i]<<" scale = "<<scale[i]<<endl;}
      }

      inter.set_mb_interaction(int_type,ff.mb_functional,sz, id, r, g, m, f, q, epsilon, sigma, nexcl, excl1, excl2, scale, displr_2, &syst.dT_2, excl_scales, prms);
      if(verb>0){
        cout<<int_type<<" interaction created\n";
        cout<<"Number of exclusions = "<<num_excl<<endl;
      }

    }// if mb

//    exit(0);

    // Before adding - should probably check for existence of the interaction in the interactions array
    if(is_new_interaction(inter)){
      if(inter.get_status()) { active_interactions.push_back(interactions.size()); }
//      cout<<"adding interaction of memory size = "<<sizeof(inter)<<endl;
      interactions.push_back(inter);
    }// if interaction is new

  }// if res
  }// if x == true

}


void listHamiltonian_MM::set_group_interactions_for_atoms
(System& syst, string int_type,vector<Group>& top_elt,vector<int>& lst1,vector<int>& lst2,ForceField& ff){
/**************************************************************
  This function determines the parameters for group (=bonds,angles,dihedrals,etc.) interactions
  corresponding to the force field and the functional form of the
  potential in use. It then creates the interaction objects
 (interactions) for all groups defined in the system. It deactivates
  the intra-fragmental group interactions.
  If the interactions were already defined - they are not included.
  Also if not all necessary parameters obtained - the interaction
  is not created.
***************************************************************/
  int g_indx[4];
  int m_indx[4];
  for(int b=0;b<top_elt.size();b++){
    int sz = top_elt[0].Group_Size;
    bool x = true;
    int at[4];
    for(int j=0;j<sz;j++){
      at[j] = top_elt[b].globAtom_Index[j];
      g_indx[j] = syst.Atoms[at[j]].globGroup_Index;
      m_indx[j] = syst.Atoms[at[j]].globMolecule_Index;
      x = (x && (is_in_vector(at[j],lst1)||is_in_vector(at[j],lst2)));
    }
    // If now x is true - then each of the atoms of the topological element
    // belongs to one of the lists
    if(x){

      //--------------- Get the interaction parameters from the force field --------------
      map<string,double> prms;      
      int res = 0;
      if(int_type=="bond"){
        double bond_order = 1.0;
        if(top_elt[b].is_Group_bond_order){ bond_order = top_elt[b].Group_bond_order;   }
        res = ff.get_bond_parameters(syst.Atoms[at[0]].Atom_ff_type, syst.Atoms[at[1]].Atom_ff_type,bond_order,prms); 
      } 
      else if(int_type=="angle"){
        double bo12,bo23;
        bo12 = bo23 = 1.0;
        int a,b;
        int coord = 2;
        if(syst.Atoms[at[1]].is_Atom_coordination){ coord = syst.Atoms[at[1]].Atom_coordination; }
        b = syst.Find_Bond(at[0],at[1]); if(b>-1){ if(syst.Bonds[b].is_Group_bond_order) { bo12 = syst.Bonds[b].Group_bond_order; } }
        b = syst.Find_Bond(at[1],at[2]); if(b>-1){ if(syst.Bonds[b].is_Group_bond_order) { bo23 = syst.Bonds[b].Group_bond_order; } }

        res = ff.get_angle_parameters(syst.Atoms[at[0]].Atom_ff_type, syst.Atoms[at[1]].Atom_ff_type, syst.Atoms[at[2]].Atom_ff_type,
                                      bo12,bo23,coord,prms);
      }
      else if(int_type=="dihedral"){
        double bo12,bo23,bo34;
        bo12 = bo23 = bo34 = 1.0;
        int b;
        b = syst.Find_Bond(at[0],at[1]); if(b>-1){ if(syst.Bonds[b].is_Group_bond_order) { bo12 = syst.Bonds[b].Group_bond_order; } }
        b = syst.Find_Bond(at[1],at[2]); if(b>-1){ if(syst.Bonds[b].is_Group_bond_order) { bo23 = syst.Bonds[b].Group_bond_order; } }
        b = syst.Find_Bond(at[2],at[3]); if(b>-1){ if(syst.Bonds[b].is_Group_bond_order) { bo34 = syst.Bonds[b].Group_bond_order; } }

        res = ff.get_dihedral_parameters(syst.Atoms[at[0]].Atom_ff_type, syst.Atoms[at[1]].Atom_ff_type,
                                         syst.Atoms[at[2]].Atom_ff_type, syst.Atoms[at[3]].Atom_ff_type,
                                         bo12,bo23,bo34,prms);
      }
      else if(int_type=="oop"){
        res = ff.get_oop_parameters(syst.Atoms[at[0]].Atom_ff_type, syst.Atoms[at[1]].Atom_ff_type,
                                    syst.Atoms[at[2]].Atom_ff_type, syst.Atoms[at[3]].Atom_ff_type, prms);
      }
      else if(int_type=="vdw"){
        // The order matters here (for small cycles it is more important to exclude 1,2-pairs
        // than 1,3-pairs and it is more important to exclude 1,3-pairs that 1,4-pairs)
        std::string excl_pairs = "no";
        if(syst.is_14pair(at[0],at[1])){ excl_pairs = "14";}
        if(syst.is_13pair(at[0],at[1])){ excl_pairs = "13";}
        if(syst.is_12pair(at[0],at[1])){ excl_pairs = "12";}

        res = ff.get_vdw_parameters(syst.Atoms[at[0]].Atom_ff_type, syst.Atoms[at[1]].Atom_ff_type, excl_pairs, prms);
      }
      else if(int_type=="elec"){
        // The order matters here (for small cycles it is more important to exclude 1,2-pairs
        // than 1,3-pairs and it is more important to exclude 1,3-pairs that 1,4-pairs)
        std::string excl_pairs = "no";
        if(syst.is_14pair(at[0],at[1])){ excl_pairs = "14";}
        if(syst.is_13pair(at[0],at[1])){ excl_pairs = "13";}
        if(syst.is_12pair(at[0],at[1])){ excl_pairs = "12";}
        double q1,q2; q1 = q2 = 0.0;
        if(syst.Atoms[at[0]].is_Atom_charge){ q1 = syst.Atoms[at[0]].Atom_charge; }
        if(syst.Atoms[at[1]].is_Atom_charge){ q2 = syst.Atoms[at[1]].Atom_charge; }

        res = ff.get_elec_parameters(syst.Atoms[at[0]].Atom_ff_type, syst.Atoms[at[1]].Atom_ff_type,
                        excl_pairs,q1,q2, syst.Atoms[at[0]].is_Atom_charge, syst.Atoms[at[1]].is_Atom_charge,prms);
//        cout<<"Atoms[at[0]].Atom_charge = "<<Atoms[at[0]].Atom_charge<<endl;
//        cout<<"Atoms[at[1]].Atom_charge = "<<Atoms[at[1]].Atom_charge<<endl;
//        cout<<"Parameters for elec exist? "<<res<<endl;
      }

      //------------------------------------------------------------------------------------
      // Create the interaction and set up the interaction if the parameters
      // exist, no interaction created otherwise
      if(res){
        Hamiltonian_MM inter;
        if(syst.is_Box){ inter.set_pbc(&syst.Box,0,0,0); }
        if(int_type=="bond"){
          inter.set_2a_interaction("bond",ff.bond_functional,
                                   syst.Atoms[at[0]].Atom_id, syst.Atoms[at[1]].Atom_id,
                                   syst.Atoms[at[0]].Atom_RB.rb_cm, syst.Atoms[at[1]].Atom_RB.rb_cm,
                                   syst.Fragments[g_indx[0]].Group_RB.rb_cm, syst.Fragments[g_indx[1]].Group_RB.rb_cm,
                                   syst.Molecules[m_indx[0]].Molecule_RB.rb_cm, syst.Molecules[m_indx[1]].Molecule_RB.rb_cm,
                                   syst.Atoms[at[0]].Atom_RB.rb_force, syst.Atoms[at[1]].Atom_RB.rb_force,
                                   syst.Atoms[at[0]].Atom_displ2, syst.Atoms[at[1]].Atom_displ2, syst.dT_2, prms);
        }
        else if(int_type=="angle"){
          inter.set_3a_interaction("angle",ff.angle_functional,
                                   syst.Atoms[at[0]].Atom_id, syst.Atoms[at[1]].Atom_id, syst.Atoms[at[2]].Atom_id,
                                   syst.Atoms[at[0]].Atom_RB.rb_cm, syst.Atoms[at[1]].Atom_RB.rb_cm, syst.Atoms[at[2]].Atom_RB.rb_cm,
                                   syst.Fragments[g_indx[0]].Group_RB.rb_cm, syst.Fragments[g_indx[1]].Group_RB.rb_cm, syst.Fragments[g_indx[2]].Group_RB.rb_cm,
                                   syst.Molecules[m_indx[0]].Molecule_RB.rb_cm, syst.Molecules[m_indx[1]].Molecule_RB.rb_cm, syst.Molecules[m_indx[2]].Molecule_RB.rb_cm,
                                   syst.Atoms[at[0]].Atom_RB.rb_force, syst.Atoms[at[1]].Atom_RB.rb_force, syst.Atoms[at[2]].Atom_RB.rb_force, prms);
        }
        else if(int_type=="dihedral"){
          inter.set_4a_interaction("dihedral",ff.dihedral_functional,
          syst.Atoms[at[0]].Atom_id, syst.Atoms[at[1]].Atom_id, syst.Atoms[at[2]].Atom_id, syst.Atoms[at[3]].Atom_id,
          syst.Atoms[at[0]].Atom_RB.rb_cm, syst.Atoms[at[1]].Atom_RB.rb_cm, syst.Atoms[at[2]].Atom_RB.rb_cm, syst.Atoms[at[3]].Atom_RB.rb_cm,
          syst.Fragments[g_indx[0]].Group_RB.rb_cm, syst.Fragments[g_indx[1]].Group_RB.rb_cm, syst.Fragments[g_indx[2]].Group_RB.rb_cm, syst.Fragments[g_indx[3]].Group_RB.rb_cm,
          syst.Molecules[m_indx[0]].Molecule_RB.rb_cm, syst.Molecules[m_indx[1]].Molecule_RB.rb_cm, syst.Molecules[m_indx[2]].Molecule_RB.rb_cm, syst.Molecules[m_indx[3]].Molecule_RB.rb_cm,
          syst.Atoms[at[0]].Atom_RB.rb_force, syst.Atoms[at[1]].Atom_RB.rb_force, syst.Atoms[at[2]].Atom_RB.rb_force, syst.Atoms[at[3]].Atom_RB.rb_force, prms);
        }
        else if(int_type=="oop"){
          inter.set_4a_interaction("oop",ff.oop_functional,
          syst.Atoms[at[0]].Atom_id, syst.Atoms[at[1]].Atom_id, syst.Atoms[at[2]].Atom_id, syst.Atoms[at[3]].Atom_id,
          syst.Atoms[at[0]].Atom_RB.rb_cm, syst.Atoms[at[1]].Atom_RB.rb_cm, syst.Atoms[at[2]].Atom_RB.rb_cm, syst.Atoms[at[3]].Atom_RB.rb_cm,
          syst.Fragments[g_indx[0]].Group_RB.rb_cm, syst.Fragments[g_indx[1]].Group_RB.rb_cm, syst.Fragments[g_indx[2]].Group_RB.rb_cm, syst.Fragments[g_indx[3]].Group_RB.rb_cm,
          syst.Molecules[m_indx[0]].Molecule_RB.rb_cm, syst.Molecules[m_indx[1]].Molecule_RB.rb_cm, syst.Molecules[m_indx[2]].Molecule_RB.rb_cm, syst.Molecules[m_indx[3]].Molecule_RB.rb_cm,
          syst.Atoms[at[0]].Atom_RB.rb_force, syst.Atoms[at[1]].Atom_RB.rb_force, syst.Atoms[at[2]].Atom_RB.rb_force, syst.Atoms[at[3]].Atom_RB.rb_force, prms);

        }
        else if(int_type=="vdw"){
          inter.set_2a_interaction("vdw",ff.vdw_functional,
                                   syst.Atoms[at[0]].Atom_id, syst.Atoms[at[1]].Atom_id,
                                   syst.Atoms[at[0]].Atom_RB.rb_cm, syst.Atoms[at[1]].Atom_RB.rb_cm,
                                   syst.Fragments[g_indx[0]].Group_RB.rb_cm, syst.Fragments[g_indx[1]].Group_RB.rb_cm,
                                   syst.Molecules[m_indx[0]].Molecule_RB.rb_cm, syst.Molecules[m_indx[1]].Molecule_RB.rb_cm,
                                   syst.Atoms[at[0]].Atom_RB.rb_force, syst.Atoms[at[1]].Atom_RB.rb_force,
                                   syst.Atoms[at[0]].Atom_displ2, syst.Atoms[at[1]].Atom_displ2, syst.dT_2, prms);
        }
        else if(int_type=="elec"){
          inter.set_2a_interaction("elec",ff.elec_functional,
                                   syst.Atoms[at[0]].Atom_id, syst.Atoms[at[1]].Atom_id,
                                   syst.Atoms[at[0]].Atom_RB.rb_cm, syst.Atoms[at[1]].Atom_RB.rb_cm,
                                   syst.Fragments[g_indx[0]].Group_RB.rb_cm, syst.Fragments[g_indx[1]].Group_RB.rb_cm,
                                   syst.Molecules[m_indx[0]].Molecule_RB.rb_cm, syst.Molecules[m_indx[1]].Molecule_RB.rb_cm,
                                   syst.Atoms[at[0]].Atom_RB.rb_force, syst.Atoms[at[1]].Atom_RB.rb_force,
                                   syst.Atoms[at[0]].Atom_displ2, syst.Atoms[at[1]].Atom_displ2, syst.dT_2, prms);
        }

        // Before adding - should probably check for existence of the interaction in the interactions array
        if(is_new_interaction(inter)){
          // If the interaction is inactive (e.g. bond inside of the rigid group) - deactivate it, but
          // keep for other possible use
          if(sz==2)      { if(!is_active(syst.Atoms[at[0]],syst.Atoms[at[1]])){  inter.deactivate();   } }
          else if(sz==3) { if(!is_active(syst.Atoms[at[0]],syst.Atoms[at[1]],syst.Atoms[at[2]])){  inter.deactivate();   } }
          else if(sz==4) { if(!is_active(syst.Atoms[at[0]],syst.Atoms[at[1]],syst.Atoms[at[2]],syst.Atoms[at[3]])){  inter.deactivate();   } }

          if(inter.get_status()) { 
/*
            if(int_type=="bond"){
            cout<<"Setting bond interaction \n";
            cout<<"indx1 = "<<at[0]<<" indx2 = "<<at[1]<<endl;
            cout<<"id1 = "<<Atoms[at[0]].Atom_id<<" id2 = "<<Atoms[at[1]].Atom_id<<endl;
            cout<<"r1 = "<<Atoms[at[0]].Atom_RB.rb_cm<<" r2 = "<<Atoms[at[1]].Atom_RB.rb_cm<<endl;
            inter.show_info();
            }

            if(int_type=="angle"){ 
            cout<<"Setting angle interaction \n";
            cout<<"indx1 = "<<at[0]<<" indx2 = "<<at[1]<<" indx3 = "<<at[2]<<endl;
            cout<<"id1 = "<<Atoms[at[0]].Atom_id<<" id2 = "<<Atoms[at[1]].Atom_id<<" id3 = "<<Atoms[at[2]].Atom_id<<endl;
            cout<<"r1 = "<<Atoms[at[0]].Atom_RB.rb_cm<<" r2 = "<<Atoms[at[1]].Atom_RB.rb_cm<<" r3 = "<<Atoms[at[2]].Atom_RB.rb_cm<<endl;
            inter.show_info();
            }
*/
            active_interactions.push_back(interactions.size()); 
          }
          interactions.push_back(inter);
        }// if interaction is new
      }// if parameters exist
    }// if the atoms of the bond belong to the lst1 and lst2
  }// for all bonds in the system
}




void listHamiltonian_MM::set_interactions_for_atoms
(System& syst, boost::python::list lst1,boost::python::list lst2,ForceField& ff,int verb, int assign_rings){
/**************************************************************
  This function sets up all interactions existing in system, 
  that is it determines the parameters for given force_fields
  of potential and relates them with the functional form and
  the atoms and groups of the system
  <verb> - defines verbosity level
***************************************************************/
  if(verb>0){ cout<<"Setting interactions with the force field "<<ff.ForceField_Name<<endl; }

  // Transform the lists of id into vectors of indexes
  vector<int> vlst1,vlst2;
  int i,id,indx,sz;

  sz = len(lst1);
  for(i=0;i<sz;i++){
    id = extract<int>(lst1[i]);
    indx = syst.get_atom_index_by_atom_id(id);
    vlst1.push_back(indx);
  }

  sz = len(lst2);
  for(i=0;i<sz;i++){
    id = extract<int>(lst2[i]);
    indx = syst.get_atom_index_by_atom_id(id);
    vlst2.push_back(indx);
  }

  // Calculate functional groups assignments
  syst.determine_functional_groups(assign_rings);

  set_atom_types(syst,vlst1,ff); 
  set_atom_types(syst,vlst2,ff);

  // Setup the interactions for all atoms - many-body interactions 
  if(ff.is_mb_functional){ set_atom_interactions_for_atoms(syst,"mb",syst.Atoms,vlst1,vlst2,ff,verb);}            // Atoms
  if(ff.is_mb_excl_functional){ set_atom_interactions_for_atoms(syst,"mb_excl",syst.Atoms,vlst1,vlst2,ff,verb);}  // Atoms - exclusions

  // Setup the interactions for differen topological groups
  if(ff.is_bond_functional){ set_group_interactions_for_atoms(syst,"bond",syst.Bonds,vlst1,vlst2,ff); }      // Bonds  
  if(ff.is_angle_functional){ set_group_interactions_for_atoms(syst,"angle",syst.Angles,vlst1,vlst2,ff); }   // Angles
  if(ff.is_dihedral_functional){ set_group_interactions_for_atoms(syst,"dihedral",syst.Dihedrals,vlst1,vlst2,ff); } // Dihedrals/Torsions
  if(ff.is_oop_functional){ set_group_interactions_for_atoms(syst,"oop",syst.Impropers,vlst1,vlst2,ff); }     // OOP
  if(ff.is_vdw_functional){ set_group_interactions_for_atoms(syst,"vdw",syst.Pairs,vlst1,vlst2,ff); }         // Vdw
  if(ff.is_elec_functional){ set_group_interactions_for_atoms(syst,"elec",syst.Pairs,vlst1,vlst2,ff); }       // Electrostatic

  // Actual pbc dimensions should be determined using R_cut parameters
//  apply_pbc_to_interactions("vdw",1,1,1);
//  apply_pbc_to_interactions("elec",1,1,1);

 //cout<<"Interactions are set: Size of the interactions = "<<sizeof(interactions)<<endl;

}


void listHamiltonian_MM::set_interactions_for_fragments
(System& syst, boost::python::list lst1,boost::python::list lst2,ForceField& ff){
/**************************************************************
  This function sets up all interactions existing in system,
  that is it determines the parameters for given force_fields
  of potential and relates them with the functional form and
  the atoms and groups of the system
***************************************************************/

  // Transform the lists of id into vectors of indexes
  vector<int> vlst1,vlst2;
  int i,id,indx,sz;

  sz = len(lst1);
  for(i=0;i<sz;i++){
    id = extract<int>(lst1[i]);
    indx = syst.get_fragment_index_by_fragment_id(id);
    vlst1.push_back(indx);
  }

  sz = len(lst2);
  for(i=0;i<sz;i++){
    id = extract<int>(lst2[i]);
    indx = syst.get_fragment_index_by_fragment_id(id);
    vlst2.push_back(indx);
  }

  // Calculate functional groups assignments
//  determine_functional_groups();

  set_fragment_types(syst,vlst1,ff);
  set_fragment_types(syst,vlst2,ff);

  // Setup the interactions for all atoms - many-body interactions
//  if(ff.is_mb_functional){ set_atom_interactions_for_atoms("mb",Atoms,vlst1,vlst2,ff);}            // Atoms


}



void listHamiltonian_MM::apply_pbc_to_interactions(System& syst, int int_type,int nx,int ny,int nz){
/**********************************************************************
  Interactions will be multiplied by the number of possible distinct
  translations of maximal degrees nx,ny and nz, that is:
  1 -> (2*nx+1)*(2*ny+1)*(2*nz+1) new interactions
  one interaction for each distinct combination of kx,ky,kz:
  nx<=kx<=nx; ny<=ky<=ny; nz<=kz<=nz;
***********************************************************************/
  if(syst.is_Box){

  int sz = interactions.size();
  for(int i=0;i<sz;i++){
    if(interactions[i].get_type()==int_type){
      if(interactions[i].is_origin()){
        for(int kx=-nx;kx<=nx;kx++){
          for(int ky=-ny;ky<=ny;ky++){
            for(int kz=-nz;kz<=nz;kz++){
              Hamiltonian_MM inter(interactions[i]);
              inter.set_pbc(&syst.Box,kx,ky,kz);
              if(is_new_interaction(inter)){ interactions.push_back(inter); }
          }// kz
        }// ky
      }// kx
    }// if is_origin
    } // if ==int_type
  }// for i
  }// Box!= NULL

}




void listHamiltonian_MM::set_respa_types(std::string s_int_type,std::string s_respa_type){
// Must be called after one or more new interactions are created

  int sz = interactions.size();
  // The following conversion should be consistent with function:
  // void Interaction::set_interaction_type_and_functional(std::string t,std::string f)
  int int_type = -1;
  if(s_int_type=="bond")         { int_type = 0; }
  else if(s_int_type=="angle")   { int_type = 1; }
  else if(s_int_type=="dihedral"){ int_type = 2; }
  else if(s_int_type=="oop")     { int_type = 3; }
  else if(s_int_type=="vdw")     { int_type = 4; }
  else if(s_int_type=="elec")    { int_type = 5; }
  else if(s_int_type=="mb")      { int_type = 6; }
  else if(s_int_type=="cg")      { int_type = 7; }
  else if(s_int_type=="mb_excl") { int_type = 8; }
  else{ cout<<"Error in set_respa_types: s_int_type = "<<s_int_type<<" is unknown\n"; exit(0); }

  int respa_type = -1;
  if(s_respa_type=="fast") { respa_type = 0; }
  else if(s_respa_type=="medium"){  respa_type = 1; }
  else if(s_respa_type=="slow"){ respa_type = 2; }
  else{ cout<<"Error in set_respa_types: s_respa_type = "<<s_respa_type<<"is unknown\n"; exit(0); }

  for(int i=0;i<sz;i++){
    interactions[i].set_respa_type(int_type,respa_type);
  }

}






}// namespace libhamiltonian_mm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian

