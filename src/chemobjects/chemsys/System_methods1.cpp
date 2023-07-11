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
  \file System_methods1.cpp
  \brief The file implements some topological functions
    
*/

#include "System.h"

/// liblibra namespace
namespace liblibra{


/// libchemobjects namespace
namespace libchemobjects{

/// libchemsys namespace
namespace libchemsys{


//================ Not-member functions =================================

void System::Generate_Connectivity_Matrix(){
/** 
  Generate a group connectivity matrix - if atoms of a bond belong to different groups, then these two groups are connected
  We do not need atomic connectivity matrix since it is too big
*/

  GroupConnMatrix = new MATRIX(Number_of_fragments,Number_of_fragments); // Group connectivity matrix
  *GroupConnMatrix = 0.0;
  is_GroupConnMatrix =1;

  int b,g1,g2;
  for(int i=0;i<Number_of_frag_bonds;i++){
    b  = Frag_bonds[i];
    g1 = Bonds[b].globAtom_Index[0];
    g2 = Bonds[b].globAtom_Index[1];

    GroupConnMatrix->M[g1*Number_of_fragments+g2]=1.0;
    GroupConnMatrix->M[g2*Number_of_fragments+g1]=1.0;
  }
}


void System::Assign_Rings(){
/**
  This function assings the rings - meaning it determines how many and which rings are present in the system. 
  Carefully - this may become a very slow procedure for large and dense graphs. 
*/

  cout<<"In System::Assign_Rings\n";
  typedef GRAPH<Atom*,Group*> Graph;
  Graph g;
  // Create a vertices of the graph
  for(int i2=0;i2<Number_of_atoms;i2++){
    VERTEX<Atom*> vrtx(&Atoms[i2]);
    g.ADD_VERTEX(vrtx);
  }
  // Create edges of the graph
  for(int i3=0;i3<Number_of_bonds;i3++){
    int indx1 = Bonds[i3].globAtom_Index[0];
    int indx2 = Bonds[i3].globAtom_Index[1];
    EDGE<Group*> edg(indx1,indx2,0,&Bonds[i3]);
    g.ADD_EDGE(edg);
  }

  vector<Path> rings;

//  std::cout<<std::endl<<"Original graph"<<std::endl;
//  int _V_ = g.show_vertices(1);
//  int _E_ = g.show_edges(1);
//  std::cout<<"_V_ = "<<_V_<<std::endl;
//  std::cout<<"_E_ = "<<_E_<<std::endl;

  // Reduce graph 
  g.REDUCE_GRAPH();

//  std::cout<<std::endl<<"Reduced Original graph"<<std::endl;
//  _V_ = g.show_vertices(1);
//  _E_ = g.show_edges(1);
//  std::cout<<"_V_ = "<<_V_<<std::endl;
//  std::cout<<"_E_ = "<<_E_<<std::endl;

  // Find the smallest rings 
  // !!! This option is commented by default - uncomment it when you need the information about rings !!!!
  g.FIND_SSSR(rings);

  // Now assign
  int indx,i,j,nrings,sz;
  nrings = rings.size();
//  cout<<"nrings = "<<nrings<<endl;
  for(i=0;i<nrings;i++){
    sz = rings[i].size();
    Group rng;
    int loc_indx = 0;
//    cout<<"Ring "<<i<<" : ";

    for(j=0;j<sz;j++){
      int indx1 = g.E[rings[i][j]].edge_vertex_id1;
      int indx2 = g.E[rings[i][j]].edge_vertex_id2;

      if(!is_in_vector(indx1,rng.globAtom_Index)){
        rng.globAtom_Index.push_back(indx1);
        rng.locAtom_Index.push_back(loc_indx++); 
        rng.Group_Size++;
        rng.globGroup_Size++;
        rng.locGroup_Size++;
      }
      if(!is_in_vector(indx2,rng.globAtom_Index)){ 
        rng.globAtom_Index.push_back(indx2);
        rng.locAtom_Index.push_back(loc_indx++);
        rng.Group_Size++;
        rng.globGroup_Size++;
        rng.locGroup_Size++;
      }

//      cout<<Atoms[indx1].Atom_id<<"  "<<Atoms[indx2].Atom_id<<"  ";

      if(!is_in_vector(sz,Atoms[indx1].Atom_ring_sizes)){
        Atoms[indx1].Atom_ring_sizes.push_back(sz);
        Atoms[indx1].is_Atom_ring_sizes = 1;
      }

      if(!is_in_vector(sz,Atoms[indx2].Atom_ring_sizes)){
        Atoms[indx2].Atom_ring_sizes.push_back(sz);
        Atoms[indx2].is_Atom_ring_sizes = 1;
      }
    }// for j
//    cout<<endl;
    rng.globGroup_Index = Number_of_rings;
    Rings.push_back(rng);
    Number_of_rings++;

  }// for i

  for(i=0;i<Number_of_atoms;i++){
    sz = Atoms[i].Atom_ring_sizes.size();
    cout<<"i= "<<i<<" number of rings this atom belong to = "<<sz<<endl;
    int min_sz=1;
    if(sz>=1){ min_sz = Atoms[i].Atom_ring_sizes[0];}
    else{
      for(j=0;j<sz-1;j++){
        min_sz=(Atoms[i].Atom_ring_sizes[j]<=Atoms[i].Atom_ring_sizes[j+1])?Atoms[i].Atom_ring_sizes[j]:Atoms[i].Atom_ring_sizes[j+1];
      }
    }      
    Atoms[i].Atom_min_ring_size = min_sz;

    cout<<"Atom i="<<i<<" belongs to the minimal ring of size = "<<Atoms[i].Atom_min_ring_size<<endl;
  }// for i


}

void System::DIVIDE_GRAPH(int IndxStart,int IndxException, vector<int>& P){
/**
  \param[in] IndxStart The index of the vertex that, together with all connected (including indirectly) to it vertices will constitute one of the subgraphs
  \param[in] IndxException The index of the vertex connected to IndxStart, which with all connected (including indirectly) to it vertices will constitute second subgraph
  \param[out] P The list of vertex indices of the first subgraph

  This function returns to the vector P indexes corresponding to the part of the graph
  given by connectivity matrix ConnMatrix. This part starts at vertex IndxStart in
  the direction opposite to IndxException
  g1 and g2 - are indexes of starting and finishing group

  This function is needed to split the connectivity graph into two parts, divided by a single bond IndxStart - IndxException
  eventually, this is to be used for rotating one part of the system w.r.t. another around the bond selected by the above
  pair of indices
*/

  int g1 = IndxStart;
  int g2 = IndxException;
  P.push_back(g1);
  for(int i=0;i<GroupConnMatrix->n_cols;i++){
    if((GroupConnMatrix->M[g1*GroupConnMatrix->n_cols+i]==1.0)&&(i!=g2)){
      int sz = P.size();
      int st = 1;
      for(int j=0;j<sz;j++){ if(P[j]==i){ st=0;}   }
      if(st==1){  DIVIDE_GRAPH(i,IndxException,P); }
    }// if
  }// for i
}


}// namespace libchemsys
}// namespace libchemobjects
}// liblibra
