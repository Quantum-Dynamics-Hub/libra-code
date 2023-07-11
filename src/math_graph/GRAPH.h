/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef GRAPH_H
#define GRAPH_H

#if defined(USING_PCH)
#include "../pch.h"
#else

#include <vector>
#include <iostream>
#include <iomanip>
#include <map>
#include <math.h>

#endif


/// liblibra namespace
namespace liblibra{

using namespace std;

/// libgraph namespace
namespace libgraph{

typedef vector<int> Path;


// ---------- Helper functions ------------

template <class DATA_TYPE> int delete_item(vector<DATA_TYPE>& data, DATA_TYPE item){

   int sz = data.size();
   for(int i=0;i<sz;i++){    if(data[i]==item) { data[i]=data[sz-1]; break;}    }
   data.pop_back();

   return 0;
}

template <class DATA_TYPE> int replace_item(vector<DATA_TYPE>& data, DATA_TYPE item_old, DATA_TYPE item_new){

   int sz = data.size();
   for(int i=0;i<sz;i++){    if(data[i]==item_old) { data[i]=item_new; break;}    }

   return 0;
}

int merge_paths(Path& result, Path p1, Path p2);
int is_included(Path& p,Path& P);
int is_included(Path& p,vector<Path>& P);
int path_xor(Path& p1, Path& p2, Path& res);
void show_path(Path& p);
void show_paths(vector<Path>& p);


//*********************************************************************

template <class VERTEX_DATA> class VERTEX;
template <class EDGE_DATA> class EDGE;
template <class VERTEX_DATA,class EDGE_DATA> class GRAPH;

//**********************************************************************

template <class VERTEX_DATA> class VERTEX{          
      public:
      // Data
      VERTEX_DATA vertex_data;       int is_vertex_data;

      // Properties         
      int vertex_index;               // Index of this vertex: corresponds to time when this vertex was created
      int vertex_id;                  // ID of this vertex: corresponds to vertex uniqueness (label property)
      int vertex_degree;              // Number of adjacent vertices
      std::string vertex_name;        // Name of this vertex
             
      vector<int> adjacent_vertices;  // Indices of this vertices.
      vector<int> adjacent_edges;     // Indices of edges which are directly 
                                      // connected to this vertex


      // Methods
      VERTEX() {
            vertex_degree=0;
            vertex_index = -1; 
            vertex_id = -1;
            vertex_name = "";
            is_vertex_data = 0;
      }
      VERTEX(int id) {
            vertex_degree=0;
            vertex_index = -1;
            vertex_id = id;
            vertex_name = "";
            is_vertex_data = 0;
      }
      VERTEX(VERTEX_DATA vd){ 
            vertex_data = vd; 
            vertex_degree=0;
            vertex_index = -1; 
            vertex_id = -1;
            vertex_name = ""; 
            is_vertex_data = 1;
      }
      
      VERTEX(const VERTEX&);         // Copy constructor

     ~VERTEX(){;;}  // Destructor

      int show_info(int);


};
template <class EDGE_DATA> class EDGE{
      public:
    
      // Data
      EDGE_DATA edge_data;          int is_edge_data;

      // Properties
      int edge_index;                // Index of this edge
      int edge_vertex1;              // Index of the first vertex
      int edge_vertex2;              // Index of the second vertex
      int edge_vertex_id1;           // ID of the first vertex
      int edge_vertex_id2;           // ID of the second vertex

      int edge_direction;            // -1  from 2 to 1;
                                     //  0  nondirected edge;
                                     //  1  from 1 to 2;
      double edge_weight;            // Cost/weight parameter of the edge

      // Methods
      EDGE(){
           edge_vertex_id1 = -1;
           edge_vertex_id2 = -1; 
           edge_vertex_id1 = 0;
           edge_vertex_id2 = 1;
           edge_direction = 0;
           edge_index = 0;
           is_edge_data = 0;
           edge_weight = 1.0;
      }
      EDGE(int id1,int id2,int dir) {
           edge_vertex1 = -1;
           edge_vertex2 = -1;
           edge_vertex_id1 = id1;
           edge_vertex_id2 = id2;

           edge_direction=dir;
           edge_index = 0;
           is_edge_data = 0; 
           edge_weight = 1.0;          
      }
      EDGE(int id1,int id2,int dir,EDGE_DATA ed) {
           edge_vertex1 = -1;
           edge_vertex2 = -1;
           edge_vertex_id1 = id1;
           edge_vertex_id2 = id2;
           edge_direction = dir;
           edge_index = 0;
           edge_data = ed;
           is_edge_data = 1; 
           edge_weight = 1;
      }
      EDGE(const EDGE&);             // Copy constructor

     ~EDGE(){;;}  // Destructor

      int set(int id1,int id2,int dir){
           edge_vertex1 = -1;
           edge_vertex2 = -1;
           edge_vertex_id1 = id1;
           edge_vertex_id2 = id2;
           edge_direction=dir;
           edge_index = 0;
           is_edge_data = 0;
     
           return 0;
      }
      int set(int id1,int id2,int dir,EDGE_DATA ed) {
           edge_vertex1 = -1;
           edge_vertex2 = -1;
           edge_vertex_id1 = id1;
           edge_vertex_id2 = id2;
           edge_direction = dir;
           edge_index = 0;
           edge_data = ed;
           is_edge_data = 1;

           return 0;
      }

      int show_info(int);
      
	  
};
//===============================================================
// GRAPH class represents the connetivity structure on data.
// It contains actual node and edge objects, so this might be
// used for any kind of data. For memory saving we indent to 
// give the complex vertex and edge data by references

template <class VERTEX_DATA,class EDGE_DATA> class GRAPH{

      // Settings
      int is_allow_parallel_edges;

      // Helper functions
      int vertex_id_to_index(int id);
      int is_vertex_id_reserved(int id); 
      int is_adjacent_vertices(int,int);
      int is_edge_exist(EDGE<EDGE_DATA>);   

      public:
 
      // Properties
      std::string graph_name;
      std::vector< VERTEX<VERTEX_DATA> > V;     int _V_; // _V_ - is a size of V array;
      std::vector< EDGE<EDGE_DATA> > E;         int _E_; // _E_ - is a size of E array;


      // Metohds
      GRAPH();
      GRAPH(int);

      GRAPH(const GRAPH&);           // Copy constructor    

     ~GRAPH(){;;}  // Destructor


      int ADD_VERTEX();
      int ADD_VERTEX(int id);
      int ADD_VERTEX(VERTEX<VERTEX_DATA>);    
      int ADD_EDGE(EDGE<EDGE_DATA>);

      int DELETE_VERTEX(int);
      int DELETE_EDGE(int);
//      int DELETE_EDGE(int,int);
     
      int REDUCE_GRAPH();
      int FIND_PATHS_FLOYD_WARSHALL(double** D,vector<Path>** P);
      int FIND_SSSR(vector<Path>&);

      int show_vertices(int);
      int show_edges(int);


};


//**************************************************************************************************
//******************** Now the declaration of methods is going *************************************
//**************************************************************************************************

// ------------------- Constructors ------------------------

template <class VERTEX_DATA,class EDGE_DATA> GRAPH<VERTEX_DATA,EDGE_DATA>::GRAPH(){
        is_allow_parallel_edges = 0;
        _V_ = 0;       
        _E_ = 0;    
}

template <class VERTEX_DATA,class EDGE_DATA> GRAPH<VERTEX_DATA,EDGE_DATA>::GRAPH(int v){
        is_allow_parallel_edges = 0;
        _V_ = v;     
        _E_ = 0;    
        for(int i=0;i<_V_;i++){   VERTEX<VERTEX_DATA> vx;  vx.vertex_index = i;  V.push_back(vx); }
}

// ------------------ Copy constructors -------------------------


template <class VERTEX_DATA> VERTEX<VERTEX_DATA>::VERTEX(const VERTEX<VERTEX_DATA>& v){
      is_vertex_data = v.is_vertex_data;
      if(is_vertex_data){
      vertex_data = v.vertex_data;
      }

      vertex_index = v.vertex_index;
      vertex_id = v.vertex_id;
      vertex_degree = v.vertex_degree;
      vertex_name = v.vertex_name;

      adjacent_vertices = v.adjacent_vertices;
      adjacent_edges    = v.adjacent_edges;

}

template <class EDGE_DATA> EDGE<EDGE_DATA>::EDGE(const EDGE<EDGE_DATA>& e){
      is_edge_data = e.is_edge_data;
      if(is_edge_data){
      edge_data = e.edge_data;
      }

      edge_index      = e.edge_index;
      edge_vertex1    = e.edge_vertex1;
      edge_vertex2    = e.edge_vertex2;
      edge_vertex_id1 = e.edge_vertex_id1;
      edge_vertex_id2 = e.edge_vertex_id2;
      edge_direction  = e.edge_direction;
      edge_weight     = e.edge_weight;

}


template <class VERTEX_DATA,class EDGE_DATA> GRAPH<VERTEX_DATA,EDGE_DATA>::GRAPH(const GRAPH<VERTEX_DATA,EDGE_DATA>& g){

      int i;
      graph_name    = g.graph_name;

      _V_           = g._V_;   V = g.V;
      _E_           = g._E_;   E = g.E;
  

}

// -------------------- Member-functions of VERTEX class --------------------
template <class VERTEX_DATA> int VERTEX<VERTEX_DATA>::show_info(int lindr){

// lindr - level of indirection
// -1  = &x
//  0  =  x
//  1  = *x

     std::cout<<"vertex( "<<vertex_index<<" ) info:"<<std::endl;
     std::cout<<" vertex_id = "<<vertex_id<<std::endl;
     std::cout<<" vertex_name = "<<vertex_name<<std::endl;
     std::cout<<" vertex_degree = "<<vertex_degree<<std::endl;
     if(vertex_degree>0){
        std::cout<<" adjacent_vertices = ";
        for(int i=0;i<vertex_degree;i++){ std::cout<<adjacent_vertices[i]<<" "; }
        std::cout<<" adjacent_edges = ";
        for(int i=0;i<vertex_degree;i++){ std::cout<<adjacent_edges[i]<<" "; }

     }
     if(is_vertex_data){
             if(lindr==-1){       std::cout<<" &vertex_data = "<<&vertex_data<<std::endl; }
        else if(lindr== 0){       std::cout<<"  vertex_data = "<<vertex_data<<std::endl; }
        else if(lindr== 1){       std::cout<<" *vertex_data = "<<*vertex_data<<std::endl; }
     }

     std::cout<<std::endl;

     return 0;
}


//---------------------- Member-functions of EDGE class --------------------
template <class EDGE_DATA> int EDGE<EDGE_DATA>::show_info(int lindr){

// lindr - level of indirection
// -1  = &x
//  0  =  x
//  1  = *x
     std::cout<<"edge( "<<edge_index<<" ) info:"<<std::endl;
     std::cout<<" ( edge_vertex1 , edge_vertex2 ) = ( "<<edge_vertex1<<" , "<<edge_vertex2<<" )"<<std::endl;
     std::cout<<" ( edge_vertex_id1 , edge_vertex_id2 ) = ( "<<edge_vertex_id1<<" , "<<edge_vertex_id2<<" )"<<std::endl;
     std::cout<<" edge_direction = "<<edge_direction<<std::endl;
     std::cout<<" edge_weight = "<<edge_weight<<std::endl;
     if(is_edge_data){
             if(lindr==-1){       std::cout<<" &edge_data = "<<&edge_data<<std::endl; }
        else if(lindr== 0){       std::cout<<"  edge_data = "<<edge_data<<std::endl; }
        else if(lindr== 1){       std::cout<<" *edge_data = "<<*edge_data<<std::endl; }
     }

     std::cout<<std::endl;

    return 0;
}

//---------------------- Various operaions on graph ---------------------
template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::is_vertex_id_reserved(int id){
  
    int res = 0; // Not reserved (= currently not in use in this graph)
    for(int i=0;i<_V_;i++){
        if(V[i].vertex_id == id) { res = 1; }
    }

    return res;
}

template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::vertex_id_to_index(int id){

    // This function return the index of the vertex whose id is equal to argument
    // If such vertex can not be found -1 is returned
    int res = -1;
    for(int i=0;i<_V_;i++){
        if(V[i].vertex_id == id) { res = i; }
    }

    return res;
}

template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::is_adjacent_vertices(int v1,int v2){

    int res = 0; // Not adjacent to each other
    int i,st1,st2;
    st1 = st2 = 0;

    for(i=0;i<V[v1].vertex_degree;i++){
        if(V[v1].adjacent_vertices[i]==v2){ st1 = 1; }
    }

    for(i=0;i<V[v2].vertex_degree;i++){
        if(V[v2].adjacent_vertices[i]==v1){ st2 = 1; }
    }

    res = st1 * st2;

    return res;
}



template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::ADD_VERTEX(){       
        VERTEX<VERTEX_DATA> vx;  vx.vertex_index = _V_;  
        int id = _V_;
        while(is_vertex_id_reserved(id)){   id++;   }
        vx.vertex_id = id;
        V.push_back(vx);   _V_++;     
        return _V_;
}

template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::ADD_VERTEX(int id){
        VERTEX<VERTEX_DATA> vx;  vx.vertex_index = _V_;
        
        if(is_vertex_id_reserved(id)){
            std::cout<<"Warning!: Can not create vertex with vertex_id = "<<id; 
            std::cout<<" Such vertex already exists\n";  
        }else{
            vx.vertex_id = id;
            V.push_back(vx);   _V_++;
        }
        return _V_;
}


template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::ADD_VERTEX(VERTEX<VERTEX_DATA> vx){   
    vx.vertex_index = _V_; 
    if(vx.vertex_id == -1){ // Not defined
        int id = _V_;
        while(is_vertex_id_reserved(id)){   id++;   }
        vx.vertex_id = id;
    }
    V.push_back(vx);   _V_++;
    return _V_;
}

template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::show_vertices(int lindr){

    for(int i=0;i<_V_;i++){     V[i].show_info(lindr);    }
    return _V_;
}

template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::show_edges(int lindr){

    for(int i=0;i<_E_;i++){     E[i].show_info(lindr);    }
    return _E_;
}


template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::is_edge_exist(EDGE<EDGE_DATA> e){

  // The edge exstence check is based on vertex_id property of vertices not on their vertex_index property

  int res = 0; // Edge does not exist
  int st11,st12,st21,st22;
 
  if(is_allow_parallel_edges){
      res = 0;
  }else{

      // Find indices of the corresponding vertices
      int v1 = vertex_id_to_index(e.edge_vertex_id1);
      int v2 = vertex_id_to_index(e.edge_vertex_id2);
 
      if((v1==-1) || (v2==-1)){
          res = 0;
      }else{

      // Check vertices for being adjacent to each other
      if(is_adjacent_vertices(v1,v2)){
                     
          st11 = (V[v1].vertex_id == e.edge_vertex_id1);
          st22 = (V[v2].vertex_id == e.edge_vertex_id2);
          st12 = (V[v1].vertex_id == e.edge_vertex_id2);
          st21 = (V[v2].vertex_id == e.edge_vertex_id1);

          if(e.edge_direction==0){ // Nondirectional edge

              res = ( ((st11==1)&&(st22==1)) || ((st12==1)&&(st21==1)) );

          }// if
          else{         // Directed edge
 
              res = ( (st11==1)&&(st22==1) );
          }// else

      }// if adjacent_vertices
      }// else


  }// else

  return res;

}

template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::ADD_EDGE(EDGE<EDGE_DATA> eg){

    int v1,v2;
    int curr_num_nodes = _V_;
    EDGE<EDGE_DATA> e;   e = eg;    
    // Find indices of the corresponding vertices
    v1 = vertex_id_to_index(e.edge_vertex_id1);
    v2 = vertex_id_to_index(e.edge_vertex_id2);  

    if(v1==-1){ ADD_VERTEX(e.edge_vertex_id1); v1 = _V_-1;}
    if(v2==-1){ ADD_VERTEX(e.edge_vertex_id2); v2 = _V_-1;}


    if(is_edge_exist(e)==0){ // Edge does not exist
    // Modify vertices

        V[v1].adjacent_vertices.push_back(v2);
        V[v1].adjacent_edges.push_back(_E_);
        V[v1].vertex_degree++;        
 
        V[v2].adjacent_vertices.push_back(v1);
        V[v2].adjacent_edges.push_back(_E_);
        V[v2].vertex_degree++;

    // Modify edges
        e.edge_vertex1 = v1;
        e.edge_vertex2 = v2;
        e.edge_vertex_id1 = e.edge_vertex_id1;
        e.edge_vertex_id2 = e.edge_vertex_id2;
        e.edge_index   = _E_;
        E.push_back(e);
        _E_++;

    }else{                   // Edge exists

        if(is_allow_parallel_edges){ // Add one more parallel edge

           // Modify vertices
           V[v1].vertex_degree++;
           V[v2].vertex_degree++;
           V[v1].adjacent_edges.push_back(_E_);
           V[v2].adjacent_edges.push_back(_E_);

           // Modify edges
           e.edge_index = _E_;
           E.push_back(e);
           _E_++;
        }else{                       // Do nothing
        }
    }

    return _E_;
}


template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::DELETE_EDGE(int e){

    int v1,v2;
    if(e>=_E_){ std::cout<<"Error: Can not delete edge with index = "<<e<<" No such edge in a graph\n";   }
    else{
     
      // Modify vertices
      v1 = E[e].edge_vertex1; 
      v2 = E[e].edge_vertex2;
      delete_item<int>(V[v1].adjacent_edges,e);
      delete_item<int>(V[v2].adjacent_edges,e);
      delete_item<int>(V[v1].adjacent_vertices,v2);
      delete_item<int>(V[v2].adjacent_vertices,v1);
      V[v1].vertex_degree--;
      V[v2].vertex_degree--;

      // Modify edges
      if(e==(_E_-1)){ ;; }
      else{
         int s1,s2;
         s1 = E[_E_-1].edge_vertex1;
         s2 = E[_E_-1].edge_vertex2;
         E[e] = E[_E_-1]; 
         E[e].edge_index = e;
         replace_item<int>(V[s1].adjacent_edges,(_E_-1),e);
         replace_item<int>(V[s2].adjacent_edges,(_E_-1),e);
        
      }

      // Modify graph
      E.pop_back();
      _E_--;
    }// else
   
    return _E_;
}


template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::DELETE_VERTEX(int v){
    
    if(v<_V_){
      // Delete all edges iteratively (because of possible index mutations)
      while(V[v].vertex_degree>0){
           this->DELETE_EDGE(V[v].adjacent_edges[0]);
      }
      // Now delte the vertex (reindex the other vertices)
      if(v==(_V_-1)){ ;; }
      else{
        // Modify vertices
        V[v] = V[_V_-1];
        V[v].vertex_index = v;

        // Modify edges
        for(int i=0;i<V[v].vertex_degree;i++){
            int e = V[v].adjacent_edges[i];
            int v1,v2;
            v1 = E[e].edge_vertex1;
            v2 = E[e].edge_vertex2;
            if(v1==(_V_-1)) { E[e].edge_vertex1 = v; replace_item<int>(V[v2].adjacent_vertices,(_V_-1),v); }
            if(v2==(_V_-1)) { E[e].edge_vertex2 = v; replace_item<int>(V[v1].adjacent_vertices,(_V_-1),v); }
        }
      }

      // Modify graph
      V.pop_back();
      _V_--;
      

    }else{
      std::cout<<"Nothing to delete. No such vertex with index = "<<v<<"\n";
    }
    return _V_;
}


template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::REDUCE_GRAPH(){    
 
   for(int i=0;i<_V_;i++){     
     if(V[i].vertex_degree <= 1){ this->DELETE_VERTEX(i); i = -1; }
   }
   
   return _V_;
}

template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::FIND_PATHS_FLOYD_WARSHALL(double** D,vector<Path>** P){


    // Step 1 - Create connectivity matrix and path matrix
//    vector<Path>** P;   - matrix of paths
//    int** D;            - distance matrix
    double INF = 1000*_V_;

    // Allocate memory
    D = new double*[_V_];
    P = new vector<Path>*[_V_];

    for(int i=0;i<_V_;i++){ 
       D[i] = new double[_V_]; 
       P[i] = new vector<Path>[_V_];

    }
    // Initialize "infinities"
    for(int i=0;i<_V_;i++){
        for(int j=0;j<_V_;j++){
            if(i==j){ D[i][j] = 0.0;}
            else{
            D[i][j] = INF;
            }
        }
    }
    // Initialize connectivity
    for(int e=0;e<_E_;e++){
        int i = E[e].edge_vertex1;
        int j = E[e].edge_vertex2;
        D[i][j] = D[j][i] = E[e].edge_weight;
        Path p;
        if(P[i][j].size()>0){ P[i][j].clear(); }
        if(P[j][i].size()>0){ P[j][i].clear(); }
        p.push_back(e);
        P[i][j].push_back(p);
        P[j][i].push_back(p);
        
    }

    // Setp 2 - Search for shortest paths using algorithm 1
    // Floyd - Warshall algorithm
    for(int k=0;k<_V_;k++){
        for(int i=0;i<_V_;i++){
            for(int j=0;j<_V_;j++){

                double dikj = D[i][k] + D[k][j];
                // New shortest path
                if(D[i][j]>dikj){
                    D[i][j] = dikj;       
                    Path res;                    
                    merge_paths(res,P[i][k][0],P[k][j][0]);
                    
                    if(P[i][j].size()>0) { P[i][j].clear(); }
                    P[i][j].push_back(res);
                    
                }// if            
                // Another shortest path
                else if(D[i][j]==dikj) {
                    Path new_path;
                    if((P[i][k].size()*P[k][j].size())>0){
                    merge_paths(new_path,P[i][k][0],P[k][j][0]);
                        if(new_path.size()>0){
                            P[i][j].push_back(new_path);
                        }
                    }

                }// else if

                
            }// for j
        }// for i
    }// for k
/*
//------- Test ---------
  cout<<"The distance matrix is: "<<endl;
  for(int i=0;i<_V_;i++){
     for(int j=0;j<_V_;j++){
         cout<<"D["<<i<<"]["<<j<<"] = "<<D[i][j]<<"  ";
     }
     cout<<endl;
  }

  cout<<"The paths matrix is: "<<endl;
  for(int i=0;i<_V_;i++){
     for(int j=0;j<_V_;j++){
         cout<<"P["<<i<<"]["<<j<<"] = "; show_paths(P[i][j]); cout<<"  ";
     }
     cout<<endl;
  }
*/
 
       
    return 0;
}


template <class VERTEX_DATA,class EDGE_DATA> int GRAPH<VERTEX_DATA,EDGE_DATA>::FIND_SSSR(vector<Path>& sssr){

// Reference on algorithm: Lee, C.J.; Kang, Y-M.; Cho, K-H.; Tai No, K. "A robust method for searching the smallest
// set of smallest rings with a path-included distance matrix" PNAS, 2009, 106, 17355-17358

    // Step 1 - Create connectivity matrix and PID matrices
    vector<Path>** P;
    vector<Path>** Pp;
    double** D;
    int** Cnum;
    double INF = 1000*_V_;

    // Allocate memory
    D = new double*[_V_];
    P = new vector<Path>*[_V_];
    Pp = new vector<Path>*[_V_];
    Cnum = new int*[_V_];

    for(int i=0;i<_V_;i++){ 
       D[i] = new double[_V_]; 
       P[i] = new vector<Path>[_V_];
       Pp[i] = new vector<Path>[_V_];
       Cnum[i] = new int[_V_];
    }
    // Initialize "infinities"
    for(int i=0;i<_V_;i++){
        for(int j=0;j<_V_;j++){
            Cnum[i][j] = 0;
            if(i==j){ D[i][j] = 0.0;}
            else{
            D[i][j] = INF;
            }
        }
    }
    // Initialize connectivity
    for(int e=0;e<_E_;e++){
        int i = E[e].edge_vertex1;
        int j = E[e].edge_vertex2;
        D[i][j] = D[j][i] = E[e].edge_weight;
        Path p;
        if(P[i][j].size()>0){ P[i][j].clear(); }
        if(P[j][i].size()>0){ P[j][i].clear(); }
        p.push_back(e);
        P[i][j].push_back(p);
        P[j][i].push_back(p);
        
    }

    // Setp 2 - Search for shortest paths using algorithm 1
    // Basically this is a Floyd - Warshall algorithm
    for(int k=0;k<_V_;k++){
        for(int i=0;i<_V_;i++){
            for(int j=0;j<_V_;j++){

                double dikj = D[i][k] + D[k][j];
                // New shortest path
                if(D[i][j]>dikj){
                    if(D[i][j]==(dikj+1.0)){
                        Pp[i][j] = P[i][j];
                    }// if
                    else{
                        if(Pp[i][j].size()>0) { Pp[i][j].clear(); }
                     // Pp <- 0 (empty set)
                    }
                    D[i][j] = dikj;       
                    Path res;                    
                    merge_paths(res,P[i][k][0],P[k][j][0]);
                    
                    if(P[i][j].size()>0) { P[i][j].clear(); }
                    P[i][j].push_back(res);
                    
                }// if            
                // Another shortest path
                else if(D[i][j]==dikj) {
                    Path new_path;
                    if((P[i][k].size()*P[k][j].size())>0){
                        merge_paths(new_path,P[i][k][0],P[k][j][0]);
                        if(new_path.size()>0){
                            P[i][j].push_back(new_path);
                        }
                    }

                }// else if
                else if(D[i][j]==(dikj-1.0)) {
                    Path new_path;
                    if((P[i][k].size()*P[k][j].size())>0){
                        merge_paths(new_path,P[i][k][0],P[k][j][0]);
                        if(new_path.size()>0){
                            Pp[i][j].push_back(new_path);
                        }
                    }

                }// else if

                
            }// for j
        }// for i
    }// for k

//------- Test ---------
/*
  for(i=0;i<_V_;i++){
     for(int j=0;j<_V_;j++){
         cout<<"D["<<i<<"]["<<j<<"] = "<<D[i][j]<<"  ";
     }
     cout<<endl;
  }

  for(i=0;i<_V_;i++){
     for(int j=0;j<_V_;j++){
         cout<<"P["<<i<<"]["<<j<<"] = "; show_paths(P[i][j]); cout<<"  ";
     }
     cout<<endl;
  }

  for(i=0;i<_V_;i++){
     for(int j=0;j<_V_;j++){
         cout<<"Pp["<<i<<"]["<<j<<"] = "; show_paths(Pp[i][j]); cout<<"  ";
     }
     cout<<endl;
  }
*/

    // Step 3 - Ring candidate search
    for(int i=0;i<_V_;i++){
        for(int j=0;j<_V_;j++){
            if( (D[i][j]==0.0)|| ( (P[i][j].size()==1) && (Pp[i][j].size()==0) ) ) {;;}
            else{
                if(Pp[i][j].size()!=0){  Cnum[i][j] = 2*int(D[i][j]) + 1; }
                else{                    Cnum[i][j] = 2*int(D[i][j]);     }                
            }
//            cout<<"Cnum["<<i<<"]["<<j<<"] = "<<Cnum[i][j]<<"  ";
        }
//        cout<<endl;
    }
    
    // Step 4 - Finding SSSR

//    vector<Path> sssr;
    if(sssr.size()>0) { sssr.clear(); }
    int n_ringidx = 0;
    int n_sssr; 
    int is_odd;
    int is_in_sssr;

    n_sssr = _E_ - _V_ + 1;

//    cout<<"Theoretical number of rings = "<<n_sssr<<endl;
    int max_cnum = 0;
    for(int i=0;i<_V_;i++){
        for(int j=0;j<_V_;j++){
            if(Cnum[i][j]>max_cnum) { max_cnum = Cnum[i][j]; }
        }
    }
    max_cnum++;
//    cout<<"Max cnum = "<<max_cnum<<endl;

    for(int cnum=1;cnum<max_cnum;cnum++){  // Effective ordering by Cnum

    for(int i=0;i<_V_;i++){
        for(int j=0;j<_V_;j++){            
            if(Cnum[i][j]==cnum){

                is_odd = Cnum[i][j] - 2*floor(0.5*Cnum[i][j]); 

                if(is_odd){  // Odd ring create

                   for(int k=0;k<Pp[i][j].size();k++){
                       Path new_path;
                       merge_paths(new_path,P[i][j][0],Pp[i][j][k]);
                       is_in_sssr = is_included(new_path,sssr);
                       
                       if(!is_in_sssr){ 
//                          cout<<"Add ring\n";
//                          show_path(new_path);
//                          cout<<endl;
                          sssr.push_back(new_path); n_ringidx++; 
                       }
                       
                   }// for k
                    
                }
                else{        // Even ring create

                   for(int k=0;k<(P[i][j].size()-1);k++){
                       Path new_path;
                       merge_paths(new_path,P[i][j][k],P[i][j][k+1]);
/*
                       cout<<"cnum = "<<cnum<<" i = "<<i<<" j = "<<j<<" k = "<<k<<endl;
                       cout<<"Merging paths";
                       show_path(P[i][j][k]);
                       cout<<"and \n";
                       show_path(P[i][j][k+1]);
                       cout<<"Resulting path is: \n";
                       show_path(new_path);
*/

                       is_in_sssr = is_included(new_path,sssr);

//                       cout<<"is it already in sssr? "<<is_in_sssr<<"\n";

                       if(!is_in_sssr){
//                          cout<<"Add ring\n";
//                          show_path(new_path);
//                          cout<<endl;
                          sssr.push_back(new_path); n_ringidx++;
                       }

                   }// for k



                }// else - even


            }// if Cnum[i][j]==cnum
        }// for j
    }// for i
    }// for cnum

    cout<<"Ring search completed. Number of rings found is "<<n_ringidx<<endl;   

    // Free memory
    for(int i=0;i<_V_;i++){
       delete [] D[i]; 
       delete [] P[i];
       delete [] Pp[i]; 
       delete [] Cnum[i];
    }

    delete [] D;
    delete [] P;
    delete [] Pp;
    delete [] Cnum;
   
    return 0;
}

}// namespace libgraph
}// namespace liblibra

#endif  // GRAPH_H
