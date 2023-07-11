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
  \file NList.cpp
  \brief The file implements auxiliary and main functions for the neighbor list construction
    
*/

#include "NList.h"

/// liblibra namespace
namespace liblibra{


/// libcell namespace
namespace libcell{



MATRIX fold_coords(MATRIX& R, MATRIX3x3& Box, std::string pbc_type){
/** 
  \param[in] R  Coordinates of all atoms: ndof x nsteps
  \param[in] box  The periodic cell in the format: (tv1, tv2, tv3)
  \param[in] pbc_type The parameter controlling the periodicity (when and if folding) of the unit cell
  Can take values: "a", "b", "c", "ab", "ac", "bc", and "abc"

*/

  int ndof = R.n_rows;
  int nat = ndof/3;
  int nsteps = R.n_cols;

  MATRIX res(ndof, nsteps); // results
  MATRIX3x3 invBox; invBox = Box.inverse();

  for(int step=0; step<nsteps; step++){

    for(int i=0;i<nat;i++){

      VECTOR r; 
      r.x = R.get(3*i, step);  r.y = R.get(3*i+1, step);  r.z = R.get(3*i+2, step);
      r = invBox*r;

      if((pbc_type=="a")||(pbc_type=="ab")||(pbc_type=="ac")||(pbc_type=="abc")) {
        r.x = r.x - floor(r.x);
      }
      if((pbc_type=="b")||(pbc_type=="ab")||(pbc_type=="bc")||(pbc_type=="abc")) {
        r.y = r.y - floor(r.y);
      }
      if((pbc_type=="c")||(pbc_type=="ac")||(pbc_type=="bc")||(pbc_type=="abc")) {
        r.z = r.z - floor(r.z);
      }
      r = Box * r;

      res.set(3*i, step, r.x);  res.set(3*i+1, step, r.y);  res.set(3*i+2, step, r.z);

    }// for i - all atoms
  }// for step - all data points

  return res;

}



void max_vector(VECTOR& t1,VECTOR& t2,VECTOR& t3,VECTOR& T){
/**
  \brief Finding the maximal-length vector

  This function finds the maximal projections on X,Y and
  Z coordinate axes. In absolute value!
  in the box H = (t1,t2,t3)

  \param[in] t1 unit cell vector (a direction)
  \param[in] t2 unit cell vector (b direction)
  \param[in] t3 unit cell vector (c direction)
  \param[out] T the resulting vector containing the projections of the maximal-length vector

*/

  VECTOR t;
  double x,y,z,maxx,maxy,maxz; 
  T = 0.0; maxx = maxy = maxz = 0.0;

  for(int n1=-1;n1<=1;n1++){
    for(int n2=-1;n2<=1;n2++){
      for(int n3=-1;n3<=1;n3++){
        t = n1*t1 + n2*t2 + n3*t3;
        x = fabs(t.x); if(x>=maxx){ maxx = x; T.x = x; }
        y = fabs(t.y); if(y>=maxy){ maxy = y; T.y = y; }
        z = fabs(t.z); if(z>=maxz){ maxz = z; T.z = z; }        
      }// for n3
    }// for n2
  }// for n1

}

VECTOR max_vector(VECTOR t1,VECTOR t2,VECTOR t3){

  VECTOR T;
  max_vector(t1,t2,t3, T);
  return T;
}

void apply_pbc(MATRIX3x3& H,int sz, VECTOR* in, VECTOR* out,vector<quartet>& T){
/**
  \brief Apply the periodic boundary conditions to the array of vectors

  This function decomposes the array of vectors <in> of size <sz>
  into array of vectors <out>, such that they all belong to the 
  box given by matrix <H> and positioned on global origin (0,0,0)
  and on the corresponding integer translation <T>, which are then
  stored. Combination of <T> and <out> contains all the information
  about original vectors <in>
  
  \param[in] H The matrix containing the parameters of the unit cell
  \param[in] sz The number of atoms to be processed
  \param[in] in The pointer to the array of input atomic coordinates (unfolded)
  \param[out] out The pointer to the array of the output atomic coordinates (folded)
  \param[out] the array of integer vectors(in units of the cell vectors) that fold each atom into the central unit cell
*/

  VECTOR t1,t2,t3,g1,g2,g3,res;
  quartet qt;
  H.get_vectors(t1,t2,t3);
  H.inverse().T().get_vectors(g1,g2,g3);

  if(T.size()>0){ T.clear(); }

  for(int i=0;i<sz;i++){
    qt.n1 = floor(g1*in[i]);
    qt.n2 = floor(g2*in[i]);
    qt.n3 = floor(g3*in[i]);

    out[i] = in[i] - (qt.n1*t1 + qt.n2*t2 + qt.n3*t3);
    T.push_back(qt);
  }

}

boost::python::list apply_pbc(MATRIX3x3 H, boost::python::list in, boost::python::list t){
  int i;
  vector<quartet> T;
  int sz1 = boost::python::len(t);
  for(i=0;i<sz1;i++){
    quartet tx;
    tx.is_central = boost::python::extract<int>(t[i][0]);
    tx.j = boost::python::extract<int>(t[i][1]);
    tx.n1 = boost::python::extract<int>(t[i][2]);
    tx.n2 = boost::python::extract<int>(t[i][3]);
    tx.n3 = boost::python::extract<int>(t[i][4]);

    T.push_back(tx);
  }

  int sz2 = boost::python::len(in);
  VECTOR* r; r = new VECTOR[sz2];
  VECTOR* rp; rp = new VECTOR[sz2];

  for(i=0;i<sz2;i++){  r[i] = boost::python::extract<VECTOR>(in[i]);  }

  apply_pbc(H, sz2, r, rp, T);

  boost::python::list res;
  for(i=0;i<sz2;i++){  res.append(rp[i]); }


  delete [] r;  
  delete [] rp;

  return res;

}


void serial_to_vector(int c,int Nx,int Ny,int Nz,int& nx,int& ny,int& nz){
/**
  \brief Mapping of a single integer index (serial) into 3 integers (vector)

  This function calculates the indexes of projections <nx>, <ny>, <nz>
  of given cell (given by its serial index <c>) in the 3-D array of such
  cells. The array is rectangular. Indexing goes as:
  nx in [0, Nx]
   ny in [0, Ny]
    nz in [0, Nz]
   c <---> nx,ny,nz
  <Nx>,<Ny> and <Nz> - are maximal number of cells in  corresponding
  dimension

  \param[in] c the serial index of the cell in the array of cells
  \param[in] Nx maximal number of cell replicas in x dimension (along a vector, actually)
  \param[in] Ny maximal number of cell replicas in y dimension (along b vector, actually)
  \param[in] Nz maximal number of cell replicas in z dimension (along c vector, actually)
  \param[out] nx the x (first) component of the integer translation vector for given cell
  \param[out] ny the y (second) component of the integer translation vector for given cell
  \param[out] nz the z (third) component of the integer translation vector for given cell

*/
/*
  double m1,m2;
  m2 = Nz;
  m1 = m2*Ny;

  nx = floor(c / m1);
  ny = floor( (c - nx*m1) / m2);
  nz = c - nx * m1 - ny * m2;
*/

// Probably it is better and safer to rewrite this mapping this way:
// Keep in mind we use integer division!!!

  nx = c/(Ny*Nz);  c = c - nx*(Ny*Nz); 
  ny = c/Nz;       
  nz = c - ny*Nz;

}

boost::python::list serial_to_vector(int c,int Nx,int Ny,int Nz){

  int nx,ny,nz;
  serial_to_vector(c,Nx,Ny,Nz,nx,ny,nz);

  boost::python::list res;
  res.append(nx);
  res.append(ny);
  res.append(nz);

  return res;
}

void serial_to_vector_symm(int c,int Nx,int Ny,int Nz,int& nx,int& ny,int& nz){
/**
  \brief Mapping of a single integer index (serial) into 3 integers (vector)

  This function calculates the indexes of projections <nx>, <ny>, <nz>
  of given cell (given by its serial index <c>) in the 3-D array of such
  cells. The array is rectangular. Indexing goes as:
  nx in [-Nx, Nx]
   ny in [-Ny, Ny]
    nz in [-Nz, Nz]
   c <---> nx,ny,nz
  <Nx>,<Ny> and <Nz> - are maximal number of cells in  corresponding
  dimension
  This function is similar to serial_to_vector functions, but is designed 
  for symmetric range of nx,ny and nz

  \param[in] c the serial index of the cell in the array of cells
  \param[in] Nx maximal number of cell replicas in x dimension (along a vector, actually)
  \param[in] Ny maximal number of cell replicas in y dimension (along b vector, actually)
  \param[in] Nz maximal number of cell replicas in z dimension (along c vector, actually)
  \param[out] nx the x (first) component of the integer translation vector for given cell
  \param[out] ny the y (second) component of the integer translation vector for given cell
  \param[out] nz the z (third) component of the integer translation vector for given cell

*/

  double m1,m2;
  m2 = (2*Nz+1);
  m1 = m2*(2*Ny+1);

  nx = floor(c / m1);
  ny = floor( (c - nx*m1) / m2);
  nz = c - nx * m1 - ny * m2;

  nx -= Nx;
  ny -= Ny;
  nz -= Nz;

}


boost::python::list serial_to_vector_symm(int c,int Nx,int Ny,int Nz){

  int nx,ny,nz;
  serial_to_vector_symm(c,Nx,Ny,Nz,nx,ny,nz);

  boost::python::list res;
  res.append(nx);
  res.append(ny);
  res.append(nz);

  return res;
}



void form_neibc(int c,vector<int>& neibc,int Nx,int Ny,int Nz,double cellx,double celly,double cellz,double Roff){  
/**
  \brief Form the neighbor cells list

  This function calculates serial indexes of the cells (neighbor cells)
  of the cells which are neighbor to the one with serial index <c>
  Cell a is defined as neighbor to cell b if at least in one dimension
  (x,y or z) the cells (any same poins) are less then <Roff> apart
  Parameters <cellx>,<celly> and <cellz> give the size of the sub-cell
  in corresponding dimension. 
  <Nx>,<Ny> and <Nz> - are maximal number of cells in  corresponding
  dimension
  Serial indexes are stored in <neibc> variable

  \param[in] c Index of the cell for which we want to construct neighbor cells list
  \param[out] neibc Indices (serial) of the cells that are neighbor to the cell c
  \param[in] Nx maximal number of cells in the a direction
  \param[in] Ny maximal number of cells in the b direction
  \param[in] Nz maximal number of cells in the c direction
  \param[in] cellx the size of the cell in a direction
  \param[in] celly the size of the cell in b direction
  \param[in] cellz the size of the cell in c direction
  \param[in] Roff cutoff distance that is used to determine of the atoms are neighbor

*/

  if(neibc.size()>0){ neibc.clear(); }

  int nx,ny,nz;
  serial_to_vector(c,Nx,Ny,Nz,nx,ny,nz);

  int ax,ay,az; // how many cells in each direction are necessary to exceed Roff in corresponding direction
/*
  ax = floor(Roff/cellx);
  ay = floor(Roff/celly);
  az = floor(Roff/cellz);
  if((Roff - ax*cellx) > 0.0) {  ax++; }
  if((Roff - ay*celly) > 0.0) {  ay++; }
  if((Roff - az*cellz) > 0.0) {  az++; }
*/
  // Rewriting the above lines in a more efficient and clear manner:
  ax = ay = az = 1;
  while(Roff > ax*cellx){ ax++; }
  while(Roff > ay*celly){ ay++; }
  while(Roff > az*cellz){ az++; }


  int lx,ly,lz,rx,ry,rz; // left and right limits in all directions

  lx = (nx-ax); lx = (lx<0)?0:lx;
  ly = (ny-ay); ly = (ly<0)?0:ly;
  lz = (nz-az); lz = (lz<0)?0:lz;

  rx = (nx+ax); rx = (rx>=Nx)?(Nx-1):rx;
  ry = (ny+ay); ry = (ry>=Ny)?(Ny-1):ry;
  rz = (nz+az); rz = (rz>=Nz)?(Nz-1):rz;

  int ix,iy,iz,nc;
  for(ix=lx;ix<=rx;ix++){
    for(iy=ly;iy<=ry;iy++){
      for(iz=lz;iz<=rz;iz++){

        nc = Nz*Ny*ix + Nz*iy + iz;
        neibc.push_back(nc);

      }// iz
    }// iy
  }// ix

}

boost::python::list form_neibc(int c,int Nx,int Ny,int Nz,double cellx,double celly,double cellz,double Roff){  

  vector<int> neibc;
  form_neibc(c, neibc, Nx, Ny, Nz, cellx, celly, cellz, Roff);
  boost::python::list res;

  for(int i=0;i<neibc.size();i++){ res.append(neibc[i]);  }

  return res;
}


void find_min_shell(VECTOR& t1,VECTOR& t2,VECTOR& t3,
                    VECTOR& g1,VECTOR& g2,VECTOR& g3,double Roff,
                    triple& minb,triple& maxb){
/**
  \brief Compute a minimal shell around given central unit cell

  This function calculates the minimal shell for the given cell shape
  which will satisfy the condition:
  n1*t1.x + n2*t2.x + n3*t3.x +/-  |maxT.x|  >= Roff
  n1*t1.y + n2*t2.y + n3*t3.y +/-  |maxT.y|  >= Roff
  n1*t1.z + n2*t2.z + n3*t3.z +/-  |maxT.z|  >= Roff
  minb,maxb - are the triples, corresponding to minimal
  (left,down) and maximal (upper,right) boundaries

  This is basically the way to construct the all the replicas of the original cell around it
  such that the most far distant point of the replicated cells is not closer than by Roff 
  to any point (even the closest one) of the original cell.

  \param[in] t1 unit cell vector in the a direction
  \param[in] t2 unit cell vector in the b direction
  \param[in] t3 unit cell vector in the c direction
  \param[in] g1 vector reciprocal to t1
  \param[in] g2 vector reciprocal to t2
  \param[in] g3 vector reciprocal to t3
  \param[in] Roff cutoff distance
  \param[out] minb is the triple containing the integer translations of the unit cell to form the minimal boundary of the super-cell
  \param[out] maxb is the triple containing the integer translations of the unit cell to form the maximal boundary of the super-cell

*/
  double d1,d2,d3;
  int n1,n2,n3,m1,m2,m3;
  minb.n1 = 0; minb.n2 = 0; minb.n3 = 0;
  maxb.n1 = 0; maxb.n2 = 0; maxb.n3 = 0;
         
// New, elegant, correct and easy way!
// Triply-nested loop to count all vertexes of the cell
  VECTOR zero(0.0,0.0,0.0);
  triple tzero; tzero.n1 = tzero.n2 = tzero.n3 = 0;
  vector<VECTOR> cell_corners(8,zero);
  vector<triple> cell_corners_int(8,tzero);
  vector<VECTOR> cube_corners(8,zero);
  int v1,v2,v3,cnt1,cnt2;

  cnt1=0;
  for(v1=0;v1<=1;v1++){
    for(v2=0;v2<=1;v2++){
      for(v3=0;v3<=1;v3++){
        cell_corners[cnt1] = (v1*t1 + v2*t2 + v3*t3);
        cell_corners_int[cnt1].n1 = v1;
        cell_corners_int[cnt1].n2 = v2;
        cell_corners_int[cnt1].n3 = v3;
        cnt1++;
      }// v3
    }// v2
  }// v1

  cnt2=0;
  VECTOR a(Roff,0.0,0.0),b(0.0,Roff,0.0),c(0.0,0.0,Roff);
  for(v1=-1;v1<=1;v1+=2){
    for(v2=-1;v2<=1;v2+=2){
      for(v3=-1;v3<=1;v3+=2){
        cube_corners[cnt2] = (v1*a + v2*b + v3*c);
        cnt2++;
      }// v3
    }// v2
  }// v1

  for(int i=0;i<cnt1;i++){
    for(int j=0;j<cnt2;j++){ 
      //
      VECTOR P = (cell_corners[i]+cube_corners[j]);
      n1 = floor(g1*P);//+cell_corners_int[i].n1; 
      n2 = floor(g2*P);//+cell_corners_int[i].n2;
      n3 = floor(g3*P);//+cell_corners_int[i].n3;

      if(n1<=minb.n1){ minb.n1 = n1;}  if(n1>=maxb.n1){ maxb.n1 = n1; }
      if(n2<=minb.n2){ minb.n2 = n2;}  if(n2>=maxb.n2){ maxb.n2 = n2; }
      if(n3<=minb.n3){ minb.n3 = n3;}  if(n3>=maxb.n3){ maxb.n3 = n3; }

    }// j
  }// i

  cell_corners.clear();
  cell_corners_int.clear();
  cube_corners.clear();
        


}

void make_nlist(int Nat,VECTOR* r,MATRIX3x3& H,
                    int maxa,int maxb,int maxc,double cellx,double celly,double cellz,
                    double Roff,vector< vector<quartet> >& nlist){
/**
  \brief Create neighbor list

  Implement an improved Verlet list method (combination of Verlet list and
  linked cell list). Many ideas according to:
  Yao, Z.; Wang, J-S.; Liu, G-R. and Cheng, M "Improved neighbor list
  algorithm in molecular simulations using cell decomposition and data
  sorting method" Computer Physics Communications 2004, 161, 27-35

  This function calculates the neighbor lists for all <Nat> atoms with
  original coordinates given by array <r> in simulation cell given by
  box <H> subject to periodic boundary conditions. The parameters <maxa>,
  <maxb> and <maxc> give the number of periodic translations of original
  simulation cell in each direction. They should be chosen such that all
  real atoms interact with all other atoms and images within <Roff> distance
  Parameters <cellx>, <celly>, <cellz> give the size of the sub-cell used
  to accelerate in computations and make them scale as O(NlogN)
  Final results will then be stored in <nlist> array

  \param[in] Nat The number of atoms in the system
  \param[in] r The pointer to the array containing the coordinates of all atoms
  \param[in] H is the matrix describing the shape and size of the unit cell
  \param[in] maxa The maximal number of periodic translation of the unit cell along a axis
  \param[in] maxb The maximal number of periodic translation of the unit cell along b axis
  \param[in] maxc The maximal number of periodic translation of the unit cell along c axis
  \param[in] cellx The size of the sub-cells in x direction
  \param[in] celly The size of the sub-cells in y direction
  \param[in] cellz The size of the sub-cells in z direction
  \param[in] Roff The cutoff distance which controls the formation of the neighbor list
  \param[in,out] nlist Indices of the neighboring sub-cells for all sub-cells

*/
  int a,b,c,i;
  VECTOR t1,t2,t3,g1,g2,g3;
  VECTOR t;//,maxT;

  // Calculate constants
  double Roff2 = Roff * Roff;

  // Get cell translation vectors and their reciprocal
  // and find biggest vector in simulation cell
  H.get_vectors(t1,t2,t3);
  H.inverse().T().get_vectors(g1,g2,g3);
//  max_vector(H,maxT);
 

  // Clean the results array
  if(nlist.size()>0){ nlist.clear(); }

  vector<VECTOR> R;
  vector<quartet> initT; // initial translations
  VECTOR* tmp; tmp = new VECTOR[Nat];

  apply_pbc(H,Nat,r,tmp,initT); // now all tmp in simulation cell, no outside


  // Bounding box is based on maximal atomic displacements
  int indx = 0;
  vector<quartet> globqt;
  VECTOR minr,maxr; minr = maxr = r[0];
  for(a=-maxa;a<=maxa;a++){
    for(b=-maxb;b<=maxb;b++){
      for(c=-maxc;c<=maxc;c++){
        for(i=0;i<Nat;i++){

          t = tmp[i] + a*t1 + b*t2 + c*t3;  R.push_back(t);
          if(t.x<=minr.x){ minr.x = t.x; } if(t.x>=maxr.x){ maxr.x = t.x; }
          if(t.y<=minr.y){ minr.y = t.y; } if(t.y>=maxr.y){ maxr.y = t.y; }
          if(t.z<=minr.z){ minr.z = t.z; } if(t.z>=maxr.z){ maxr.z = t.z; }
 
          quartet locqt; 
          locqt.j = i; locqt.n1 = a; locqt.n2 = b; locqt.n3 = c;
          globqt.push_back(locqt);
          indx++;
        }// for i
      }// for c
    }// for b
  }// for a
  int sz = indx; // number of "complex atoms" = atom position + PBC translation

  // Number of partitions in corresponding direction
  VECTOR maxdr; maxdr = maxr - minr;
  int Nx = (floor(maxdr.x/cellx)+1);
  int Ny = (floor(maxdr.y/celly)+1);
  int Nz = (floor(maxdr.z/cellz)+1);
  int Ncells = Nx*Ny*Nz;
  cout<<"Number of X partitions = "<<Nx<<endl;
  cout<<"Number of Y partitions = "<<Ny<<endl;
  cout<<"Number of Z partitions = "<<Nz<<endl;
  cout<<"Number of cells = "<<Ncells<<endl;


  vector<int> dummy;
  vector<quartet> qdummy;
  vector<int> at2cell_indx(sz,-1);            // is a complex index of the given atom 
  vector< vector<int> > cell2at(Ncells,dummy);// cell2at[i] - contains complex indexes of all atoms in cell i 
  nlist = std::vector< vector<quartet> >(Nat,qdummy);


  // Calculate neighbors of each cell (sub-cell)
  // we use serial indexes of both central cell and its neighbors
  vector< vector<int> > neibc;                // indexes of neighboring cells for given cell index
  for(c=0;c<Ncells;c++){
    vector<int> neibc_c;
    form_neibc(c,neibc_c,Nx,Ny,Nz,cellx,celly,cellz,Roff);
    neibc.push_back(neibc_c);
    // Check indexes of neighbor cells
    for(int nc=0;nc<neibc_c.size();nc++){
      if(neibc_c[nc]>=Ncells){  cout<<"Hey! Potential problem\n"; }
    }
  }


  // Calculate mappings between atom indexes and cell (sub-cell) indexes
  for(i=0;i<sz;i++){
    VECTOR diff = R[i] - minr;  // position of 
    triple trp;
    trp.n1 = floor(diff.x/cellx);
    trp.n2 = floor(diff.y/celly);
    trp.n3 = floor(diff.z/cellz);
    c = Nz*Ny*trp.n1 + Nz*trp.n2 + trp.n3;

    at2cell_indx[i] = c;
    cell2at[c].push_back(i);
  }

  
  int cc = (2*maxc+1)*(2*maxb+1)*maxa + (2*maxc+1)*maxb + maxc; // index of central cell

  for(int at_indx1=0;at_indx1<Nat;at_indx1++){
    int i1  = Nat*cc+at_indx1; //complex index of real atom at_indx1
    int ci1 = at2cell_indx[i1]; // complex index of cell to which atom i belongs

    for(c=0;c<neibc[ci1].size();c++){
      int ci2 = neibc[ci1][c]; // one of the neighboring cell of cell l

      for(a=0;a<cell2at[ci2].size();a++){ // iterations over atoms in cell ci2
        int i2 = cell2at[ci2][a];
        int at_indx2 = i2 % Nat;
        VECTOR dR; dR = R[i1] - R[i2];

        if(fabs(dR.x)<=Roff){
          if(fabs(dR.y)<=Roff){
            if(fabs(dR.z)<=Roff){
              if((dR.x*dR.x+dR.y*dR.y+dR.z*dR.z)<=Roff2){
                quartet qt; qt = globqt[i2];
                qt.n1 += (initT[at_indx1].n1 - initT[at_indx2].n1);
                qt.n2 += (initT[at_indx1].n2 - initT[at_indx2].n2);
                qt.n3 += (initT[at_indx1].n3 - initT[at_indx2].n3);
                nlist[at_indx1].push_back(qt);
              }
            }//zik
          }//yik
        }//xik
      }// for k

    }// for j
  }// for i


  delete [] tmp;

}

void make_nlist_auto(int Nat,VECTOR* r,MATRIX3x3& H,
                     double cellx,double celly,double cellz,
                     double Roff,vector< vector<quartet> >& nlist){
/**
  \brief Create neighbor list (no manual maximal translation parameters)

  Implement improved Verlet list method (combination of Verlet list and
  linked cell list). Many ideas according to:
  Yao, Z.; Wang, J-S.; Liu, G-R. and Cheng, M "Improved neighbor list
  algorithm in molecular simulations using cell decomposition and data
  sorting method" Computer Physics Communications 2004, 161, 27-35

  This function calculates the neighbor lists for all <Nat> atoms with
  original coordinates given by array <r> in simulation cell given by
  box <H> subject to periodic boundary conditions. The parameters <maxa>,
  <maxb> and <maxc> give the number of periodic translations of original
  simulation cell in each direction. In this version they are determined
  automatically from the shape of the simulation cell such that all
  real atoms interact with all other atoms and images within <Roff> distance
  Parameters <cellx>, <celly>, <cellz> give the size of the sub-cell used
  to accelerate in computations and make them scale as O(NlogN)
  Final results will then be stored in <nlist> array

  \param[in] Nat The number of atoms in the system
  \param[in] r The pointer to the array containing the coordinates of all atoms
  \param[in] H is the matrix describing the shape and size of the unit cell
  \param[in] cellx The size of the sub-cells in x direction
  \param[in] celly The size of the sub-cells in y direction
  \param[in] cellz The size of the sub-cells in z direction
  \param[in] Roff The cutoff distance which controls the formation of the neighbor list
  \param[in,out] nlist Indices of the neighboring sub-cells for all sub-cells

************************************************************************/

  VECTOR t1,t2,t3,g1,g2,g3;
  VECTOR t;//,maxT;

  // Calculate constants
  double Roff2 = Roff * Roff;

  // Get cell translation vectors and their reciprocal
  // and find biggest vector in simulation cell
  H.get_vectors(t1,t2,t3);
  H.inverse().T().get_vectors(g1,g2,g3);
//  max_vector(H,maxT);
 
  vector<VECTOR> R;
  vector<quartet> initT; // initial translations
  VECTOR* tmp; tmp = new VECTOR[Nat];

  apply_pbc(H,Nat,r,tmp,initT); // now all tmp in simulation cell, no outside

  triple minb,maxb;
  find_min_shell(t1,t2,t3,g1,g2,g3,Roff,minb,maxb);
//  cout<<"a in range ["<<minb.n1<<","<<maxb.n1<<"]\n";
//  cout<<"b in range ["<<minb.n2<<","<<maxb.n2<<"]\n";
//  cout<<"c in range ["<<minb.n3<<","<<maxb.n3<<"]\n";

  int a,b,c,i;

  // Bounding box is based on maximal atomic displacements
  int indx = 0;
  vector<quartet> globqt;
  VECTOR minr,maxr; minr = maxr = r[0];
  for(a=minb.n1;a<=maxb.n1;a++){
    for(b=minb.n2;b<=maxb.n2;b++){
      for(c=minb.n3;c<=maxb.n3;c++){
        for(i=0;i<Nat;i++){

          t = tmp[i] + a*t1 + b*t2 + c*t3;  R.push_back(t);
          if(t.x<=minr.x){ minr.x = t.x; } if(t.x>=maxr.x){ maxr.x = t.x; }
          if(t.y<=minr.y){ minr.y = t.y; } if(t.y>=maxr.y){ maxr.y = t.y; }
          if(t.z<=minr.z){ minr.z = t.z; } if(t.z>=maxr.z){ maxr.z = t.z; }
 
          quartet locqt; 
          locqt.j = i; locqt.n1 = a; locqt.n2 = b; locqt.n3 = c;
          globqt.push_back(locqt);
          indx++;
        }// for i
      }// for c
    }// for b
  }// for a
  int sz = indx; // number of "complex atoms" = atom position + PBC translation

  // Number of partitions in corresponding direction
  VECTOR maxdr; maxdr = maxr - minr;
  int Nx = (floor(maxdr.x/cellx)+1);
  int Ny = (floor(maxdr.y/celly)+1);
  int Nz = (floor(maxdr.z/cellz)+1);
  int Ncells = Nx*Ny*Nz;
//  cout<<"Number of X partitions = "<<Nx<<endl;
//  cout<<"Number of Y partitions = "<<Ny<<endl;
//  cout<<"Number of Z partitions = "<<Nz<<endl;
//  cout<<"Number of cells = "<<Ncells<<endl;
//  exit(0);

  vector<int> dummy;
  vector<quartet> qdummy;
  vector<int> at2cell_indx(sz,-1);            // is a complex index of the given atom 
  vector< vector<int> > cell2at(Ncells,dummy);// cell2at[i] - contains complex indexes of all atoms in cell i 

  // Clean the results array
  if(nlist.size()>0){ nlist.clear(); }
  nlist = std::vector< vector<quartet> >(Nat,qdummy);


  // Calculate neighbors of each cell (sub-cell)
  // we use serial indexes of both central cell and its neighbors
  vector< vector<int> > neibc(Ncells,dummy);     // indexes of neighboring cells for given cell index
  for(c=0;c<Ncells;c++){
//    vector<int> neibc_c;
    form_neibc(c,neibc[c],Nx,Ny,Nz,cellx,celly,cellz,Roff);
//    neibc.push_back(neibc_c);
    // Check indexes of neighbor cells
    //for(int nc=0;nc<neibc_c.size();nc++){
    //  if(neibc_c[nc]>=Ncells){  cout<<"Hey! Potential problem\n"; }
    //}
  }


  // Calculate mappings between atom indexes and cell (sub-cell) indexes
  for(i=0;i<sz;i++){
    VECTOR diff = R[i] - minr;  // position of 
    triple trp;
    trp.n1 = floor(diff.x/cellx);
    trp.n2 = floor(diff.y/celly);
    trp.n3 = floor(diff.z/cellz);
    c = Nz*Ny*trp.n1 + Nz*trp.n2 + trp.n3;

    at2cell_indx[i] = c;
    cell2at[c].push_back(i);
  }

  
  int cc = (maxb.n3-minb.n3+1)*(maxb.n2-minb.n2+1)*(-minb.n1) + (maxb.n3-minb.n3+1)*(-minb.n2) + (-minb.n3); // index of central cell

  for(int at_indx1=0;at_indx1<Nat;at_indx1++){
    int i1  = Nat*cc+at_indx1; //complex index of real atom at_indx1
    int ci1 = at2cell_indx[i1]; // complex index of cell to which atom i belongs
    int sz1 = neibc[ci1].size();// number of neighbor cells of the cell with index ci1

    int newsize = 0;
    for(c=0;c<sz1;c++){
      int ci2 = neibc[ci1][c]; // one of the neighboring cell of cell l
      int sz2 = cell2at[ci2].size();// number of atoms in the cell with index ci2
      newsize += sz2; // total possible(max) number of neighbor complex atoms
    }
    if(nlist[at_indx1].capacity()<=newsize){  nlist[at_indx1].reserve(newsize); }


    for(c=0;c<sz1;c++){
      int ci2 = neibc[ci1][c]; // one of the neighboring cell of cell l
      int sz2 = cell2at[ci2].size();// number of atoms in the cell with index ci2
     
      for(a=0;a<sz2;a++){ // iterations over atoms in cell ci2
        int i2 = cell2at[ci2][a]; // complex intex of atom a of the cell ci2
        int at_indx2 = i2 % Nat;  // real index of atom, corresponding to the atom with complex index i2
        if(at_indx2>=at_indx1){
        VECTOR dR = R[i1] - R[i2];

        double modR = dR.x*dR.x;
        if(modR<=Roff2){
          modR += dR.y*dR.y;
          if(modR<=Roff2){
            modR += dR.z*dR.z;
            if(modR<=Roff2){
              quartet qt; qt = globqt[i2];
              qt.n1 += (initT[at_indx1].n1 - initT[at_indx2].n1);
              qt.n2 += (initT[at_indx1].n2 - initT[at_indx2].n2);
              qt.n3 += (initT[at_indx1].n3 - initT[at_indx2].n3);
              qt.is_central = 0;
              if(qt.n1==0 && qt.n2==0 && qt.n3==0) { qt.is_central = 1; }
              nlist[at_indx1].push_back(qt);
            }//zik
          }//yik
        }//xik
        }// at_indx2>=at_indx1
      }// for a

    }// for c
  }// for at_indx1


  delete [] tmp;

}


double energy(int Nat,VECTOR* r,MATRIX3x3& H,vector< vector<quartet> >& nlist){
/**
  \brief Auxiliary function to test neighbor list

  The function computes the sum of distances from all atoms to all their replicas

  \param[in] Nat The number of atoms
  \param[in] r The pointer to the array of atomic coordinates
  \param[in] H Matrix describing the size and shape of the simulation cell
  \param[in] nlist The neighbor list for all atoms
 
  The function returs the computed sum - to serve a the quantitative descriptor of summation 
*/

  VECTOR t1,t2,t3,t;
 
  // Get cell translation vectors and their reciprocal
  // and find biggest vector in simulation cell
  H.get_vectors(t1,t2,t3);

  double res = 0.0 ;
  double cnt = 0.0;
  for(int i=0;i<Nat;i++){
    int sz = nlist[i].size();
    for(int k=0;k<sz;k++){
      int j = nlist[i][k].j;
      res += (r[i] - (r[j] + nlist[i][k].n1 * t1 + nlist[i][k].n2 * t2 + nlist[i][k].n3 * t3) ).length();
      cnt += 1.0;   
    }// for k
  }// for i

  return (res/cnt);
}


void bruteforce(int Nat,VECTOR* r,MATRIX3x3& H,int maxa,int maxb,int maxc,double Roff,vector< vector<quartet> >& nlist){
/**
  \brief The function for the brute force formation of the neighbor list

  \param[in] Nat The number of atoms
  \param[in] r The pointer to the array of atomic coordinates
  \param[in] H Matrix describing the size and shape of the simulation cell
  \param[in] maxa The maximal number of periodic translation of the unit cell along a axis
  \param[in] maxb The maximal number of periodic translation of the unit cell along b axis
  \param[in] maxc The maximal number of periodic translation of the unit cell along c axis
  \param[in] cellx The size of the sub-cells in x direction
  \param[in] celly The size of the sub-cells in y direction
  \param[in] cellz The size of the sub-cells in z direction
  \param[in] Roff The cutoff distance which controls the formation of the neighbor list
  \param[in,out] nlist Indices of the neighboring sub-cells for all sub-cells
 
*/


  VECTOR t1,t2,t3,g1,g2,g3;
  VECTOR t;
 
  // Calculate constants
  double Roff2 = Roff * Roff;

  // Get cell translation vectors and their reciprocal
  // and find biggest vector in simulation cell
  H.get_vectors(t1,t2,t3);
  H.inverse().T().get_vectors(g1,g2,g3);

  // Clean the results array
  if(nlist.size()>0){ nlist.clear(); }

  vector<quartet> qdummy;
  nlist = std::vector< vector<quartet> >(Nat,qdummy);


  for(int i=0;i<Nat;i++){
    for(int j=0;j<Nat;j++){

      for(int a=-maxa;a<=maxa;a++){
        for(int b=-maxb;b<=maxb;b++){
          for(int c=-maxc;c<=maxc;c++){

            t = r[i] - (r[j] + a*t1 + b*t2 + c*t3); 
            if(t.length2()<Roff2){
              quartet qt;
              qt.j = j;
              qt.n1 = a; qt.n2 = b; qt.n3 = c;
              nlist[i].push_back(qt);

            }// if
            
          }// for c
        }// for b
      }// for a
   
    }// for j
  }// for i

}


}// namespace libcell
}// liblibra
