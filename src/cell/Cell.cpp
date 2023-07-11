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
  \file Cell.cpp
  \brief The file implements classes and functions for periodic calculations
    
*/

#include "Cell.h"

/// liblibra namespace
namespace liblibra{


/// libcell namespace
namespace libcell{


void Cell::brute_force(VECTOR& rij, int degree, vector<triple>& res,triple& central_translation){
/**
  \brief Brute force generation of the neighbor list for a given pair of atoms

  \param[in] rij  Vector connecting the two atoms for which we want to construct neighbor list
  \param[in] degree The maximal number of unit cells to consider in all directions (usually, if the
             cell is large and the cutoff distance is not, it may be sufficient to have degree=1). Degree = 0 
             implies only the original unit cell with no periodic re[licas. Set degree to a larger value, if
             the simulation cell is small.
  \param[out] res The list of triples with each triple describing the integer translations (for given specific pair of atoms)
             of the original cell needed to account for all neighbors
  \param[out] central_translation The triple that makes the two atoms to be minimally separated (be in the central cell). For instance, if the two
              atoms are near the opposite sides of the box, the 1 box translation will put them together (nearby). That translation is 
              then the central translation.
    
*/


  double Roff2 = Roff * Roff;
  VECTOR r;

  if(res.size()>0){  res.clear(); }    
  central_translation.n1 = -degree;
  central_translation.n2 = -degree;
  central_translation.n3 = -degree;

  r = rij;
  central_translation.n1 = 0;
  central_translation.n2 = 0;
  central_translation.n3 = 0;
  double min_dist = r.length2();
  
  for(int n1=-degree;n1<=degree;n1++){
    for(int n2=-degree;n2<=degree;n2++){
      for(int n3=-degree;n3<=degree;n3++){

        r = (rij - n1*t1 - n2*t2 - n3*t3);
        double d = r.length2(); 
        if(d<=Roff2){
          triple t;
          t.n1 = n1;
          t.n2 = n2;
          t.n3 = n3;          
          if(d<=min_dist){ central_translation = t; min_dist = d;}
          res.push_back(t);
        }
      }
    }
  }

}


void Cell::calculate(VECTOR& rij,vector<triple>& res,triple& central_translation){
/**
  \brief Function for generation of the neighbor list for a given pair of atoms. This is more efficient version than the brute
         force. Use the brute force only for verification.

  \param[in] rij  Vector connecting the two atoms for which we want to construct neighbor list
  \param[out] res The list of triples with each triple describing the integer translations (for given specific pair of atoms)
             of the original cell needed to account for all neighbors
  \param[out] central_translation The triple that makes the two atoms to be minimally separated (be in the central cell). For instance, if the two
              atoms are near the opposite sides of the box, the 1 box translation will put them together (nearby). That translation is 
              then the central translation.
    
*/

//  int db = 0;
//  if(fabs(rij.x-3.5)<0.001 && fabs(rij.y+4.033)<0.001 && fabs(rij.z-0.747)<0.001){ db = 1; }
//  if(db){ cout<<"In Cell::calculate function\n"; }

  double nm,np,delta,det,det1,det2;
  int n1,n2,n3;
  int N1_min,N1_max;
  int N2_min,N2_max;
  int N3_min,N3_max;

  central_translation.n1 = floor(g1*rij + 0.5);
  central_translation.n2 = floor(g2*rij + 0.5);
  central_translation.n3 = floor(g3*rij + 0.5);


  //Pre-calculate the parameters, which depend on rij
  A  = t1*rij;
  B  = t2*rij;
  C  = t3*rij;
  C2 = C*C;
  D = Roff*Roff - rij.length2();
  N = (B - T23*C/c2);
  Q = (A - T13*C/c2);
  R = (D+C2/c2);
  Yb = (Q*Xa + M*N);
  Yc =-(R*Xa + N*N);

  // Clean results:
  if(res.size()>0){ res.clear(); }

  det = Yb*Yb - Ya*Yc;
  if(det<0.0 && fabs(det)<1e-5){ det = 0.0; }
//  if(db){ cout<<"det = "<<det; }
  if(det >= 0.0){
    delta = sqrt(det);
    nm = (Yb - delta)/Ya;
    np = (Yb + delta)/Ya;
    N1_min = ceil(nm);
    N1_max = floor(np);
//    if(db){ cout<<" N1_min="<<N1_min<<" N1_max="<<N1_max<<endl; }

    // Now run loop in n1
    for(n1=N1_min;n1<=N1_max;n1++){
      Xb = n1*M+N;
      Xc = n1*n1*P - 2.0*Q*n1 - R;
      det1 = Xb*Xb - Xa*Xc;
      if(det1<0.0 && fabs(det1)<1e-5){ det1 = 0.0; }
      delta = sqrt(det1);
      nm = (Xb - delta)/Xa;
      np = (Xb + delta)/Xa;
      N2_min = ceil(nm);
      N2_max = floor(np);
//      if(db){cout<<"  n1="<<n1<<" N2_min="<<N2_min<<" N2_max="<<N2_max<<" delta^2="<<(Xb*Xb - Xa*Xc)<<" nm="<<nm<<" np="<<np<<" Xa="<<Xa<<endl; }
   
      // Now run loop in n2
      for(n2=N2_min;n2<=N2_max;n2++){
        Zb = (C - n1*T13 - n2*T23);
        Zc = (n1*n1*a2 - 2.0*A*n1 + n2*n2*b2 - 2.0*B*n2 - D + 2.0*n1*n2*T12);        
        det2 = Zb*Zb - Za*Zc;
        if(det2<0.0 && fabs(det2)<1e-5){ det2 = 0.0; }
        delta = sqrt(det2);
        nm = (Zb - delta)/Za;
        np = (Zb + delta)/Za;
        N3_min = ceil(nm);
        N3_max = floor(np);
//        if(db){ cout<<"    n2="<<n2<<" N3_min="<<N3_min<<" N3_max="<<N3_max<<" delta^2 = "<<(Zb*Zb - Za*Zc)<<" nm="<<nm<<" np="<<np<<" Za = "<<Za<<endl; }

        // Now run loop in n3
        for(n3=N3_min;n3<=N3_max;n3++){

          triple t;
          t.n1 = n1;
          t.n2 = n2;
          t.n3 = n3;
          if( (n1==central_translation.n1) &&
              (n2==central_translation.n2) &&
              (n3==central_translation.n3)
            ){ t.is_central = 1; }
          else{t.is_central = 0; }

          res.push_back(t);

        }// for n3
      }// for n2
    }// for n1
  }// if det>=0.0

//  central_translation.n1 = floor(g1*rij + 0.5);
//  central_translation.n2 = floor(g2*rij + 0.5);
//  central_translation.n3 = floor(g3*rij + 0.5);

}

void Cell::update_vlist(int sz,VECTOR* r,vector< vector<quartet> >& at_neib, vector<triple>& central_translation){
/**
  \brief Function to update the Verlet list

  In this case, we supply the coordinates of all atoms and for all atom we find an integer unit cell translation vector
  that corresponds to the cells that may contain atoms within the cutoff range of the considered atom.

  \param[in] sz Size of the atoms array = the number of atoms.
  \param[in] r  Pointer to the array of the atomic coordinates
  \param[out] at_neib The list of neighbor cell translations for each atom. Such translations are chosen to ensure that there are
              some atoms (replicas) in the translated cell that(atoms) are within the cutoff distance from the original (central) one.
  \param[out] central_translation The triple that makes the two atoms to be minimally separated (be in the central cell). For instance, if the two
              atoms are near the opposite sides of the box, the 1 box translation will put them together (nearby). That translation is 
              then the central translation.
    
*/


// Unconditional Update - in all cases:

  // Clean result holders
  if(at_neib.size()>0){ at_neib.clear(); }
  if(central_translation.size()>0) { central_translation.clear(); }

  // Double loop for all pairs (distinct)
  // This update procedure scales as O(N^2)
  for(int i=0;i<sz;i++){
    vector<quartet> vq;
    for(int j=i;j<sz;j++){
      VECTOR rij; rij = r[i] - r[j];
//      if(i==27 && j==79){
//        cout<<"r[27]="<<r[i]<<" r[79]="<<r[j]<<" rij = "<<rij<<endl;
//      }
      vector<triple> res; triple ct; calculate(rij,res,ct);

      if(res.size()>0){
        for(int k=0;k<res.size();k++){
          quartet q;
          q.j = j;
          q.n1 = res[k].n1;
          q.n2 = res[k].n2;
          q.n3 = res[k].n3;
          q.is_central = res[k].is_central;
          vq.push_back(q);
        }
      }
//      if(i==27){
//        cout<<"j="<<j<<" res.size()="<<res.size()<<" vq.size()="<<vq.size()<<endl;
//      }

      central_translation.push_back(ct);   
    }// for j
    at_neib.push_back(vq);
  }// for i

}




}//namespace libcell
}// liblibra









