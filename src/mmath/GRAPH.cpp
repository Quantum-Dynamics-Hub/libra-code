#include "GRAPH.h"
#include <vector>
#include <algorithm>

using namespace std;

namespace libmmath{


int merge_paths(Path& result, Path p1, Path p2){

   int sz1 = p1.size();
   int sz2 = p2.size();
   //cout<<"merge_paths of sizes "<<sz1<<" "<<sz2<<endl;

   if((sz1>0)&&(sz2>0)){
       result = p1;      
       for(int i=0;i<sz2;i++){
           result.push_back(p2[i]);
       }
   }
   else{
       if(sz1>0){ result = p1; }
       if(sz2>0){ result = p2; }
   }

   return 0;
}
/*
int is_included(Path& p,vector<Path>& P){
// Check if the path p is included in set of paths P
// Returns 0 (not included) or 1 (included)
   int res = 0;
   int set_size = P.size();
   int sz = p.size();
   for(int i=0;i<set_size;i++){
       // Compares sizes
       if(res==0){
       if(sz==P[i].size()){
           int are_equal = 1;
           // Check direct order
           for(int j=0;j<sz;j++){
               if(p[j]!=P[i][j]) { are_equal = 0; }
           }
           if(are_equal){
               res = 1;
           }else{
               // Check reverse order
               are_equal = 1;
               for(int j=0;j<sz;j++){
                   if(p[sz-1-j]!=P[i][j]) { are_equal = 0; }
               }// for j
               res = are_equal;
           }// else
       }// if sz==P[i].size()    
       }// if res==0
   }// for i
  
   return res;
}
*/
int path_xor(Path& p1, Path& p2, Path& res){
// res = p1 XOR p2

    if(res.size()>0) { res.clear(); }
    int sz1 = p1.size();
    int sz2 = p2.size();

    for(int i=0;i<sz1;i++){
        vector<int>::iterator it;
        it = find(p2.begin(),p2.end(),p1[i]);

        if(it==p2.end()) { res.push_back(p1[i]); }
    }

    for(int i=0;i<sz2;i++){
        vector<int>::iterator it;
        it = find(p1.begin(),p1.end(),p2[i]);

        if(it==p1.end()) { res.push_back(p2[i]); }
    }

    return 0;

}
int is_included(Path& p,Path& P){
// Check if set p is included in set P (or equal to it)

   int res = 1;
   int sz1,sz2;
   sz1 = p.size();
   sz2 = P.size();  
   
   for(int i=0;i<sz1;i++){
       vector<int>::iterator it;
       it = find(P.begin(), P.end(), p[i]);
       if(it!=P.end()) { res *= 1; }
       else{             res *= 0; }       
   }

  /*
   else{
       for(int i=0;i<sz2;i++){
           vector<int>::iterator it;
           it = find(p.begin(), p.end(), P[i]);
           if(it!=p.end()) { res = 1; }
       }
   }// else

*/
  
   return res;
}

int is_included(Path& p,vector<Path>& P){
// Check if any of P[i] or (P[i] XOR P[j]) is included in p
    int res = 0;
    for(int i=0;i<P.size();i++){
        if(is_included(P[i],p)) { res = 1; }
        for(int j=0;j<P.size();j++){
            Path tmp;
            path_xor(P[i],P[j],tmp);
            if(tmp.size()>0){
                if(is_included(tmp,p)) { res = 1; }
            }
        }
    }
    
    return res;
}


void show_path(Path& p){

   cout<<"(";
   for(int i=0;i<p.size();i++){
   cout<<p[i];
   }
   cout<<")";

}

void show_paths(vector<Path>& p){

   cout<<"[ ";
   for(int i=0;i<p.size();i++){
   show_path(p[i]);
   cout<<" ";
   }
   cout<<"]";

}

}// namespace libmmath

//====================================================

