#include "Nuclear.h"

namespace libdyn{
namespace libnuclear{


Nuclear::Nuclear(){  
// cout<<"Nuclear ctor\n";
 nnucl = 0; 
}

Nuclear::Nuclear(int _nnucl){  
//    cout<<"Nuclear ctor\n";
    nnucl = _nnucl;
    mass = vector<double>(nnucl,2000.0); 
    q = vector<double>(nnucl,0.0);
    p = vector<double>(nnucl,0.0);
    f = vector<double>(nnucl,0.0);
    ctyp = vector<int>(nnucl,0);  
}

Nuclear::Nuclear(const Nuclear& ob){
//  cout<<"Nuclear cctor\n";
  nnucl = ob.nnucl;
  mass = ob.mass;
  q = ob.q; 
  p = ob.p;
  f = ob.f;
  ctyp = ob.ctyp;
}

Nuclear::~Nuclear(){
//  cout<<"Nuclear destructor\n";
  if(mass.size()>0){  mass.clear(); }
  if(q.size()>0){  q.clear(); }
  if(p.size()>0){  p.clear(); }
  if(f.size()>0){  f.clear(); }
  if(ctyp.size()>0){  ctyp.clear(); }
}


void Nuclear::propagate_p(int i,double dt){  p[i] += dt*f[i];  }
void Nuclear::propagate_p(double dt){   for(int i=0;i<nnucl;i++){  p[i] += dt*f[i];   }  }
void Nuclear::propagate_p(double dt,vector<int>& active){  
  int sz = active.size();  for(int i=0;i<sz;i++){ int I = active[i]; p[I] += dt*f[I];   }
}

void Nuclear::scale_p(int i,double scl){  p[i] *= scl;  }
void Nuclear::scale_p(double scl){  for(int i=0;i<nnucl;i++){  p[i] *= scl;   }  }
void Nuclear::scale_p(double scl,vector<int>& active){  
  int sz = active.size(); for(int i=0;i<sz;i++){ int I = active[i]; p[I] *= scl;   } 
}



void Nuclear::propagate_q(int i,double dt){  q[i] += dt*p[i]/mass[i];  }
void Nuclear::propagate_q(double dt){  for(int i=0;i<nnucl;i++){  q[i] += dt*p[i]/mass[i];   }  }
void Nuclear::propagate_q(double dt,vector<int>& active){  
  int sz = active.size();  for(int i=0;i<sz;i++){ int I = active[i]; q[I] += dt*p[I]/mass[I];   }
}

void Nuclear::scale_q(int i,double scl){  q[i] *= scl;  }
void Nuclear::scale_q(double scl){  for(int i=0;i<nnucl;i++){  q[i] *= scl;   }  }
void Nuclear::scale_q(double scl,vector<int>& active){  
  int sz = active.size(); for(int i=0;i<sz;i++){ int I = active[i]; q[I] *= scl;   } 
}




}// namespace libnuclear
}// namespace libdyn


