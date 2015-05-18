#include "MO.h"
#include "MOAO.h"

//---------------------- Members of MO class -------------------------
void MO::show_info(){

  std::cout<<"MO properties:"<<std::endl;

  if(is_spin){ std::cout<<"spin = "<<spin<<std::endl; }
  if(is_expansion_size){ std::cout<<"expansion_size = "<<expansion_size<<" AOs"<<std::endl;

      for(int i=0;i<expansion_size;i++){

      std::cout<<"coefficients["<<i<<"] = "<<coefficients[i]<<std::endl;
      std::cout<<"ao["<<i<<"] = "<<std::endl;   ao[i].show_info();

      }// for i
  }

  std::cout<<std::endl;

}
/*
void MO::set(object at){

    set_value(is_spin, spin,  at,"spin");

}
*/

double MO::Evaluate(VECTOR& pos){

    double res;

    for(int i=0;i<expansion_size;i++){
        res += coefficients[i] * ao[i].Evaluate(pos);
    }

    return res;
}

void MO::normalize(){

  // Normalization for entire contraction
  // MO = N*summ(c_i * ao[i]) , where ao - are normalized
  //         i
  // Calculate the norm
  double res = 0.0;
  VECTOR dIdA,dIdB;
  for(int i=0;i<expansion_size;i++){
    for(int j=0;j<expansion_size;j++){
      res += coefficients[i] * coefficients[j] * OVERLAP_INTEGRAL(ao[i],ao[j],0,dIdA,dIdB);
    }
  }

  // Normalize
  res = (1.0/sqrt(res));
  cout<<"Normalization constant = "<<res<<endl;
  cout<<"expansion_size = "<<expansion_size<<endl;
  // Now scale contraction coefficents c[i] -> c'[i] = N*c[i], so MO' is normalized
  // int(MO'*MO'dr) = 1
  for(int i=0;i<expansion_size;i++){
    coefficients[i] *= res;//(res*ao[i].norm());
  }

}



double MO::norm2(){

  // Calculate the normalization of entire contraction
  // MO = summ(c_i * ao[i]) , where prim - are not normalized
  //        i
  // Calculate the norm
  double res = 0.0;
  VECTOR dIdA,dIdB;
  for(int i=0;i<expansion_size;i++){
    for(int j=0;j<expansion_size;j++){
      res += coefficients[i] * coefficients[j] * OVERLAP_INTEGRAL(ao[i],ao[j],0,dIdA,dIdB);
    }
  }

  return (1.0/res);

}

MO MO::operator+(MO ob){

  MO res;
  if(   expansion_size &&    is_ao &&    is_coefficients ){
    res.expansion_size = expansion_size;
    for(int i=0;i<expansion_size;i++){
      res.ao.push_back(ao[i]);
      res.coefficients.push_back(coefficients[i]);
    }
    res.is_expansion_size = 1;
    res.is_ao = 1;
    res.is_coefficients = 1;

  }

  if(ob.expansion_size && ob.is_ao && ob.is_coefficients ){
    res.expansion_size += ob.expansion_size;
    for(int i=0;i<ob.expansion_size;i++){
      res.ao.push_back(ob.ao[i]);
      res.coefficients.push_back(ob.coefficients[i]);
    }
    res.is_expansion_size = 1;
    res.is_ao = 1;
    res.is_coefficients = 1;

  }
  return res;

}

void MO::operator+=(MO ob){

  if(ob.expansion_size && ob.is_ao && ob.is_coefficients){ 
    expansion_size += ob.expansion_size; is_expansion_size = 1;

    for(int i=0;i<ob.expansion_size;i++){
      ao.push_back(ob.ao[i]);
      coefficients.push_back(ob.coefficients[i]);
    }// for i
    is_expansion_size = 1;
    is_ao = 1;
    is_coefficients = 1;

  }// if

}

/*
MO MO::operator=(MO ob){

  if(ob.is_spin){ spin = ob.spin; is_spin = 1; }
  if(ob.expansion_size){ expansion_size = ob.expansion_size; is_expansion_size = 1; } 
  if(ob.is_ao){ ao = ob.ao; is_ao = 1; }
  if(ob.is_coefficients){ coefficients = ob.coefficients; is_coefficients = 1; }

  return *this;
}
*/


MO operator*(const double& f,  const MO& m1){

  MO res;
  res = m1;
  for(int i=0;i<res.expansion_size;i++){
    res.coefficients[i] *= f;
  }
  return res;
}




