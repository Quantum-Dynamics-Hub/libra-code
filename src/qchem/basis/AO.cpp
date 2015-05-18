#include "AO.h"
#include "MOAO.h"

//------------------ Members of  AO class -------------------------

AO::AO(const AO& g){
  // Default parameters
                                is_x_exp = 0;
                                is_y_exp = 0;
                                is_z_exp = 0;
                                is_element = 0;
                                is_ao_shell = 0;
                                is_ao_shell_type = 0;
                                is_ao_name = 0;
  expansion_size = 0;           is_expansion_size = 1;
                                is_primitives = 0;
                                is_s_primitives = 0;
                                is_coefficients = 0;
                                is_at_indx = 0;

  if(g.is_x_exp)   {  x_exp = g.x_exp; is_x_exp = 1; }
  if(g.is_y_exp)   {  y_exp = g.y_exp; is_y_exp = 1; }
  if(g.is_z_exp)   {  z_exp = g.z_exp; is_z_exp = 1; }
  if(g.is_element) {  element = g.element; is_element = 1; }
  if(g.is_ao_shell){  ao_shell = g.ao_shell; is_ao_shell = 1; }
  if(g.is_ao_shell_type){  ao_shell_type = g.ao_shell_type; is_ao_shell_type = 1; }
  if(g.is_ao_name) {  ao_name = g.ao_name; is_ao_name = 1; }
  if(g.is_expansion_size){  expansion_size = g.expansion_size; is_expansion_size = 1;  }
  if(g.is_primitives) {
    primitives.clear();
    for(int i=0;i<g.primitives.size();i++){ primitives.push_back(g.primitives[i]);  }
    is_primitives = 1;
  }
  if(g.is_s_primitives) {
    s_primitives.clear();
    for(int i=0;i<g.s_primitives.size();i++){ s_primitives.push_back(g.s_primitives[i]);  }
    is_s_primitives = 1;
  }

  if(g.is_coefficients){
    coefficients.clear();
    for(int i=0;i<g.coefficients.size();i++){ coefficients.push_back(g.coefficients[i]);  }
    is_coefficients = 1;
  }
  if(g.is_at_indx){  at_indx = g.at_indx; is_at_indx = 1; }

}

// Assignment operator
AO& AO::operator=(const AO& g){

   // Default parameters
                                is_x_exp = 0;
                                is_y_exp = 0;
                                is_z_exp = 0;
                                is_element = 0;
                                is_ao_shell = 0;
                                is_ao_shell_type = 0;
                                is_ao_name = 0;
  expansion_size = 0;           is_expansion_size = 1;
                                is_primitives = 0;
                                is_s_primitives = 0;
                                is_coefficients = 0;
                                is_at_indx = 0;

  if(g.is_x_exp)   {  x_exp = g.x_exp; is_x_exp = 1; }
  if(g.is_y_exp)   {  y_exp = g.y_exp; is_y_exp = 1; }
  if(g.is_z_exp)   {  z_exp = g.z_exp; is_z_exp = 1; }
  if(g.is_element) {  element = g.element; is_element = 1; }
  if(g.is_ao_shell){  ao_shell = g.ao_shell; is_ao_shell = 1; }
  if(g.is_ao_shell_type){  ao_shell_type = g.ao_shell_type; is_ao_shell_type = 1; }
  if(g.is_ao_name) {  ao_name = g.ao_name; is_ao_name = 1; }
  if(g.is_expansion_size){  expansion_size = g.expansion_size; is_expansion_size = 1;  }
  if(g.is_primitives) { 
    primitives.clear(); 
    for(int i=0;i<g.primitives.size();i++){ primitives.push_back(g.primitives[i]); }
//    for(int i=0;i<g.primitives.size();i++){ primitives[i] = g.primitives[i]; }
    is_primitives = 1; 
  }
  if(g.is_s_primitives) { 
    s_primitives.clear(); 
    for(int i=0;i<g.s_primitives.size();i++){ s_primitives.push_back(g.s_primitives[i]); }
//    for(int i=0;i<g.s_primitives.size();i++){ s_primitives[i] = g.s_primitives[i]; }
    is_s_primitives = 1; 
  }

  if(g.is_coefficients){
    coefficients.clear();
    for(int i=0;i<g.coefficients.size();i++){ coefficients.push_back(g.coefficients[i]); }
//    for(int i=0;i<g.coefficients.size();i++){ coefficients[i] = g.coefficients[i]; }
    is_coefficients = 1; 
  }
  if(g.is_at_indx){  at_indx = g.at_indx; is_at_indx = 1; }

  return *this;

}


void AO::show_info(){

  std::cout<<"AO properties:"<<std::endl;

  if(is_element){ std::cout<<"element = "<<element<<std::endl; }
  if(is_ao_shell){std::cout<<"ao_shell = "<<ao_shell<<std::endl;}
  if(is_ao_shell_type){std::cout<<"ao_shell_type = "<<ao_shell_type<<std::endl;}
  if(is_ao_name){ std::cout<<"ao_name = "<<ao_name<<std::endl; }
  if(is_x_exp)   {std::cout<<"x_exp = "<<x_exp<<" unitless"<<std::endl;   }
  if(is_y_exp)   {std::cout<<"y_exp = "<<y_exp<<" unitless"<<std::endl;   }
  if(is_z_exp)   {std::cout<<"z_exp = "<<z_exp<<" unitless"<<std::endl;   }
  if(is_expansion_size){ std::cout<<"expansion_size = "<<expansion_size<<" primitives"<<std::endl;
  
    for(int i=0;i<expansion_size;i++){
 
      std::cout<<"coefficients["<<i<<"] = "<<coefficients[i]<<std::endl;
      if(is_primitives){  std::cout<<"primitives["<<i<<"] = "<<std::endl;   primitives[i].show_info();   }
      else if(is_s_primitives){ std::cout<<"s_primitives["<<i<<"] = "<<std::endl;   s_primitives[i].show_info(); }

    }// for i
  } 
  if(is_at_indx) {std::cout<<"at_indx = "<<at_indx<<endl; }

  std::cout<<std::endl;

}

double AO::Evaluate(VECTOR& pos){

  double res = 0.0;

//  cout<<"In AO:Evaluate\n";
  for(int i=0;i<expansion_size;i++){
    res += coefficients[i] * primitives[i].Evaluate(pos); 
//    cout<<"i = "<<i<<"  res = "<<res<<endl;
  }
//  cout<<"End of AO:Evaluate\n";

  return res;
}

void AO::normalize(){

  // Normalization for entire contraction
  // AO = N*summ(c_i * prim[i]) , where prim - are not normalized
  //         i
  // Calculate the norm


//  cout<<"in AO::normalize()\n";
//  cout<<"is_primitives = "<<is_primitives<<" is_s_primitives = "<<is_s_primitives<<endl;

  double res = 0.0;
  VECTOR dIdA,dIdB;
  for(int i=0;i<expansion_size;i++){
    for(int j=0;j<expansion_size;j++){
      if(is_primitives){
        res += coefficients[i] * coefficients[j] * OVERLAP_INTEGRAL(primitives[i],primitives[j],1,dIdA,dIdB);
      }
      else if(is_s_primitives){
        res += coefficients[i] * coefficients[j] * STO_OVERLAP_INTEGRAL(s_primitives[i],s_primitives[j],1,dIdA,dIdB);
      }
    }
  }

//  exit(0);

  // Normalize
  res = (1.0/sqrt(res));
//  cout<<"Normalization constant = "<<res<<endl;
//  cout<<"expansion_size = "<<expansion_size<<endl;

  // Now scale contraction coefficents c[i] -> c'[i] = N*c[i], so AO' is normalized
  // int(AO'*AO'dr) = 1
  for(int i=0;i<expansion_size;i++){
    if(is_primitives){    coefficients[i] *= (res*primitives[i].norm()); }
    else if(is_s_primitives){ coefficients[i] *= (res*s_primitives[i].norm()); }
  }

}

double AO::norm2(){

  // Calculate the normalization of entire contraction
  // AO = N*summ(c_i * prim[i]) , where prim - are not normalized
  //         i
  // Calculate the norm
  double res = 0.0; 
  VECTOR dIdA,dIdB; 
  for(int i=0;i<expansion_size;i++){
    for(int j=0;j<expansion_size;j++){
      if(is_primitives){
        res += coefficients[i] * coefficients[j] * OVERLAP_INTEGRAL(primitives[i],primitives[j],0,dIdA,dIdB);
      }
      else if(is_s_primitives){
        res += coefficients[i] * coefficients[j] * STO_OVERLAP_INTEGRAL(s_primitives[i],s_primitives[j],0,dIdA,dIdB);
      }

    }
  }

  return res;
}

double AO::norm(){

  // Calculate the normalization of entire contraction
  // AO = N*summ(c_i * prim[i]) , where prim - are not normalized
  //         i
  // Calculate the norm
  double res = 0.0;
  VECTOR dIdA,dIdB;
  for(int i=0;i<expansion_size;i++){
    for(int j=0;j<expansion_size;j++){
      if(is_primitives){
        res += coefficients[i] * coefficients[j] * OVERLAP_INTEGRAL(primitives[i],primitives[j],0,dIdA,dIdB);
      }
      else if(is_s_primitives){
        res += coefficients[i] * coefficients[j] * STO_OVERLAP_INTEGRAL(s_primitives[i],s_primitives[j],0,dIdA,dIdB);
      }

    }
  }

  return sqrt(res);
}

void AO::move(const VECTOR& R){
// Move all primitives by vector R

  for(int i=0;i<expansion_size;i++){
    *(primitives[i].G_center) += R;
  }

}

void AO::set_position(const VECTOR& R){
// Set the position (coordinate) of all primitives to be the vector R

  for(int i=0;i<expansion_size;i++){
    *(primitives[i].G_center) = R;
  }

}



AO AO::operator+(AO ob){

// Concatenating AOs

  AO res;
  if(is_at_indx && ob.is_at_indx){
    if(at_indx==ob.at_indx){  res.at_indx = at_indx; is_at_indx = 1; }
    else{ cout<<"Error: Can not concatenate AOs centered on different atomic species\n"; exit(0); }
  }

  if(   expansion_size &&    (is_primitives||is_s_primitives) &&    is_coefficients ){
    res.expansion_size = expansion_size;
    for(int i=0;i<expansion_size;i++){
      if(is_primitives){     res.primitives.push_back(primitives[i]); }
      else if(is_s_primitives){ res.s_primitives.push_back(s_primitives[i]); }
      
      res.coefficients.push_back(coefficients[i]);
    }
    res.is_expansion_size = 1;
    res.is_primitives = is_primitives;
    res.is_s_primitives = is_s_primitives;
    res.is_coefficients = 1;

  }

  if(ob.expansion_size && (ob.is_primitives||ob.is_s_primitives) && ob.is_coefficients ){
    res.expansion_size += ob.expansion_size;
    for(int i=0;i<ob.expansion_size;i++){
      if(ob.is_primitives) { res.primitives.push_back(ob.primitives[i]); }
      else if(ob.is_s_primitives){ res.s_primitives.push_back(ob.s_primitives[i]); }
      res.coefficients.push_back(ob.coefficients[i]);
    }
    res.is_expansion_size = 1;
    res.is_primitives = ob.is_primitives;
    res.is_s_primitives = ob.is_s_primitives;
    res.is_coefficients = 1;

  }
  return res;

}

void AO::operator+=(AO ob){

  if(is_at_indx && ob.is_at_indx){
    if(at_indx==ob.at_indx){}
    else{ cout<<"Error: Can not concatenate AOs centered on different atomic species\n"; exit(0); }
  }

  if(ob.expansion_size && (ob.is_primitives||ob.is_s_primitives) && ob.is_coefficients){ 
    expansion_size += ob.expansion_size;

    for(int i=0;i<ob.expansion_size;i++){
      if(ob.is_primitives){       primitives.push_back(ob.primitives[i]); }
      else if(ob.is_s_primitives){       s_primitives.push_back(ob.s_primitives[i]); }
      coefficients.push_back(ob.coefficients[i]);
    }// for i

    is_expansion_size = 1; 
    is_primitives = ob.is_primitives;
    is_s_primitives = ob.is_s_primitives;
    is_coefficients = 1;
  }// if

}

/*
AO AO::operator=(AO ob){

  if(ob.is_spin){ spin = ob.spin; is_spin = 1; }
  if(ob.expansion_size){ expansion_size = ob.expansion_size; is_expansion_size = 1; } 
  if(ob.is_ao){ ao = ob.ao; is_ao = 1; }
  if(ob.is_coefficients){ coefficients = ob.coefficients; is_coefficients = 1; }

  return *this;
}
*/



AO operator*(const double& f,  const AO& m1){

  AO res;
  res = m1;
  for(int i=0;i<res.expansion_size;i++){
    res.coefficients[i] *= f;
  }
  return res;
}


