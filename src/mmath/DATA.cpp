#include "DATA.h"

//=========================== DATA class ================================

namespace libmmath{

DATA::DATA(vector<double> d){

//  Data = d;
   if(Data.size()>0){ Data.clear(); }
   for(int i=0;i<d.size();i++){
   Data.push_back(d[i]);
   }


  is_ave = 0;
  is_var = 0;
  is_sd  = 0;
  is_se  = 0;
  is_mse = 0;
  is_rmse= 0;

  is_min = 0;
  is_min_indx = 0;
  is_max = 0;
  is_max_indx = 0;

  scale_factor = 1.0;  is_scale_factor = 1;
  shift_amount = 0.0;  is_shift_amount = 1;

}

DATA::DATA(int sz,double* d){

  if(Data.size()>0){
     Data.clear();
  }
  for(int i=0;i<sz;i++){
     Data.push_back(d[i]);
  }

  is_ave = 0;
  is_var = 0;
  is_sd  = 0;
  is_se  = 0;
  is_mse = 0;
  is_rmse= 0;

  is_min = 0;
  is_min_indx = 0;
  is_max = 0;
  is_max_indx = 0;

  scale_factor = 1.0;  is_scale_factor = 1;
  shift_amount = 0.0;  is_shift_amount = 1;


}

DATA::DATA(boost::python::list obj){

   int sz = len(obj);
   double tmp;

   if(Data.size()>0){
      Data.clear();
   }

   for(int i=0;i<sz;i++){
      tmp = boost::python::extract<double>(obj[i]);
      Data.push_back(tmp);
   }

  is_ave = 0;
  is_var = 0;
  is_sd  = 0;
  is_se  = 0;
  is_mse = 0;
  is_rmse= 0;

  is_min = 0;
  is_min_indx = 0;
  is_max = 0;
  is_max_indx = 0;

  scale_factor = 1.0;  is_scale_factor = 1;
  shift_amount = 0.0;  is_shift_amount = 1;


}


DATA::~DATA(){

   if(Data.size()>0) {  Data.clear(); }

}

DATA::DATA(const DATA& d){

// Copy constructor

  is_ave = 0;
  is_var = 0;
  is_sd  = 0;
  is_se  = 0;
  is_mse = 0;
  is_rmse= 0;

  is_min = 0;
  is_min_indx = 0;
  is_max = 0;
  is_max_indx = 0;

  scale_factor = 1.0;  is_scale_factor = 1;
  shift_amount = 0.0;  is_shift_amount = 1;




   if(Data.size()>0){ Data.clear(); }
   for(int i=0;i<d.Data.size();i++){
   Data.push_back(d.Data[i]);
   }

   if(d.is_ave){ ave = d.ave; is_ave = 1; }
   if(d.is_var){ var = d.var; is_var = 1; }
   if(d.is_sd) { sd  = d.sd;  is_sd = 1;  }
   if(d.is_se) { se  = d.se;  is_se = 1;  }
   if(d.is_mse){ mse = d.mse; is_mse = 1; }
   if(d.is_rmse){rmse= d.rmse;is_rmse=1;  }

   if(d.is_min){ min = d.min; is_min = 1; }
   if(d.is_min_indx){min_indx = d.min_indx; is_min_indx = 1; }
   if(d.is_max){ max = d.max; is_max = 1; }
   if(d.is_max_indx){max_indx = d.max_indx; is_max_indx = 1; }

   if(d.is_scale_factor){ scale_factor = d.scale_factor; is_scale_factor = 1; }
   if(d.is_shift_amount){ shift_amount = d.shift_amount; is_shift_amount = 1; }

}

DATA& DATA::operator=(const DATA& d){

  is_ave = 0;
  is_var = 0;
  is_sd  = 0;
  is_se  = 0;
  is_mse = 0;
  is_rmse= 0;

  is_min = 0;
  is_min_indx = 0;
  is_max = 0;
  is_max_indx = 0;

  scale_factor = 1.0;  is_scale_factor = 1;
  shift_amount = 0.0;  is_shift_amount = 1;




   if(Data.size()>0){ Data.clear(); }
   for(int i=0;i<d.Data.size();i++){
   Data.push_back(d.Data[i]);
   }


   if(d.is_ave){ ave = d.ave; is_ave = 1; }
   if(d.is_var){ var = d.var; is_var = 1; }
   if(d.is_sd) { sd  = d.sd;  is_sd = 1;  }
   if(d.is_se) { se  = d.se;  is_se = 1;  }
   if(d.is_mse){ mse = d.mse; is_mse = 1; }
   if(d.is_rmse){rmse= d.rmse;is_rmse=1;  }

   if(d.is_min){ min = d.min; is_min = 1; }
   if(d.is_min_indx){min_indx = d.min_indx; is_min_indx = 1; }
   if(d.is_max){ max = d.max; is_max = 1; }
   if(d.is_max_indx){max_indx = d.max_indx; is_max_indx = 1; }

   if(d.is_scale_factor){ scale_factor = d.scale_factor; is_scale_factor = 1; }
   if(d.is_shift_amount){ shift_amount = d.shift_amount; is_shift_amount = 1; }

   return *this;
}
int operator == (const DATA& d1, const DATA& d2){

   int res = (d1.Data.size()==d2.Data.size());

   if(res){  for(int i=0;i<d1.Data.size();i++){  if(d1.Data[i]!=d2.Data[i]) { res = 0; }  }    }
   return res;
}

int operator != (const DATA& d1, const DATA& d2){

   int res = (d1.Data.size()!=d2.Data.size());

   if(!res){ for(int i=0;i<d1.Data.size();i++){  if(d1.Data[i]!=d2.Data[i]){ res = 1; }   }    }
   return res;
}



int DATA::PutData(vector<double>& out){

  out = Data;

  return 0;
}

int DATA::GetData(vector<double>& inp){

  Data = inp;

  return 0;
}
int DATA::Calculate_Estimators(double& Ave, double& Var, double& Sd,double& Se, double& Mse, double& Mae,double& Rmse){

  ave = 0.0;
  var = 0.0;
  mse = 0.0;
  mae = 0.0;

  double tmp;
  double sz = Data.size();

  if(sz>1){


  for(int i=0;i<sz;i++){
     ave += Data[i];
  }
  ave = (ave/sz);

  for(i=0;i<sz;i++){
     tmp = (Data[i]-ave);
     var += tmp*tmp;
     mae += fabs(tmp);
  }

  mse = var;
  var = (var/(sz-1.0));
  mse = (mse/sz);
  mae = (mae/sz);

  sd  = sqrt(var);
  se  = sd/sqrt(sz);

  rmse = sqrt(mse);

  Ave  = ave;
  Var  = var;
  Sd   = sd;
  Se   = se;
  Mse  = mse;
  Mae  = mae;
  Rmse = rmse;

  is_ave = 1;
  is_var = 1;
  is_sd  = 1;
  is_se  = 1;
  is_mse = 1;
  is_mae = 1;
  is_rmse= 1;

  }else{
   std::cout<<"Warning: Data size is less than 2. Can not calculate statistics for such data\n";
  }


  return 0;
}

int DATA::Calculate_MiniMax(double& Min_Val,int& Min_Indx,double& Max_Val,int& Max_Indx){

  int sz = Data.size();


  min = max = Data[0];
  min_indx = max_indx = 0;

  for(int i=1;i<sz;i++){
      if(Data[i]<min) { min = Data[i]; min_indx = i; }
      if(Data[i]>max) { max = Data[i]; max_indx = i; }
  }

  Min_Val  = min;
  Min_Indx = min_indx;

  Max_Val  = max;
  Max_Indx = max_indx;

  is_min      = 1;
  is_max      = 1;
  is_min_indx = 1;
  is_max_indx = 1;

  return 0;
}

int DATA::Calculate_Distribution(vector<double>& Interval,vector<double>& Density,vector<double>& Cumulative){

  int sz = Data.size();
  int szi= Interval.size()-1; // Number of sub-intervals

  // Clear results
  if(Density.size()>0){ Density.clear();}
  if(Cumulative.size()>0){ Cumulative.clear();}

  // Init results
  for(int i=0;i<szi;i++){
     Density.push_back(0.0);
     Cumulative.push_back(0.0);
  }

  // Count
  for(i=0;i<sz;i++){

     for(int j=0;j<szi;j++){
         if( (Interval[j]<=Data[i])&&(Data[i]<Interval[j+1]) ){
              Density[j]++;
              for(int k=j;k>0;k--){
              Cumulative[k]++;
              }// for k
         }// if
     }// for j
  }// for i

  return 0;
}
int DATA::LinearTransformData(double sc_fact,double sh){

  int sz = Data.size();
  scale_factor = sc_fact*scale_factor;
  shift_amount = sc_fact*shift_amount + sh;

  for(int i=0;i<sz;i++){
      Data[i] = Data[i]*sc_fact + sh;
  }

  // Change corresponding estimators
  if(is_ave){  ave = sc_fact * ave + sh; }
  if(is_var){  var = sc_fact * sc_fact * var;      }
  if(is_sd) {  sd  = sc_fact * sd;   }
  if(is_se) {  se  = sc_fact * se;   }
  if(is_mse){  mse = sc_fact * sc_fact * mse; }
  if(is_rmse){ rmse= sc_fact * rmse; }
  if(is_min&&is_max&&is_min_indx&&is_max_indx){
     if(sc_fact>=0.0) { min = min*sc_fact + sh;
                        max = max*sc_fact + sh;
                      }
     else{              double tmp;
                        int itmp;
                        tmp = min;
                        itmp = min_indx;
                        min = sc_fact*max + sh;
                        max = sc_fact*tmp + sh;
                        min_indx = max_indx;
                        max_indx = itmp;
     }// else
  }// if


  return 0;
}
int DATA::ScaleData(double sc_fact){

  this->LinearTransformData(sc_fact,0.0);

  return 0;
}

int DATA::ScaleData(double low, double upp){

  if(is_min&&is_max){
  }else{
   double a,b;
   int i,j;
   this->Calculate_MiniMax(a,i,b,j);
  }
  if(min==max){
    // Do nothing since in this case this operation is not defined
  }else{
    // Do rescaling and shifting
    double sc_factor = ((upp-low)/(max-min));
    double sh = (low - min*sc_factor);

    this->LinearTransformData(sc_factor,sh);

  }

  return 0;
}

int DATA::ShiftData(double sh){

  this->LinearTransformData(1.0,sh);

  return 0;
}
int DATA::NormalizeData(){

  double Ave,Var,Sd,Se,Mse,Mae,Rmse;
  if(is_ave&&is_sd){
  }else{
   this->Calculate_Estimators(Ave,Var,Sd,Se,Mse,Mae,Rmse);
  }

  if(Sd!=0.0){
  this->LinearTransformData((1.0/Sd),(-Ave/Sd));
  }

  return 0;
}


/*
int DATA::Lin_Regression(int flag,double& a,double& b, double& erra, double& errb){

  return 0;
}
*/



}// namespace libmmath
