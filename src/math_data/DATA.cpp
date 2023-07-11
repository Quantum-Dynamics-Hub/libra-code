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

#include "DATA.h"

//=========================== DATA class ================================

/// liblibra namespace
namespace liblibra{

/// libdata namespace
namespace libdata{

DATA::DATA(vector<double> d){

//  Data = d;
  if(Data.size()>0){ Data.clear(); }
  for(int i=0;i<d.size();i++){  Data.push_back(d[i]);  }

  is_ave = 0;
  is_var = 0;
  is_sd  = 0;
  is_se  = 0;
  is_mse = 0;
  is_rmse= 0;

  is_min_val = 0;
  is_min_indx = 0;
  is_max_val = 0;
  is_max_indx = 0;

  scale_factor = 1.0;  is_scale_factor = 1;
  shift_amount = 0.0;  is_shift_amount = 1;

}

DATA::DATA(int sz,double* d){

  if(Data.size()>0){   Data.clear(); }
  for(int i=0;i<sz;i++){  Data.push_back(d[i]); }

  is_ave = 0;
  is_var = 0;
  is_sd  = 0;
  is_se  = 0;
  is_mse = 0;
  is_rmse= 0;

  is_min_val = 0;
  is_min_indx = 0;
  is_max_val = 0;
  is_max_indx = 0;

  scale_factor = 1.0;  is_scale_factor = 1;
  shift_amount = 0.0;  is_shift_amount = 1;

}

DATA::DATA(boost::python::list obj){

   int sz = len(obj);
   double tmp;

   if(Data.size()>0){    Data.clear(); }

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

  is_min_val = 0;
  is_min_indx = 0;
  is_max_val = 0;
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

  is_min_val = 0;
  is_min_indx = 0;
  is_max_val = 0;
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

  if(d.is_min_val){ min_val = d.min_val; is_min_val = 1; }
  if(d.is_min_indx){min_indx = d.min_indx; is_min_indx = 1; }
  if(d.is_max_val){ max_val = d.max_val; is_max_val = 1; }
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

  is_min_val = 0;
  is_min_indx = 0;
  is_max_val = 0;
  is_max_indx = 0;

  scale_factor = 1.0;  is_scale_factor = 1;
  shift_amount = 0.0;  is_shift_amount = 1;

  if(Data.size()>0){ Data.clear(); }
  for(int i=0;i<d.Data.size();i++){  Data.push_back(d.Data[i]); }

  if(d.is_ave){ ave = d.ave; is_ave = 1; }
  if(d.is_var){ var = d.var; is_var = 1; }
  if(d.is_sd) { sd  = d.sd;  is_sd = 1;  }
  if(d.is_se) { se  = d.se;  is_se = 1;  }
  if(d.is_mse){ mse = d.mse; is_mse = 1; }
  if(d.is_rmse){rmse= d.rmse;is_rmse=1;  }

  if(d.is_min_val){ min_val = d.min_val; is_min_val = 1; }
  if(d.is_min_indx){min_indx = d.min_indx; is_min_indx = 1; }
  if(d.is_max_val){ max_val = d.max_val; is_max_val = 1; }
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
  int i;
  ave = 0.0;
  var = 0.0;
  mse = 0.0;
  mae = 0.0;

  double tmp;
  double sz = Data.size();

  if(sz>1){

    for(i=0;i<sz;i++){  ave += Data[i]; }    
    ave = ave/sz;

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

    Ave  = ave;         is_ave = 1; 
    Var  = var;         is_var = 1; 
    Sd   = sd;          is_sd  = 1; 
    Se   = se;          is_se  = 1; 
    Mse  = mse;         is_mse = 1; 
    Mae  = mae;         is_mae = 1; 
    Rmse = rmse;        is_rmse= 1; 

  }else{  std::cout<<"Warning: Data size is less than 2. Can not calculate statistics for such data\n";  }

  return 0;
}

int DATA::Calculate_Estimators(){
  double Ave, Var, Sd, Se, Mse, Mae, Rmse;

  return Calculate_Estimators(Ave, Var, Sd, Se, Mse, Mae, Rmse);
}


int DATA::Calculate_MiniMax(double& Min_Val,int& Min_Indx,double& Max_Val,int& Max_Indx){

  int sz = Data.size();


  min_val = max_val = Data[0];
  min_indx = max_indx = 0;

  for(int i=1;i<sz;i++){
      if(Data[i]<min_val) { min_val = Data[i]; min_indx = i; }
      if(Data[i]>max_val) { max_val = Data[i]; max_indx = i; }
  }

  Min_Val  = min_val;
  Min_Indx = min_indx;

  Max_Val  = max_val;
  Max_Indx = max_indx;

  is_min_val  = 1;
  is_max_val  = 1;
  is_min_indx = 1;
  is_max_indx = 1;

  return 0;
}


int DATA::Calculate_MiniMax(){
  double Min_Val,Max_Val;
  int Min_Indx, Max_Indx;

  return Calculate_MiniMax(Min_Val, Min_Indx, Max_Val, Max_Indx);  

}

int DATA::Calculate_Distribution(vector<double>& Interval,vector<double>& Density,vector<double>& Cumulative){
// Cumulative(x) - is the probability to find the point in the interval (Intrval[0], x)
  int i,j;
  int sz = Data.size();
  int szi= Interval.size(); // Number of sub-intervals
  double nrm = (1.0/double(sz));

  // Clear results
  if(Density.size()>0){ Density.clear();}
  if(Cumulative.size()>0){ Cumulative.clear();}

  // Init results
  for(i=0;i<szi;i++){
    Density.push_back(0.0);
    Cumulative.push_back(0.0);
  }

  // Count
  for(j=0;j<szi;j++){
    for(i=0;i<sz;i++){
      if(Data[i]<=Interval[j]){
        Cumulative[j] += nrm;
      }
    }// for i
  }// for j

  // Probability density
  Density[0] = 0.0;
  for(j=1;j<szi;j++){
    Density[j] = (Cumulative[j] - Cumulative[j-1])/(Interval[j] - Interval[j-1]);
  }  

  return 0;
}

boost::python::list DATA::Calculate_Distribution(boost::python::list Interval){

  int i;

  // Convert input list to vector
  int sz = boost::python::len(Interval);
  vector<double> interval(sz,0.0);
  for(i=0;i<sz;i++){ 
    interval[i] = boost::python::extract<double>(Interval[i]);
  }

  // Perform computations
  vector<double> density, cumulant;
  Calculate_Distribution(interval, density, cumulant);

  // Convert output vectors to lists
  boost::python::list Density;
  boost::python::list Cumulant;
  sz = density.size();

  for(i=0;i<sz;i++){
    Density.append(density[i]);
    Cumulant.append(cumulant[i]);
  }

  boost::python::list res;
  res.append(Density);
  res.append(Cumulant);

  return res;

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
  if(is_min_val && is_max_val && is_min_indx && is_max_indx){
     if(sc_fact>=0.0) { min_val = min_val*sc_fact + sh;
                        max_val = max_val*sc_fact + sh;
                      }
     else{              double tmp;
                        int itmp;
                        tmp = min_val;
                        itmp = min_indx;
                        min_val = sc_fact*max_val + sh;
                        max_val = sc_fact*tmp + sh;
                        min_indx = max_indx;
                        max_indx = itmp;
     }// else
  }// if


  return 0;
}



int DATA::invLinearTransformData(){
/**

   1. shift by  -shift_amount,   reset shift to 0.0
   2. scale by 1.0/scale_factor  reset scale to 1.0

*/

  int sz = Data.size();

  double sc_fact = scale_factor;
  double sh = shift_amount;

  scale_factor = 1.0;
  shift_amount = 0.0;

  for(int i=0;i<sz;i++){
      Data[i] = (Data[i] - sh )/sc_fact;
  }

  // Change corresponding estimators
  if(is_ave){  ave = (ave - sh)/sc_fact; }
  if(is_var){  var = var / (sc_fact * sc_fact);      }
  if(is_sd) {  sd  = sd / sc_fact;   }
  if(is_se) {  se  = se / sc_fact;   }
  if(is_mse){  mse = mse / (sc_fact * sc_fact); }
  if(is_rmse){ rmse= rmse / sc_fact; }

  if(is_min_val && is_max_val && is_min_indx && is_max_indx){
     if(sc_fact>=0.0) { min_val = (min_val - sh)/sc_fact;
                        max_val = (max_val - sh)/sc_fact;
                      }
     else{              double tmp;
                        int itmp;
                        tmp = min_val;
                        itmp = min_indx;
                        min_val = (max_val - sh) / sc_fact;
                        max_val = (tmp - sh) / sc_fact;
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

  if(is_min_val && is_max_val){
  }else{
   double a,b;
   int i,j;
   this->Calculate_MiniMax(a,i,b,j);
  }
  if(min_val==max_val){
    // Do nothing since in this case this operation is not defined
  }else{
    // Do rescaling and shifting
    double sc_factor = ((upp-low)/(max_val-min_val));
    double sh = (low - min_val*sc_factor);

    this->LinearTransformData(sc_factor,sh);

  }

  return 0;
}

int DATA::ShiftData(double sh){

  this->LinearTransformData(1.0,sh);

  return 0;
}
int DATA::NormalizeData(){
  /*
  double Ave,Var,Sd,Se,Mse,Mae,Rmse;
  if(is_ave&&is_sd){
  }else{
   this->Calculate_Estimators(Ave,Var,Sd,Se,Mse,Mae,Rmse);
  }

  if(Sd!=0.0){
  this->LinearTransformData((1.0/Sd),(-Ave/Sd));
  }
  */
  this->Calculate_Estimators();
  this->LinearTransformData( 1.0/sd, (-ave/sd) );


  return 0;
}


/*
int DATA::Lin_Regression(int flag,double& a,double& b, double& erra, double& errb){

  return 0;
}
*/


}// namespace libdata
}// namespace liblibra
