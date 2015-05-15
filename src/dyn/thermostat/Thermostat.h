#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include "../../mmath/libmmath.h"

namespace libdyn{
namespace libthermostat{


class Thermostat{

  //--------- Auxiliary internal functions -------------
  void init_variables();// Initializes variables
  void copy_content(const Thermostat&); // Copies the content which is defined
  void extract_dictionary(boost::python::dict);


  // Dynamic variables and internal parameters for Nose-Hoover chain thermostat
  vector<double> s_t;           int s_t_size;           // Nose-Hoover chain thermostat positions coupled to tr. deg. of freedom
  vector<double> s_r;           int s_r_size;           // Nose-Hoover chain thermostat positions coupled to rot. deg. of freedom
  vector<double> s_b;           int s_b_size;           // Nose-Hoover chain thermostat positions coupled to barostat
  vector<double> ksi_t;         int ksi_t_size;         // Nose-Hoover chain thermostat velocities coupled to tr. deg. of freedom
  vector<double> ksi_r;         int ksi_r_size;         // Nose-Hoover chain thermostat velocities coupled to rot. deg. of freedom
  vector<double> ksi_b;         int ksi_b_size;         // Nose-Hoover chain thermostat velocities coupled to barostat
  vector<double> G_t;           int G_t_size;           // Nose-Hoover chain thermostat forces coupled to tr. deg. of freedom
  vector<double> G_r;           int G_r_size;           // Nose-Hoover chain thermostat forces coupled to rot. deg. of freedom
  vector<double> G_b;           int G_b_size;           // Nose-Hoover chain thermostat forces coupled to barostat
  vector<double> Q_t;           int Q_t_size;           // Nose-Hoover chain thermostat masses coupled to tr. deg. of freedom
  vector<double> Q_r;           int Q_r_size;           // Nose-Hoover chain thermostat masses coupled rot. deg. of freedom
  vector<double> Q_b;           int Q_b_size;           // Nose-Hoover chain thermostat masses coupled to barostat

  double Nf_t;                  int is_Nf_t;            // Number of translational degrees of freedom in system
                                                        // to which the thermostat is applied
  double Nf_r;                  int is_Nf_r;            // Number of rotational degrees of freedom in system
                                                        // to which the thermostat is applied
  double Nf_b;                  int is_Nf_b;            // Number of barostat degrees of freedom in system
                                                        // to which the thermostat is applied
                                                        // = 1 for isotropic dilation, = 9 for flexible cell


public:  
  // Dynamic variables and internal parameters for Nose-Poincare thermostat
  double s_var;                 int is_s_var;           // Thermostat variable
  double Ps;                    int is_Ps;              // Momentum conjugate to s_var
  double Q;                     int is_Q;               // Thermostat mass
//  double Nf_t;                  int is_Nf_t;            // Number of translational degrees of freedom in system
                                                        // to which the thermostat is applied
//  double Nf_r;                  int is_Nf_r;            // Number of rotational degrees of freedom in system
                                                        // to which the thermostat is applied
//  double Nf_b;                  int is_Nf_b;            // Number of barostat degrees of freedom in system
                                                        // to which the thermostat is applied
                                                        // = 1 for isotropic dilation, = 9 for flexible cell

  // Input parameters
  double NHC_size;              int is_NHC_size;        // Length (size) of the Nose-Hoover thermostat chain
  double nu_therm;              int is_nu_therm;        // Thermostat frequency parameter (time-scale)
  double Temperature;           int is_Temperature;     // Target temperature of the thermostat
  std::string thermostat_type;  int is_thermostat_type; // Type of the thermostat in use: Nose-Poincare, Nose-Hoover, etc.



  //----------- Basic class operations ---------------------------
  // Defined in Thermostat.cpp
  Thermostat();                   // constructor
  Thermostat(boost::python::dict);
  Thermostat(const Thermostat&);  // copy-constructor
 ~Thermostat();                   // destructor

  Thermostat& operator=(const Thermostat&); // assignment operator
  void show_info();
  void set(object);
  void save(boost::property_tree::ptree& pt,std::string path);
  void load(boost::property_tree::ptree& pt,std::string path,int& status);


  //-------------- Getters, setters -----------------------------
  // Defined in Thermostat_aux.cpp
  double energy();
  void set_Nf_t(double nf_t){ Nf_t = nf_t; is_Nf_t = 1; }
  void set_Nf_r(double nf_r){ Nf_r = nf_r; is_Nf_r = 1; }
  void set_Nf_b(double nf_b){ Nf_b = nf_b; is_Nf_b = 1; }
  void set_Nf_t(int nf_t){ Nf_t = nf_t; is_Nf_t = 1; }
  void set_Nf_r(int nf_r){ Nf_r = nf_r; is_Nf_r = 1; }
  void set_Nf_b(int nf_b){ Nf_b = nf_b; is_Nf_b = 1; }

  double get_Nf_t(){ if(is_Nf_t){ return Nf_t; } else{ std::cout<<"Error: Nf_t is not defined\n"; exit(1); } }
  double get_Nf_r(){ if(is_Nf_r){ return Nf_r; } else{ std::cout<<"Error: Nf_r is not defined\n"; exit(1); } }
  double get_Nf_b(){ if(is_Nf_b){ return Nf_b; } else{ std::cout<<"Error: Nf_b is not defined\n"; exit(1); } }

  double get_s_var(){ return s_var; }
  double get_ksi_t(){ if(ksi_t_size>0) {return ksi_t[0];} else{ return 0.0; /*std::cout<<"Error: ksi_t vector is of zero size\n"; exit(1);*/ } }
  double get_ksi_r(){ if(ksi_r_size>0) {return ksi_r[0];} else{ return 0.0; /*std::cout<<"Error: ksi_r vector is of zero size\n"; exit(1);*/ } }
  double get_ksi_b(){ if(ksi_b_size>0) {return ksi_b[0];} else{ return 0.0; /*std::cout<<"Error: ksi_b vector is of zero size\n"; exit(1);*/ } }

  //------------------- Methods ---------------------------------
  // Defined in Thermostat_methods.cpp
  void propagate_sPs(double);
  void propagate_Ps(double);
  double vel_scale(double);
  double ang_vel_scale(double);
  void update_thermostat_forces(double, double, double);
  void update_thermostat_forces(double, double, double, int);
  void init_nhc();
  void propagate_nhc(double,double, double, double);


};

void save(boost::property_tree::ptree& pt,std::string path,vector<Thermostat>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Thermostat>& vt,int& status);


}// namespace libthermostat
}// namespace libdyn



#endif // THERMOSTAT_H
