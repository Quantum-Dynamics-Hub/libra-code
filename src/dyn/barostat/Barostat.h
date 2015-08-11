#ifndef BAROSTAT_H
#define BAROSTAT_H

#include "../../mmath/libmmath.h"

namespace libdyn{
namespace libbarostat{

class Barostat{

  //--------- Auxiliary internal functions -------------
  void init_variables();// Initializes variables
  void copy_content(const Barostat&); // Copies the content which is defined
  void extract_dictionary(boost::python::dict);

  double Nf_t;                  int is_Nf_t;            // Number of translational degrees of freedom in system
                                                        // to which the barostat is applied
  double Nf_r;                  int is_Nf_r;            // Number of rotational degrees of freedom in system
                                                        // to which the barostat is applied
  double Nf_b;                  int is_Nf_b;            // Number of barostat degrees of freedom in system
                                                        // to which the barostat is applied
                                                        // = 1 for isotropic dilation, = 9 for flexible cell


public:
  // Dynamic variables and internal parameters for Anderson barostat
  MATRIX3x3 ksi_eps;            int is_ksi_eps;         // Barostat "velocity" tensor
  MATRIX3x3 G_eps;              int is_G_eps;           // Barostat force tensor
  double eps_iso;               int is_eps_iso;         // Isotropic barostat position
  double ksi_eps_iso;           int is_ksi_eps_iso;     // Isotropic barostat momentum
  double G_eps_iso;             int is_G_eps_iso;       // Isotropic barostat force
  double Wg;                    int is_Wg;              // Barostat mass
//  double Nf_t;                  int is_Nf_t;            // Number of translational degrees of freedom in system
                                                        // to which the barostat is applied
//  double Nf_r;                  int is_Nf_r;            // Number of rotational degrees of freedom in system
                                                        // to which the barostat is applied
//  double Nf_b;                  int is_Nf_b;            // Number of barostat degrees of freedom in system
                                                        // to which the barostat is applied
                                                        // = 1 for isotropic dilation, = 9 for flexible cell

  // Input parametes
  double nu_baro;               int is_nu_baro;         // Barostat frequency parameter
  double Pressure;              int is_Pressure;        // Target pressure of barostat
  std::string barostat_type;    int is_barostat_type;   // Type of barostat in use: Andersen, Parrinello-Rahman, Berendsen


  //----------- Basic class operations ---------------------------
  // Defined in Barostat.cpp
  Barostat();                   // constructor
  Barostat(boost::python::dict);
  Barostat(const Barostat&);    // copy-constructor
 ~Barostat();                   // destructor

  Barostat& operator=(const Barostat&); // assignment operator
  void show_info();
  void set(object);
  void save(boost::property_tree::ptree& pt,std::string path);
  void load(boost::property_tree::ptree& pt,std::string path,int& status);


  //-------------- Getters, setters -----------------------------
  // Defined in Barostat_aux.cpp
  void set_Nf_t(double nf_t){ Nf_t = nf_t; is_Nf_t = 1; }
  void set_Nf_r(double nf_r){ Nf_r = nf_r; is_Nf_r = 1; }
  void set_Nf_b(double nf_b){ Nf_b = nf_b; is_Nf_b = 1; }
  void set_Nf_t(int nf_t){ Nf_t = nf_t; is_Nf_t = 1; }
  void set_Nf_r(int nf_r){ Nf_r = nf_r; is_Nf_r = 1; }
  void set_Nf_b(int nf_b){ Nf_b = nf_b; is_Nf_b = 1; }

  double get_Nf_t(){ if(is_Nf_t){ return Nf_t; } else{ std::cout<<"Error: Nf_t is not defined\n"; exit(1); } }
  double get_Nf_r(){ if(is_Nf_r){ return Nf_r; } else{ std::cout<<"Error: Nf_r is not defined\n"; exit(1); } }
  double get_Nf_b(){ if(is_Nf_b){ return Nf_b; } else{ std::cout<<"Error: Nf_b is not defined. Verify that either NPT or NPT_FLEX ensemble is used\n"; exit(1); } }

  double ekin_baro();
  void apply_barostat_force(double);
  void scale_velocity(double);
  void propagate_velocity(double dt,double ksi_b); // exp(iL*dt) = (G_eps/Wg - ksi_b * ksi_eps)d/dksi_eps

  //------------------- Methods ---------------------------------
  // Defined in Barostat_methods.cpp
  void update_barostat_forces(double,double,double, double);
  void update_barostat_forces(double,double,double, MATRIX3x3&);
  void init(double);
  MATRIX3x3 pos_scale(double);
  MATRIX3x3 vpos_scale(double);
  MATRIX3x3 vel_scale(double,double);
  MATRIX3x3 ang_vel_scale(double,double);
  void cool();


};

void save(boost::property_tree::ptree& pt,std::string path,vector<Barostat>& vt);
void load(boost::property_tree::ptree& pt,std::string path,vector<Barostat>& vt,int& status);



}// namespace libbarostat
}// namespace libdyn



#endif // BAROSTAT_H
